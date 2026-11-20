#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/rdkit_extensions/atomistic_conversions.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/helm/monomer_coordgen.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/rdkit_extensions/molops.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <set>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/json.hpp>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/Depictor/DepictUtils.h>
#include <rdkit/GraphMol/Depictor/RDDepictor.h>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/ROMol.h>

namespace schrodinger
{
namespace rdkit_extensions
{

namespace
{

using ResidueKey = std::tuple<std::string, int, std::string>;
using MonomerKey = std::pair<std::string, int>;
using ResidueToMonomerMap = std::map<ResidueKey, std::tuple<std::string, int>>;
constexpr double BACKBONE_BOND_SCALE = 1.0;
constexpr double SIDECHAIN_BOND_SCALE = 0.7;

struct ResidueGeometry {
    std::vector<unsigned int> atom_ids;
    std::vector<RDGeom::Point3D> positions;
    RDGeom::Point3D center;
    RDGeom::Point3D size;
    std::optional<unsigned int> monomer_idx;
};

using ResidueMap = std::map<ResidueKey, ResidueGeometry>;
using AtomResidueMap = std::unordered_map<unsigned int, ResidueKey>;
using ResidueKeys = std::set<ResidueKey>;
using DisulfidePartners = std::map<ResidueKey, ResidueKey>;
using LinearCounts = std::map<ResidueKey, size_t>;

bool has_atom_name(const RDKit::Atom* atom, std::string_view name)
{
    const auto* info =
        dynamic_cast<const RDKit::AtomPDBResidueInfo*>(atom->getMonomerInfo());
    if (!info)
        return false;
    auto atom_name = info->getName();
    atom_name.erase(std::remove(atom_name.begin(), atom_name.end(), ' '),
                    atom_name.end());
    return atom_name == name;
}

std::optional<size_t> atom_offset(const ResidueGeometry& residue,
                                  unsigned int atom_id)
{
    const auto found = std::ranges::lower_bound(residue.atom_ids, atom_id);
    if (found == residue.atom_ids.end() || *found != atom_id)
        return std::nullopt;
    return found - residue.atom_ids.begin();
}

std::optional<unsigned int> named_atom(const RDKit::ROMol& mol,
                                       const ResidueGeometry& residue,
                                       std::string_view name)
{
    for (auto atom_id : residue.atom_ids) {
        if (has_atom_name(mol.getAtomWithIdx(atom_id), name))
            return atom_id;
    }
    return std::nullopt;
}

struct PeptideLink {
    ResidueKey next;
    RDGeom::Point3D axis;
};

using PeptideLinks = std::map<ResidueKey, PeptideLink>;

// These are 2D drawings, not cis/trans configurations of the peptide bond.
// C(prev) and N(next) belong to the adjacent residues; R is this residue's
// side chain. All sketches enter along the same downward C(prev)-N bond;
// the whole drawing can be mirrored vertically.
//
// ALL_TRANS: the backbone stays in an extended zigzag with 120 degrees
// at CA; the side chain points on the opposite side from the previous residue.
//           R
//           |
// C(prev)  CA       N(next)
//       \  /  \     /
//        N      C
//               ||
//               O
//
// LINEAR_CIS_N_CA: N-CA is horizontal; 150 degrees at CA, so the side
// chain stays on the same side as in the preceding residue.
//                O
//                ||
// C(prev)       C
//        \      / \
//         N----CA   N(next)
//              |
//              R
//
// CIS_N_CA: reflect after N-CA. With the incoming bond drawn as above, the
// backbone bends upward and this residue's R points below it; retain 120
// degrees at CA.
//           O       N(next)
//            \\     /
//             C
//             |
// C(prev)     CA
//       \    /  \
//        N       R
//
// CIS_CA_CO: reflect after CA-C. With the same incoming bond, the backbone
// bends downward and this residue's R points above it; retain 120 degrees
// at CA.
//           R         O
//           |        //
// C(prev)   CA      /
//       \   /  \    /
//        N      C
//               |
//             N(next)
enum class ResidueShape {
    ALL_TRANS,
    LINEAR_CIS_N_CA,
    CIS_N_CA,
    CIS_CA_CO,
};

ResidueShape initial_residue_shape(const RDKit::ROMol& mol,
                                   const ResidueGeometry& residue)
{
    // Read the residue name from its N atom. Start proline and D-proline in
    // LINEAR_CIS_N_CA; all other residues start in ALL_TRANS.
    const auto nitrogen = named_atom(mol, residue, "N");
    const auto* info =
        nitrogen ? dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                       mol.getAtomWithIdx(*nitrogen)->getMonomerInfo())
                 : nullptr;
    const auto proline = info && (info->getResidueName() == "PRO" ||
                                  info->getResidueName() == "DPR");
    return proline ? ResidueShape::LINEAR_CIS_N_CA : ResidueShape::ALL_TRANS;
}

RDGeom::Point3D rotate_around(const RDGeom::Point3D& point,
                              const RDGeom::Point3D& center, double angle)
{
    const auto p = point - center;
    const auto cosine = std::cos(angle), sine = std::sin(angle);
    return center + RDGeom::Point3D(p.x * cosine - p.y * sine,
                                    p.x * sine + p.y * cosine, 0);
}

RDGeom::Point3D reflect_across_bond(const RDGeom::Point3D& point,
                                    const RDGeom::Point3D& start,
                                    const RDGeom::Point3D& end)
{
    const auto bond = end - start;
    const auto length = std::hypot(bond.x, bond.y);
    if (length == 0)
        return point;
    const auto axis = bond / length;
    const auto p = point - start;
    return start + axis * (2 * (p.x * axis.x + p.y * axis.y)) - p;
}

void orient_peptide_sidechain(const RDKit::ROMol& mol,
                              const ResidueGeometry& residue,
                              unsigned int nitrogen, unsigned int alpha,
                              unsigned int carbonyl,
                              const RDGeom::Point3D& outward,
                              RDKit::Conformer& output)
{
    // Find the side-chain atoms attached to CA. For a ring that closes at N
    // (proline), reflect them across N-CA if that points the ring outward.
    // Otherwise, rotate them toward outward and scale their bond lengths.
    const auto nitrogen_position = output.getAtomPos(nitrogen);
    const auto alpha_position = output.getAtomPos(alpha);
    const auto in_residue = [&](unsigned int atom_id) {
        return std::ranges::binary_search(residue.atom_ids, atom_id);
    };
    std::vector<unsigned int> pending;
    std::optional<unsigned int> sidechain_anchor =
        named_atom(mol, residue, "CB");
    for (const auto* neighbor : mol.atomNeighbors(mol.getAtomWithIdx(alpha))) {
        const auto atom_id = neighbor->getIdx();
        if (neighbor->getAtomicNum() == 1) {
            continue;
        }
        if (atom_id == nitrogen || atom_id == carbonyl ||
            !in_residue(atom_id)) {
            continue;
        }
        pending.push_back(atom_id);
        if (!sidechain_anchor) {
            sidechain_anchor = atom_id;
        }
    }
    std::set<unsigned int> sidechain_atoms;
    bool closes_to_nitrogen = false;
    while (!pending.empty()) {
        const auto atom_id = pending.back();
        pending.pop_back();
        if (!sidechain_atoms.insert(atom_id).second) {
            continue;
        }
        for (const auto* neighbor :
             mol.atomNeighbors(mol.getAtomWithIdx(atom_id))) {
            const auto neighbor_id = neighbor->getIdx();
            if (neighbor_id == nitrogen) {
                closes_to_nitrogen = true;
            } else if (neighbor_id != alpha && neighbor_id != carbonyl &&
                       in_residue(neighbor_id)) {
                pending.push_back(neighbor_id);
            }
        }
    }
    if (sidechain_anchor && sidechain_atoms.contains(*sidechain_anchor)) {
        const auto anchor = output.getAtomPos(*sidechain_anchor);
        if (closes_to_nitrogen) {
            const auto reflected =
                reflect_across_bond(anchor, nitrogen_position, alpha_position);
            const auto original_score =
                (anchor.x - alpha_position.x) * outward.x +
                (anchor.y - alpha_position.y) * outward.y;
            const auto reflected_score =
                (reflected.x - alpha_position.x) * outward.x +
                (reflected.y - alpha_position.y) * outward.y;
            if (reflected_score > original_score) {
                for (auto atom_id : sidechain_atoms) {
                    output.setAtomPos(
                        atom_id,
                        reflect_across_bond(output.getAtomPos(atom_id),
                                            nitrogen_position, alpha_position));
                }
            }
        } else {
            const auto direction = anchor - alpha_position;
            const auto sidechain_angle = std::atan2(outward.y, outward.x) -
                                         std::atan2(direction.y, direction.x);
            for (auto atom_id : sidechain_atoms) {
                const auto rotated =
                    rotate_around(output.getAtomPos(atom_id), alpha_position,
                                  sidechain_angle);
                output.setAtomPos(atom_id,
                                  alpha_position + (rotated - alpha_position) *
                                                       (SIDECHAIN_BOND_SCALE /
                                                        BACKBONE_BOND_SCALE));
            }
        }
    }
}

bool place_residue_by_bond(const ResidueGeometry& residue,
                           unsigned int source_start_atom,
                           unsigned int source_end_atom,
                           const RDGeom::Point3D& target_start,
                           const RDGeom::Point3D& target_end,
                           RDKit::Conformer& output)
{
    // Match one bond in RDKit's residue drawing to its target bond in the
    // backbone layout. Apply the same translation, rotation, and uniform
    // scale to every atom, bringing the side chain along with the backbone.
    const auto source_start =
        residue.positions[*atom_offset(residue, source_start_atom)];
    const auto source_end =
        residue.positions[*atom_offset(residue, source_end_atom)];
    const auto source_bond = source_end - source_start;
    const auto target_bond = target_end - target_start;
    const auto source_length = std::hypot(source_bond.x, source_bond.y);
    if (source_length == 0) {
        return false;
    }
    const auto scale = std::hypot(target_bond.x, target_bond.y) / source_length;
    const auto angle = std::atan2(target_bond.y, target_bond.x) -
                       std::atan2(source_bond.y, source_bond.x);
    const auto cosine = std::cos(angle), sine = std::sin(angle);
    for (size_t i = 0; i < residue.atom_ids.size(); ++i) {
        const auto p = (residue.positions[i] - source_start) * scale;
        output.setAtomPos(residue.atom_ids[i],
                          target_start +
                              RDGeom::Point3D(p.x * cosine - p.y * sine,
                                              p.x * sine + p.y * cosine, 0));
    }
    return true;
}

void place_terminal_amide_nitrogen(const RDKit::ROMol& mol,
                                   const ResidueGeometry& residue,
                                   unsigned int backbone_nitrogen,
                                   unsigned int carbonyl,
                                   const RDGeom::Point3D& alpha_position,
                                   const RDGeom::Point3D& carbonyl_position,
                                   const RDGeom::Point3D& oxygen_position,
                                   RDKit::Conformer& output)
{
    // After placing C and O in the backbone layout, put a terminal C(=O)NH2
    // nitrogen in the remaining trigonal direction and move its H atoms with
    // it.
    const RDKit::Atom* amide = nullptr;
    for (const auto* neighbor :
         mol.atomNeighbors(mol.getAtomWithIdx(carbonyl))) {
        const auto id = neighbor->getIdx();
        if (id == backbone_nitrogen || neighbor->getAtomicNum() != 7 ||
            !std::ranges::binary_search(residue.atom_ids, id) ||
            mol.getBondBetweenAtoms(carbonyl, id)->getBondType() !=
                RDKit::Bond::SINGLE) {
            continue;
        }
        const auto terminal = std::ranges::all_of(
            mol.atomNeighbors(neighbor), [&](const auto* attached) {
                return attached->getIdx() == carbonyl ||
                       attached->getAtomicNum() == 1;
            });
        if (terminal) {
            amide = neighbor;
            break;
        }
    }
    if (!amide) {
        return;
    }

    const auto ca = alpha_position - carbonyl_position;
    const auto oxygen = oxygen_position - carbonyl_position;
    const auto ca_length = std::hypot(ca.x, ca.y);
    const auto oxygen_length = std::hypot(oxygen.x, oxygen.y);
    if (ca_length == 0 || oxygen_length == 0) {
        return;
    }
    const auto opposite = ca / ca_length + oxygen / oxygen_length;
    const auto opposite_length = std::hypot(opposite.x, opposite.y);
    if (opposite_length == 0) {
        return;
    }

    const auto target =
        carbonyl_position - opposite * (RDDepict::BOND_LEN *
                                        SIDECHAIN_BOND_SCALE / opposite_length);
    const auto displacement = target - output.getAtomPos(amide->getIdx());
    output.setAtomPos(amide->getIdx(), target);
    for (const auto* attached : mol.atomNeighbors(amide)) {
        if (attached->getAtomicNum() == 1) {
            output.setAtomPos(attached->getIdx(),
                              output.getAtomPos(attached->getIdx()) +
                                  displacement);
        }
    }
}

void place_alpha_hydrogen(const RDKit::ROMol& mol,
                          const ResidueGeometry& residue, unsigned int alpha,
                          const RDGeom::Point3D& nitrogen_position,
                          const RDGeom::Point3D& alpha_position,
                          const RDGeom::Point3D& carbonyl_position,
                          RDKit::Conformer& output)
{
    // Moving C and orienting the side chain can leave a single explicit H on
    // CA in the wrong position. Put it perpendicular to the N-C backbone,
    // opposite CB; leave residues with two H atoms to their RDKit depiction.
    std::vector<unsigned int> hydrogens;
    for (const auto* neighbor : mol.atomNeighbors(mol.getAtomWithIdx(alpha))) {
        if (neighbor->getAtomicNum() == 1) {
            hydrogens.push_back(neighbor->getIdx());
        }
    }
    if (hydrogens.size() != 1) {
        return;
    }

    const auto backbone = carbonyl_position - nitrogen_position;
    const auto backbone_length = std::hypot(backbone.x, backbone.y);
    if (backbone_length == 0) {
        return;
    }
    auto direction = RDGeom::Point3D(-backbone.y / backbone_length,
                                     backbone.x / backbone_length, 0);
    if (const auto beta = named_atom(mol, residue, "CB")) {
        const auto beta_direction = output.getAtomPos(*beta) - alpha_position;
        if (direction.x * beta_direction.x + direction.y * beta_direction.y >
            0) {
            direction *= -1;
        }
    }
    output.setAtomPos(hydrogens.front(),
                      alpha_position + direction * (RDDepict::BOND_LEN *
                                                    SIDECHAIN_BOND_SCALE));
}

bool place_peptide_backbone_residue(const RDKit::ROMol& mol,
                                    const ResidueGeometry& residue,
                                    const RDGeom::Point3D& target_nitrogen,
                                    const RDGeom::Point3D& target_alpha,
                                    const RDGeom::Point3D& target_carbonyl,
                                    const RDGeom::Point3D& target_oxygen,
                                    const RDGeom::Point3D& outward,
                                    RDKit::Conformer& output)
{
    // Place a residue with an N-CA-C backbone, regardless of its side chain.
    // Align N-CA with the chosen layout, orient the side chain outward, then
    // place the carbonyl and terminal groups.
    const auto nitrogen = named_atom(mol, residue, "N");
    const auto alpha = named_atom(mol, residue, "CA");
    const auto carbonyl = named_atom(mol, residue, "C");
    if (!nitrogen || !alpha || !carbonyl ||
        !mol.getBondBetweenAtoms(*nitrogen, *alpha) ||
        !mol.getBondBetweenAtoms(*alpha, *carbonyl)) {
        return false;
    }

    if (!place_residue_by_bond(residue, *nitrogen, *alpha, target_nitrogen,
                               target_alpha, output)) {
        return false;
    }

    orient_peptide_sidechain(mol, residue, *nitrogen, *alpha, *carbonyl,
                             outward, output);

    output.setAtomPos(*carbonyl, target_carbonyl);
    if (const auto oxygen = named_atom(mol, residue, "O")) {
        output.setAtomPos(*oxygen, target_oxygen);
    }
    place_terminal_amide_nitrogen(mol, residue, *nitrogen, *carbonyl,
                                  target_alpha, target_carbonyl, target_oxygen,
                                  output);
    if (const auto oxygen_terminal = named_atom(mol, residue, "OXT")) {
        output.setAtomPos(*oxygen_terminal,
                          target_carbonyl +
                              (target_alpha - target_nitrogen) *
                                  (SIDECHAIN_BOND_SCALE / BACKBONE_BOND_SCALE));
    }

    return true;
}

bool place_terminal_amine(const RDKit::ROMol& mol,
                          const ResidueGeometry& residue,
                          const RDGeom::Point3D& target_nitrogen,
                          RDKit::Conformer& output)
{
    // Place an [am] endcap at the next peptide N position. It has no CA or
    // carbonyl C, so move its N and attached H atoms together from RDKit's
    // residue drawing.
    const auto nitrogen = named_atom(mol, residue, "N");
    if (!nitrogen || named_atom(mol, residue, "CA") ||
        named_atom(mol, residue, "C")) {
        return false;
    }
    if (std::ranges::any_of(residue.atom_ids, [&](auto atom_id) {
            return atom_id != *nitrogen &&
                   mol.getAtomWithIdx(atom_id)->getAtomicNum() != 1;
        })) {
        return false;
    }
    const auto source_nitrogen =
        residue.positions[*atom_offset(residue, *nitrogen)];
    for (size_t i = 0; i < residue.atom_ids.size(); ++i) {
        output.setAtomPos(residue.atom_ids[i],
                          target_nitrogen +
                              (residue.positions[i] - source_nitrogen) *
                                  SIDECHAIN_BOND_SCALE);
    }
    return true;
}

bool place_n_terminal_carbonyl_cap(const RDKit::ROMol& mol,
                                   const ResidueGeometry& residue,
                                   const RDGeom::Point3D& target_carbonyl,
                                   const RDGeom::Point3D& target_oxygen,
                                   RDKit::Conformer& output)
{
    // [ac] has a carbonyl C but no N-CA pair. Anchor its C=O at the
    // position preceding the first peptide nitrogen.
    const auto carbonyl = named_atom(mol, residue, "C");
    const auto oxygen = named_atom(mol, residue, "O");
    if (!carbonyl || !oxygen || named_atom(mol, residue, "N") ||
        named_atom(mol, residue, "CA")) {
        return false;
    }
    return place_residue_by_bond(residue, *carbonyl, *oxygen, target_carbonyl,
                                 target_oxygen, output);
}

void clear_molfile_bond_stereo(RDKit::ROMol& mol)
{
    // Remove Molfile bond stereo properties tied to the old coordinates
    // before replacing the conformer.
    for (auto* bond : mol.bonds()) {
        if (bond->hasProp(RDKit::common_properties::_MolFileBondStereo)) {
            bond->clearProp(RDKit::common_properties::_MolFileBondStereo);
        }
        if (bond->hasProp(RDKit::common_properties::_MolFileBondCfg)) {
            bond->clearProp(RDKit::common_properties::_MolFileBondCfg);
        }
    }
}

bool has_complete_residue_info(const RDKit::ROMol& mol)
{
    // Require PDB residue info on every heavy atom. An H may lack its own
    // info if it is bonded to an atom that has it; require at least one tag.
    bool found_residue = false;
    for (const auto* atom : mol.atoms()) {
        if (dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                atom->getMonomerInfo())) {
            found_residue = true;
            continue;
        }
        if (atom->getAtomicNum() != 1) {
            return false;
        }
        bool has_residue_neighbor = false;
        for (const auto* neighbor : mol.atomNeighbors(atom)) {
            if (dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                    neighbor->getMonomerInfo())) {
                has_residue_neighbor = true;
                break;
            }
        }
        if (!has_residue_neighbor) {
            return false;
        }
    }
    return found_residue;
}

ResidueMap get_residue_geometry(const RDKit::ROMol& mol)
{
    // Group atoms by PDB residue, including H atoms without their own tag.
    // Draw each residue as an RDKit fragment and record its 2D coordinates,
    // center, and size for the atomistic layout.
    ResidueMap residues;
    for (const auto* atom : mol.atoms()) {
        const auto* info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());
        if (!info) {
            if (atom->getAtomicNum() == 1) {
                continue;
            }
            throw std::invalid_argument(
                "Monomer layout requires PDB residue info on every heavy atom");
        }
        residues[{info->getChainId(), info->getResidueNumber(),
                  info->getInsertionCode()}]
            .atom_ids.push_back(atom->getIdx());
    }

    // PDB hydrogens may lack residue info. Keep them with their bonded residue.
    for (const auto* atom : mol.atoms()) {
        if (atom->getAtomicNum() != 1 ||
            dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                atom->getMonomerInfo())) {
            continue;
        }
        bool assigned = false;
        for (const auto* neighbor : mol.atomNeighbors(atom)) {
            const auto* info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                neighbor->getMonomerInfo());
            if (info) {
                residues[{info->getChainId(), info->getResidueNumber(),
                          info->getInsertionCode()}]
                    .atom_ids.push_back(atom->getIdx());
                assigned = true;
                break;
            }
        }
        if (!assigned) {
            throw std::invalid_argument(
                "Monomer layout requires residue info for isolated hydrogens");
        }
    }

    for (auto& [key, residue] : residues) {
        std::ranges::sort(residue.atom_ids);
        auto fragment = ExtractMolFragment(mol, residue.atom_ids, false);
        compute2DCoords(*fragment, {}, false);
        const auto& conf = fragment->getConformer();
        auto min_x = conf.getAtomPos(0).x, max_x = min_x;
        auto min_y = conf.getAtomPos(0).y, max_y = min_y;
        RDGeom::Point3D center(0, 0, 0);
        for (unsigned int i = 0; i < fragment->getNumAtoms(); ++i) {
            const auto p = conf.getAtomPos(i);
            residue.positions.push_back(p);
            center += p;
            min_x = std::min(min_x, p.x);
            max_x = std::max(max_x, p.x);
            min_y = std::min(min_y, p.y);
            max_y = std::max(max_y, p.y);
        }
        const auto nitrogen = named_atom(mol, residue, "N");
        const auto carbonyl = named_atom(mol, residue, "C");
        const auto alpha = named_atom(mol, residue, "CA");
        residue.center =
            nitrogen && carbonyl && alpha
                ? (residue.positions[*atom_offset(residue, *nitrogen)] +
                   residue.positions[*atom_offset(residue, *carbonyl)]) /
                      2
                : center / fragment->getNumAtoms();
        residue.size = {max_x - min_x + RDDepict::BOND_LEN,
                        max_y - min_y + RDDepict::BOND_LEN, 0};
    }
    return residues;
}

double prepare_monomer_layout(RDKit::ROMol& monomer_mol, ResidueMap& residues)
{
    // Match atomistic residues to monomer beads, generate the bead layout,
    // and return a uniform scale that gives their atomistic drawings room.
    const auto mapping_text =
        monomer_mol.getProp<std::string>(ST_TO_MOL_MONOMERIC_MAPPING);
    const auto residue_to_monomer = boost::json::value_to<ResidueToMonomerMap>(
        boost::json::parse(mapping_text));
    std::map<MonomerKey, unsigned int> monomer_indices;
    for (const auto* atom : monomer_mol.atoms()) {
        monomer_indices[{get_polymer_id(atom),
                         static_cast<int>(get_residue_number(atom))}] =
            atom->getIdx();
    }
    for (auto& [key, residue] : residues) {
        const auto mapped = residue_to_monomer.find(key);
        if (mapped == residue_to_monomer.end()) {
            continue;
        }
        const auto& [polymer_id, residue_number] = mapped->second;
        const auto found = monomer_indices.find({polymer_id, residue_number});
        if (found == monomer_indices.end()) {
            continue;
        }
        residue.monomer_idx = found->second;
        monomer_mol.getAtomWithIdx(*residue.monomer_idx)
            ->setProp(MONOMER_ITEM_SIZE,
                      std::array<double, 3>{residue.size.x, residue.size.y, 0});
    }
    compute_monomer_mol_coords(monomer_mol);
    const auto& monomer_conf = monomer_mol.getConformer();

    // Some monomer layout paths use fixed center spacing. Enlarge that
    // spacing uniformly so expanded atomistic residues cannot overlap.
    double layout_scale = 1;
    for (auto first = residues.begin(); first != residues.end(); ++first) {
        if (!first->second.monomer_idx) {
            continue;
        }
        for (auto second = residues.begin(); second != first; ++second) {
            if (!second->second.monomer_idx) {
                continue;
            }
            const auto delta =
                monomer_conf.getAtomPos(*first->second.monomer_idx) -
                monomer_conf.getAtomPos(*second->second.monomer_idx);
            const auto distance = std::hypot(delta.x, delta.y);
            if (distance == 0) {
                continue;
            }
            const auto radius = [](const auto& residue) {
                return std::hypot(residue.size.x, residue.size.y) / 2;
            };
            const auto required_distance =
                radius(first->second) + radius(second->second);
            layout_scale = std::max(layout_scale, required_distance / distance);
        }
    }

    return layout_scale;
}

AtomResidueMap map_atoms_to_residues(const ResidueMap& residues)
{
    // Build a reverse lookup from each atom index to its PDB residue.
    AtomResidueMap atom_residue;
    for (const auto& [key, residue] : residues) {
        for (auto idx : residue.atom_ids) {
            atom_residue.emplace(idx, key);
        }
    }

    return atom_residue;
}

DisulfidePartners find_disulfide_partners(const RDKit::ROMol& mol,
                                          const AtomResidueMap& atom_residue)
{
    // Find S-S bonds between distinct residues and record each partner in
    // both directions for later side-chain placement.
    DisulfidePartners partners;
    for (const auto* bond : mol.bonds()) {
        const auto* first = bond->getBeginAtom();
        const auto* second = bond->getEndAtom();
        if (first->getAtomicNum() != 16 || second->getAtomicNum() != 16) {
            continue;
        }
        const auto first_residue = atom_residue.find(first->getIdx());
        const auto second_residue = atom_residue.find(second->getIdx());
        if (first_residue == atom_residue.end() ||
            second_residue == atom_residue.end() ||
            first_residue->second == second_residue->second) {
            continue;
        }
        partners.emplace(first_residue->second, second_residue->second);
        partners.emplace(second_residue->second, first_residue->second);
    }
    return partners;
}

struct PeptideTopology {
    PeptideLinks links;
    ResidueKeys incoming;
};

PeptideTopology find_peptide_links(const RDKit::ROMol& mol,
                                   const RDKit::ROMol& monomer_mol,
                                   const RDKit::Conformer& monomer_conf,
                                   const ResidueMap& residues,
                                   const AtomResidueMap& atom_residue)
{
    // Find C-N peptide bonds between mapped residues in the same peptide
    // chain. Record each link's direction in the monomer layout and which
    // residues have an incoming peptide bond.
    PeptideLinks peptide_links;
    ResidueKeys incoming_peptide_links;
    for (const auto* bond : mol.bonds()) {
        const auto* begin = bond->getBeginAtom();
        const auto* end = bond->getEndAtom();
        const RDKit::Atom* carbonyl = nullptr;
        const RDKit::Atom* nitrogen = nullptr;
        if (has_atom_name(begin, "C") && has_atom_name(end, "N")) {
            carbonyl = begin;
            nitrogen = end;
        } else if (has_atom_name(end, "C") && has_atom_name(begin, "N")) {
            carbonyl = end;
            nitrogen = begin;
        } else {
            continue;
        }
        const auto carbonyl_key = atom_residue.find(carbonyl->getIdx());
        const auto nitrogen_key = atom_residue.find(nitrogen->getIdx());
        if (carbonyl_key == atom_residue.end() ||
            nitrogen_key == atom_residue.end() ||
            carbonyl_key->second == nitrogen_key->second) {
            continue;
        }
        const auto& from = residues.at(carbonyl_key->second);
        const auto& to = residues.at(nitrogen_key->second);
        if (!from.monomer_idx || !to.monomer_idx) {
            continue;
        }
        const auto* from_monomer =
            monomer_mol.getAtomWithIdx(*from.monomer_idx);
        const auto* to_monomer = monomer_mol.getAtomWithIdx(*to.monomer_idx);
        if (getChainType(*from_monomer) != ChainType::PEPTIDE ||
            getChainType(*to_monomer) != ChainType::PEPTIDE ||
            get_polymer_id(from_monomer) != get_polymer_id(to_monomer)) {
            continue;
        }
        const auto delta = monomer_conf.getAtomPos(*to.monomer_idx) -
                           monomer_conf.getAtomPos(*from.monomer_idx);
        if (delta.x == 0 && delta.y == 0) {
            continue;
        }
        const auto distance = std::hypot(delta.x, delta.y);
        const RDGeom::Point3D axis(delta.x / distance, delta.y / distance, 0);
        if (peptide_links
                .emplace(carbonyl_key->second,
                         PeptideLink{nitrogen_key->second, axis})
                .second) {
            incoming_peptide_links.insert(nitrogen_key->second);
        }
    }

    return {std::move(peptide_links), std::move(incoming_peptide_links)};
}

struct PeptideTurn {
    size_t first_residue;
    size_t residue_count;
    double turn_side;
    RDGeom::Point3D incoming_axis;
    RDGeom::Point3D incoming_frame;
    RDGeom::Point3D outgoing_frame;
};

void open_turn_angles(const RDKit::ROMol& mol, const ResidueMap& residues,
                      const std::vector<ResidueKey>& chain,
                      const PeptideTurn& turn, double target_turn_angle,
                      size_t spread_residue_count, RDKit::Conformer& output)
{
    // Adjust a turn to target_turn_angle (radians) across the requested
    // residues. Rotate only downstream atoms so an S-S bridge cannot pull
    // the other arm into the turn.
    if (turn.first_residue >= chain.size() || spread_residue_count == 0 ||
        spread_residue_count > chain.size() - turn.first_residue) {
        return;
    }
    const auto atom = [&](size_t residue, std::string_view name) {
        return *named_atom(mol, residues.at(chain[residue]), name);
    };
    // A neighboring proline or cysteine may have LINEAR_CIS_N_CA geometry,
    // so its N-C direction is not a reliable measure of the straight arm.
    const auto& first = turn.incoming_frame;
    const auto& last = turn.outgoing_frame;
    const auto actual = std::atan2(first.x * last.y - first.y * last.x,
                                   first.x * last.x + first.y * last.y);
    const auto correction =
        std::remainder(target_turn_angle - actual, 2 * std::numbers::pi);
    // Applying small rotations keeps each backbone corner close to its
    // original direction.
    constexpr size_t bonds_per_residue = 3; // N-CA, CA-C, C-N(next)
    const auto bonds_in_turn =
        (spread_residue_count - 1) * bonds_per_residue + 1;
    const auto angle = correction / bonds_in_turn;
    for (size_t step = 0; step < bonds_in_turn; ++step) {
        const auto residue_index =
            turn.first_residue + step / bonds_per_residue;
        const auto bond = step % bonds_per_residue;
        const auto pivot = atom(residue_index, bond == 0   ? "N"
                                               : bond == 1 ? "CA"
                                                           : "C");
        const auto downstream = bond == 0   ? atom(residue_index, "CA")
                                : bond == 1 ? atom(residue_index, "C")
                                            : atom(residue_index + 1, "N");
        // Move later peptide residues and the branch beyond this bond within
        // the current residue. Traversing the whole graph would follow an
        // S-S bridge back to the other arm and fold the chain onto itself.
        std::set<unsigned int> moving;
        for (size_t later = residue_index + 1; later < chain.size(); ++later) {
            moving.insert(residues.at(chain[later]).atom_ids.begin(),
                          residues.at(chain[later]).atom_ids.end());
        }
        if (bond != 2) {
            const auto& local_atoms =
                residues.at(chain[residue_index]).atom_ids;
            std::vector<unsigned int> pending{downstream};
            while (!pending.empty()) {
                const auto id = pending.back();
                pending.pop_back();
                if (id == pivot ||
                    !std::ranges::binary_search(local_atoms, id) ||
                    !moving.insert(id).second) {
                    continue;
                }
                for (const auto* neighbor :
                     mol.atomNeighbors(mol.getAtomWithIdx(id))) {
                    pending.push_back(neighbor->getIdx());
                }
            }
        }
        const auto center = output.getAtomPos(pivot);
        for (const auto id : moving) {
            output.setAtomPos(
                id, rotate_around(output.getAtomPos(id), center, angle));
        }
    }
}

std::vector<PeptideTurn>
find_peptide_turns(const std::vector<ResidueKey>& chain,
                   const PeptideLinks& links)
{
    // Find three- or four-residue turns connecting nearly opposite straight
    // arms. Require each monomer link within the turn to bend smoothly in
    // the same direction.
    std::vector<PeptideTurn> turns;
    if (chain.size() < 6) {
        return turns;
    }
    const auto axis = [&](size_t index) -> const RDGeom::Point3D& {
        return links.at(chain[index]).axis;
    };
    const auto dot = [](const auto& a, const auto& b) {
        return a.x * b.x + a.y * b.y;
    };
    const auto cross = [](const auto& a, const auto& b) {
        return a.x * b.y - a.y * b.x;
    };
    // Include the last residue of the incoming arm and the first of the
    // returning arm; one or two residues lie between them.
    for (size_t first = 1; first + 3 < chain.size(); ++first) {
        const auto& incoming = axis(first - 1);
        // A short upper arm has only one link before the bend. When a
        // second link exists, require both to point in the same direction.
        if (first > 1 &&
            dot(axis(first - 2), incoming) < std::cos(std::numbers::pi / 8)) {
            continue;
        }
        const auto turn_side = cross(incoming, axis(first)) < 0 ? -1.0 : 1.0;
        for (size_t count : {size_t{3}, size_t{4}}) {
            if (first + count >= chain.size()) {
                continue;
            }
            const auto& outgoing = axis(first + count - 1);
            if (dot(incoming, outgoing) > -std::cos(std::numbers::pi / 6)) {
                continue;
            }
            bool smooth_bend = true;
            for (size_t link = first; link < first + count; ++link) {
                const auto turn = std::atan2(cross(axis(link - 1), axis(link)),
                                             dot(axis(link - 1), axis(link)));
                if (turn_side * turn < std::numbers::pi / 8 ||
                    turn_side * turn > 2 * std::numbers::pi / 3) {
                    smooth_bend = false;
                    break;
                }
            }
            if (!smooth_bend) {
                continue;
            }
            turns.push_back({first, count, turn_side, incoming, {}, {}});
            first += count - 1;
            break;
        }
    }
    return turns;
}

void mirror_frame_across_bond(const RDGeom::Point3D& bond,
                              RDGeom::Point3D& along, RDGeom::Point3D& normal)
{
    // Update only the directions of bonds still to be drawn. No placed atom
    // needs to be reflected when a bend chooses CIS_N_CA or CIS_CA_CO.
    const RDGeom::Point3D origin(0, 0, 0);
    along = reflect_across_bond(along, origin, bond);
    normal = reflect_across_bond(normal, origin, bond);
}

ResidueShape straight_shape_facing_disulfide(
    const ResidueKey& key, ResidueShape initial_shape,
    const DisulfidePartners& partners, const ResidueMap& residues,
    const RDKit::Conformer& monomer_conf, const RDGeom::Point3D& frame_normal,
    double previous_side)
{
    // On a straight arm, use LINEAR_CIS_N_CA when an ALL_TRANS cysteine's
    // side chain would point away from its disulfide partner. Account for
    // earlier shape choices that may have reversed the side-chain pattern.
    if (initial_shape != ResidueShape::ALL_TRANS) {
        return initial_shape;
    }
    const auto partner = partners.find(key);
    if (partner == partners.end()) {
        return initial_shape;
    }
    const auto& residue = residues.at(key);
    const auto& other = residues.at(partner->second);
    if (!residue.monomer_idx || !other.monomer_idx) {
        return initial_shape;
    }
    const auto toward_partner = monomer_conf.getAtomPos(*other.monomer_idx) -
                                monomer_conf.getAtomPos(*residue.monomer_idx);
    const auto trans_outward = frame_normal * -previous_side;
    const auto points_away = trans_outward.x * toward_partner.x +
                                 trans_outward.y * toward_partner.y <
                             0;
    return points_away ? ResidueShape::LINEAR_CIS_N_CA
                       : ResidueShape::ALL_TRANS;
}

bool bend_contains_proline(const RDKit::ROMol& mol, const ResidueMap& residues,
                           const std::vector<ResidueKey>& chain,
                           const PeptideTurn& turn)
{
    // Check whether this turn contains proline or D-proline, whose ring
    // requires a turn at N-CA rather than at CA-C.
    for (size_t i = turn.first_residue;
         i < turn.first_residue + turn.residue_count; ++i) {
        if (initial_residue_shape(mol, residues.at(chain[i])) ==
            ResidueShape::LINEAR_CIS_N_CA) {
            return true;
        }
    }
    return false;
}

struct PeptideChain {
    std::vector<ResidueKey> residues;
    std::optional<ResidueKey> n_terminal_cap;
};

PeptideChain find_peptide_chain(const RDKit::ROMol& mol,
                                const ResidueMap& residues,
                                const PeptideLinks& links,
                                const ResidueKey& start)
{
    // Follow peptide links from the chain start. If the first residue has
    // C and O but no backbone N, save it as an N-terminal carbonyl cap and
    // remove it from the amino acids drawn by the backbone loop.
    PeptideChain chain;
    ResidueKeys seen;
    auto current = start;
    while (seen.insert(current).second) {
        chain.residues.push_back(current);
        const auto link = links.find(current);
        if (link == links.end()) {
            break;
        }
        current = link->second.next;
    }
    if (chain.residues.size() < 2) {
        return chain;
    }
    const auto& first = residues.at(chain.residues.front());
    if (!named_atom(mol, first, "N") && named_atom(mol, first, "C") &&
        named_atom(mol, first, "O")) {
        chain.n_terminal_cap = chain.residues.front();
        chain.residues.erase(chain.residues.begin());
    }
    return chain;
}

struct PeptideShapeState {
    size_t next_turn = 0;
    size_t bend_end = 0;
    std::optional<ResidueShape> bend_shape;
};

ResidueShape choose_peptide_shape(
    const RDKit::ROMol& mol, const RDKit::Conformer& monomer_conf,
    const ResidueMap& residues, const DisulfidePartners& disulfide_partners,
    const std::vector<ResidueKey>& chain, std::vector<PeptideTurn>& turns,
    size_t index, double previous_side, const RDGeom::Point3D& previous_outward,
    const RDGeom::Point3D& frame_along, const RDGeom::Point3D& frame_normal,
    PeptideShapeState& state)
{
    // Use the same CIS_N_CA or CIS_CA_CO shape for every residue in a bend.
    // Choose it from the bend direction relative to the preceding residue's
    // side chain. On straight sections, proline and disulfides may change
    // the side-chain phase.
    if (index == state.bend_end) {
        state.bend_shape.reset();
    }
    if (state.next_turn < turns.size() &&
        index == turns[state.next_turn].first_residue) {
        auto& turn = turns[state.next_turn++];
        turn.incoming_frame = frame_along;
        const auto cross = turn.incoming_axis.x * previous_outward.y -
                           turn.incoming_axis.y * previous_outward.x;
        // Opposite sides use CIS_N_CA; the same side uses CIS_CA_CO.
        state.bend_shape = cross * turn.turn_side < 0 ? ResidueShape::CIS_N_CA
                                                      : ResidueShape::CIS_CA_CO;
        state.bend_end = index + turn.residue_count;
        // Proline can use only CIS_N_CA in a bend because its ring closes
        // at N. Keep that shape for the entire bend.
        if (state.bend_shape == ResidueShape::CIS_CA_CO &&
            bend_contains_proline(mol, residues, chain, turn)) {
            state.bend_shape = ResidueShape::CIS_N_CA;
        }
    }
    const auto initial_shape =
        initial_residue_shape(mol, residues.at(chain[index]));
    auto straight_shape = straight_shape_facing_disulfide(
        chain[index], initial_shape, disulfide_partners, residues, monomer_conf,
        frame_normal, previous_side);
    if (!state.bend_shape && state.next_turn < turns.size() &&
        index + 1 == turns[state.next_turn].first_residue &&
        straight_shape == ResidueShape::ALL_TRANS) {
        const auto& turn = turns[state.next_turn];
        if (bend_contains_proline(mol, residues, chain, turn)) {
            const auto trans_outward = frame_normal * -previous_side;
            const auto cross = turn.incoming_axis.x * trans_outward.y -
                               turn.incoming_axis.y * trans_outward.x;
            // If the proline bend would require CIS_CA_CO, use
            // LINEAR_CIS_N_CA here to reverse the side-chain phase before
            // the bend. The bend can then use CIS_N_CA throughout.
            if (cross * turn.turn_side >= 0) {
                straight_shape = ResidueShape::LINEAR_CIS_N_CA;
            }
        }
    }
    auto shape = state.bend_shape.value_or(straight_shape);
    // A terminal cysteine needs no outgoing peptide bond to keep straight.
    if (index + 1 == chain.size() && !state.bend_shape &&
        shape == ResidueShape::LINEAR_CIS_N_CA &&
        initial_shape == ResidueShape::ALL_TRANS) {
        shape = ResidueShape::CIS_N_CA;
    }
    return shape;
}

bool place_peptide_n_terminal_cap(
    const RDKit::ROMol& mol, const ResidueMap& residues, const ResidueKey& cap,
    const RDGeom::Point3D& first_nitrogen, const RDGeom::Point3D& along,
    double first_side, bool first_linear_cis_n_ca, RDKit::Conformer& output)
{
    // Match the cap's C-N angle to the first residue's side-chain phase.
    const auto bond_length = RDDepict::BOND_LEN * BACKBONE_BOND_SCALE;
    const auto run = bond_length * std::sqrt(3.0) / 2;
    const auto rise = bond_length / 2;
    const auto small_bond_length = RDDepict::BOND_LEN * SIDECHAIN_BOND_SCALE;
    const RDGeom::Point3D normal(-along.y, along.x, 0);
    const auto cap_run = first_linear_cis_n_ca ? rise : run;
    const auto cap_rise = first_linear_cis_n_ca ? run : rise;
    const auto carbonyl_position =
        first_nitrogen - along * cap_run + normal * (first_side * cap_rise);
    const auto oxygen_position =
        carbonyl_position + normal * (first_side * small_bond_length);
    return place_n_terminal_carbonyl_cap(
        mol, residues.at(cap), carbonyl_position, oxygen_position, output);
}

struct PeptideBackbonePositions {
    RDGeom::Point3D alpha;
    RDGeom::Point3D carbonyl;
    RDGeom::Point3D oxygen;
    RDGeom::Point3D next_nitrogen;
    RDGeom::Point3D outward;
};

PeptideBackbonePositions next_peptide_backbone_positions(
    ResidueShape shape, double side, const RDGeom::Point3D& nitrogen,
    RDGeom::Point3D& frame_along, RDGeom::Point3D& frame_normal)
{
    // Calculate CA, carbonyl C, O, and the next residue's N from this N.
    // Update the frame after N-CA or CA-C when the shape calls for a turn.
    const auto bond_length = RDDepict::BOND_LEN * BACKBONE_BOND_SCALE;
    const auto small_bond_length = RDDepict::BOND_LEN * SIDECHAIN_BOND_SCALE;
    const auto run = bond_length * std::sqrt(3.0) / 2;
    const auto rise = bond_length / 2;
    const auto linear = shape == ResidueShape::LINEAR_CIS_N_CA;
    const auto alpha =
        nitrogen + (linear ? frame_along * bond_length
                           : frame_along * run + frame_normal * (side * rise));
    if (shape == ResidueShape::CIS_N_CA) {
        mirror_frame_across_bond(alpha - nitrogen, frame_along, frame_normal);
    }
    const auto outward = frame_normal * side;
    const auto carbonyl =
        alpha + frame_along * run - frame_normal * (side * rise);
    if (shape == ResidueShape::CIS_CA_CO) {
        mirror_frame_across_bond(carbonyl - alpha, frame_along, frame_normal);
    }
    const auto oxygen = carbonyl - frame_normal * (side * small_bond_length);
    const auto next_nitrogen =
        carbonyl + frame_along * run + frame_normal * (side * rise);
    return {alpha, carbonyl, oxygen, next_nitrogen, outward};
}

ResidueKeys place_peptide_chains(
    const RDKit::ROMol& mol, const RDKit::Conformer& monomer_conf,
    const ResidueMap& residues, const PeptideTopology& topology,
    const DisulfidePartners& disulfide_partners, double layout_scale,
    const ResidueKeys& flipped_chain_starts,
    LinearCounts& nonproline_linear_counts, RDKit::Conformer& output)
{
    // Draw each peptide chain along the monomer path. Choose the residue
    // shapes and side-chain phase, place backbone atoms and endcaps, then
    // distribute any remaining correction across four-residue turns.
    ResidueKeys placed_peptides;
    const auto& peptide_links = topology.links;
    const auto& incoming_peptide_links = topology.incoming;
    const auto peptide_bond_length = RDDepict::BOND_LEN * BACKBONE_BOND_SCALE;
    const auto run = peptide_bond_length * std::sqrt(3.0) / 2;
    // A residue with no incoming peptide link starts a chain.
    for (const auto& [start, first_link] : peptide_links) {
        if (incoming_peptide_links.contains(start)) {
            continue;
        }
        auto peptide_chain =
            find_peptide_chain(mol, residues, peptide_links, start);
        const auto& chain = peptide_chain.residues;
        if (chain.size() < 2 && !peptide_chain.n_terminal_cap) {
            continue;
        }
        const auto& first_residue = residues.at(chain.front());

        // Anchor the first N to the monomer layout and use the first link
        // to set the initial backbone direction.
        const auto along = first_link.axis;
        auto frame_along = along;
        RDGeom::Point3D frame_normal(-along.y, along.x, 0);
        size_t placed_count = 0;
        auto turns = find_peptide_turns(chain, peptide_links);
        const auto first_nitrogen_position =
            monomer_conf.getAtomPos(*first_residue.monomer_idx) * layout_scale -
            along * run;
        auto nitrogen_position = first_nitrogen_position;
        double previous_side =
            flipped_chain_starts.contains(start) ? 1.0 : -1.0;
        double first_side = 0;
        bool first_linear_cis_n_ca = false;
        size_t chain_linear_count = 0;
        RDGeom::Point3D previous_outward(0, 0, 0);
        PeptideShapeState shape_state;
        // Walk the chain in peptide order. Each shape determines both the
        // current residue's side-chain side and the next backbone frame.
        for (size_t i = 0; i < chain.size(); ++i) {
            const auto& residue = residues.at(chain[i]);
            if (!named_atom(mol, residue, "N")) {
                break;
            }
            const auto initial_shape = initial_residue_shape(mol, residue);
            const auto shape = choose_peptide_shape(
                mol, monomer_conf, residues, disulfide_partners, chain, turns,
                i, previous_side, previous_outward, frame_along, frame_normal,
                shape_state);
            // Count non-proline phase changes when comparing starting phases.
            if (initial_shape == ResidueShape::ALL_TRANS &&
                shape == ResidueShape::LINEAR_CIS_N_CA) {
                ++chain_linear_count;
            }
            const auto linear_cis_n_ca = shape == ResidueShape::LINEAR_CIS_N_CA;
            const auto side = linear_cis_n_ca ? previous_side : -previous_side;
            if (i == 0) {
                first_side = side;
                first_linear_cis_n_ca = linear_cis_n_ca;
            }
            const auto positions = next_peptide_backbone_positions(
                shape, side, nitrogen_position, frame_along, frame_normal);
            const auto placed = place_peptide_backbone_residue(
                mol, residue, nitrogen_position, positions.alpha,
                positions.carbonyl, positions.oxygen, positions.outward,
                output);
            // A final amine has N but no CA or carbonyl; place it separately.
            if (!placed && !(i + 1 == chain.size() &&
                             place_terminal_amine(mol, residue,
                                                  nitrogen_position, output))) {
                break;
            }
            previous_side = side;
            previous_outward = positions.outward;
            // Save the outgoing direction for a possible angle correction.
            if (shape_state.bend_shape && i + 1 == shape_state.bend_end) {
                turns[shape_state.next_turn - 1].outgoing_frame = frame_along;
            }
            placed_peptides.insert(chain[i]);
            ++placed_count;
            if (!placed) {
                break;
            }
            nitrogen_position = positions.next_nitrogen;
        }

        // Position an N-terminal carbonyl cap after its first peptide N.
        if (peptide_chain.n_terminal_cap && placed_count != 0) {
            if (place_peptide_n_terminal_cap(
                    mol, residues, *peptide_chain.n_terminal_cap,
                    first_nitrogen_position, along, first_side,
                    first_linear_cis_n_ca, output)) {
                placed_peptides.insert(*peptide_chain.n_terminal_cap);
            }
        }

        // Four-residue bends need their angles adjusted to complete
        // a 180-degree turn between the two straight arms.
        if (placed_count == chain.size()) {
            for (const auto& turn : turns) {
                if (turn.residue_count == 4) {
                    open_turn_angles(mol, residues, chain, turn,
                                     turn.turn_side * std::numbers::pi,
                                     turn.residue_count, output);
                }
            }
        }
        nonproline_linear_counts[start] = chain_linear_count;
    }

    return placed_peptides;
}

std::optional<RDGeom::Point3D>
free_trigonal_direction(const RDGeom::Point3D& center,
                        const RDGeom::Point3D& first_neighbor,
                        std::optional<RDGeom::Point3D> second_neighbor,
                        std::optional<RDGeom::Point3D> avoid)
{
    // With two bonds, use the third trigonal direction. With one bond,
    // choose the 120-degree direction farther from the atom to avoid.
    const auto first = first_neighbor - center;
    const auto first_length = std::hypot(first.x, first.y);
    if (first_length == 0) {
        return std::nullopt;
    }
    const auto direction = first / first_length;
    if (second_neighbor) {
        const auto second = *second_neighbor - center;
        const auto second_length = std::hypot(second.x, second.y);
        if (second_length == 0) {
            return std::nullopt;
        }
        const auto opposite = direction + second / second_length;
        const auto length = std::hypot(opposite.x, opposite.y);
        if (length == 0) {
            return std::nullopt;
        }
        return -opposite / length;
    }
    const RDGeom::Point3D origin(0, 0, 0);
    const auto left =
        rotate_around(direction, origin, 2 * std::numbers::pi / 3);
    if (!avoid) {
        return left;
    }
    const auto right =
        rotate_around(direction, origin, -2 * std::numbers::pi / 3);
    const auto away = *avoid - center;
    return right.x * away.x + right.y * away.y <
                   left.x * away.x + left.y * away.y
               ? right
               : left;
}

bool place_unnamed_terminal_residue(const RDKit::ROMol& mol,
                                    const ResidueGeometry& residue,
                                    const ResidueMap& residues,
                                    const AtomResidueMap& atom_residue,
                                    const ResidueKeys& placed_peptides,
                                    RDKit::Conformer& output)
{
    // Find an unnamed terminal group's bond to an already placed peptide.
    // Place a single H or N from the free direction at that bond; for a
    // carbonyl cap, align its C=O with the peptide backbone.
    for (auto atom_id : residue.atom_ids) {
        const auto* terminal_atom = mol.getAtomWithIdx(atom_id);
        for (const auto* peptide_atom : mol.atomNeighbors(terminal_atom)) {
            const auto found = atom_residue.find(peptide_atom->getIdx());
            if (found == atom_residue.end() ||
                !placed_peptides.contains(found->second)) {
                continue;
            }
            const auto& peptide = residues.at(found->second);

            // A one-atom H cap attaches to the first peptide N.
            if (residue.atom_ids.size() == 1 &&
                terminal_atom->getAtomicNum() == 1 &&
                has_atom_name(peptide_atom, "N")) {
                const auto alpha = named_atom(mol, peptide, "CA");
                if (!alpha) {
                    continue;
                }
                const auto nitrogen_position =
                    output.getAtomPos(peptide_atom->getIdx());
                const auto carbonyl = named_atom(mol, peptide, "C");
                std::optional<RDGeom::Point3D> avoid;
                if (carbonyl) {
                    avoid = output.getAtomPos(*carbonyl);
                }
                const auto outward = free_trigonal_direction(
                    nitrogen_position, output.getAtomPos(*alpha), std::nullopt,
                    avoid);
                if (!outward) {
                    continue;
                }
                output.setAtomPos(
                    atom_id,
                    nitrogen_position +
                        *outward * (RDDepict::BOND_LEN * SIDECHAIN_BOND_SCALE));
                return true;
            }

            // A one-atom N cap fills the third direction at peptide C=O.
            if (residue.atom_ids.size() == 1 &&
                terminal_atom->getAtomicNum() == 7 &&
                has_atom_name(peptide_atom, "C")) {
                const auto alpha = named_atom(mol, peptide, "CA");
                const auto oxygen = named_atom(mol, peptide, "O");
                if (!alpha || !oxygen) {
                    continue;
                }
                const auto carbonyl_position =
                    output.getAtomPos(peptide_atom->getIdx());
                const auto outward = free_trigonal_direction(
                    carbonyl_position, output.getAtomPos(*alpha),
                    output.getAtomPos(*oxygen), std::nullopt);
                if (!outward) {
                    continue;
                }
                output.setAtomPos(
                    atom_id,
                    carbonyl_position +
                        *outward * (RDDepict::BOND_LEN * BACKBONE_BOND_SCALE));
                return true;
            }

            // An inline SMILES carbonyl cap has no PDB C/O names. Identify
            // its C by the bond to peptide N and its O by the double bond.
            if (terminal_atom->getAtomicNum() != 6 ||
                has_atom_name(terminal_atom, "C") ||
                !has_atom_name(peptide_atom, "N")) {
                continue;
            }
            const auto next_alpha = named_atom(mol, peptide, "CA");
            const auto next_carbonyl = named_atom(mol, peptide, "C");
            if (!next_alpha || !next_carbonyl) {
                continue;
            }
            for (const auto* oxygen : mol.atomNeighbors(terminal_atom)) {
                if (oxygen->getAtomicNum() != 8 ||
                    !std::ranges::binary_search(residue.atom_ids,
                                                oxygen->getIdx()) ||
                    mol.getBondBetweenAtoms(atom_id, oxygen->getIdx())
                            ->getBondType() != RDKit::Bond::DOUBLE) {
                    continue;
                }
                const auto nitrogen_position =
                    output.getAtomPos(peptide_atom->getIdx());
                auto along =
                    output.getAtomPos(*next_carbonyl) - nitrogen_position;
                const auto length = std::hypot(along.x, along.y);
                if (length == 0) {
                    continue;
                }
                along /= length;
                const RDGeom::Point3D normal(-along.y, along.x, 0);
                const auto n_ca =
                    output.getAtomPos(*next_alpha) - nitrogen_position;
                const auto side =
                    n_ca.x * normal.x + n_ca.y * normal.y >= 0 ? 1.0 : -1.0;
                const auto bond_length =
                    RDDepict::BOND_LEN * BACKBONE_BOND_SCALE;
                // The cap's C occupies the same free direction as a
                // one-atom H cap, with a full backbone bond length.
                const auto cap_direction = free_trigonal_direction(
                    nitrogen_position, output.getAtomPos(*next_alpha),
                    std::nullopt, output.getAtomPos(*next_carbonyl));
                if (!cap_direction) {
                    continue;
                }
                const auto target_carbonyl =
                    nitrogen_position + *cap_direction * bond_length;
                const auto target_oxygen =
                    target_carbonyl +
                    normal * (side * RDDepict::BOND_LEN * SIDECHAIN_BOND_SCALE);
                if (!place_residue_by_bond(residue, atom_id, oxygen->getIdx(),
                                           target_carbonyl, target_oxygen,
                                           output)) {
                    return false;
                }
                // Keep the fragment's carbon chain behind the C-N bond.
                for (const auto* neighbor : mol.atomNeighbors(terminal_atom)) {
                    if (neighbor->getAtomicNum() != 6 ||
                        !std::ranges::binary_search(residue.atom_ids,
                                                    neighbor->getIdx())) {
                        continue;
                    }
                    const auto branch =
                        output.getAtomPos(neighbor->getIdx()) - target_carbonyl;
                    if (branch.x * along.x + branch.y * along.y > 0) {
                        for (auto residue_atom : residue.atom_ids) {
                            output.setAtomPos(
                                residue_atom,
                                reflect_across_bond(
                                    output.getAtomPos(residue_atom),
                                    target_carbonyl, target_oxygen));
                        }
                    }
                    break;
                }
                return true;
            }
        }
    }
    return false;
}

void place_remaining_residues(const RDKit::ROMol& mol,
                              const RDKit::ROMol& monomer_mol,
                              const RDKit::Conformer& monomer_conf,
                              const ResidueMap& residues,
                              const AtomResidueMap& atom_residue,
                              const ResidueKeys& placed_peptides,
                              double layout_scale, RDKit::Conformer& output)
{
    // Place residues that the peptide backbone pass did not handle. Attach
    // unnamed terminal groups to placed peptides; otherwise rotate each
    // residue toward its neighbors in the monomer layout. Put groups with
    // no monomer mapping beside the layout.
    double fallback_x = 0;
    for (const auto* atom : monomer_mol.atoms()) {
        const auto p = monomer_conf.getAtomPos(atom->getIdx());
        fallback_x =
            std::max(fallback_x, p.x * layout_scale + MONOMER_BOND_LENGTH);
    }
    for (const auto& [key, residue] : residues) {
        if (placed_peptides.contains(key)) {
            continue;
        }
        if (place_unnamed_terminal_residue(mol, residue, residues, atom_residue,
                                           placed_peptides, output)) {
            continue;
        }
        RDGeom::Point3D target;
        if (residue.monomer_idx) {
            target =
                monomer_conf.getAtomPos(*residue.monomer_idx) * layout_scale;
        } else {
            // Conversion omits some PDB groups, notably waters. Keep their
            // atoms in the output, placing those groups beside the layout.
            target = {fallback_x + residue.size.x / 2, 0, 0};
            fallback_x += residue.size.x + MONOMER_BOND_LENGTH;
        }

        // Match the direction of each inter-residue attachment to its
        // neighbor's direction in the monomer layout, as in the atomistic
        // monomer preview's circular mean of rotation angles.
        double cos_sum = 0, sin_sum = 0;
        if (residue.monomer_idx) {
            for (size_t i = 0; i < residue.atom_ids.size(); ++i) {
                const auto* atom = mol.getAtomWithIdx(residue.atom_ids[i]);
                const auto source = residue.positions[i] - residue.center;
                if (source.x == 0 && source.y == 0) {
                    continue;
                }
                for (const auto* neighbor : mol.atomNeighbors(atom)) {
                    const auto neighbor_key =
                        atom_residue.find(neighbor->getIdx());
                    if (neighbor_key == atom_residue.end() ||
                        neighbor_key->second == key) {
                        continue;
                    }
                    const auto& neighbor_residue =
                        residues.at(neighbor_key->second);
                    if (!neighbor_residue.monomer_idx) {
                        continue;
                    }
                    const auto direction =
                        monomer_conf.getAtomPos(*neighbor_residue.monomer_idx) -
                        target;
                    if (direction.x == 0 && direction.y == 0) {
                        continue;
                    }
                    const auto angle = std::atan2(direction.y, direction.x) -
                                       std::atan2(source.y, source.x);
                    cos_sum += std::cos(angle);
                    sin_sum += std::sin(angle);
                }
            }
        }
        const auto angle = std::atan2(sin_sum, cos_sum);
        const auto cosine = std::cos(angle), sine = std::sin(angle);
        for (size_t i = 0; i < residue.atom_ids.size(); ++i) {
            const auto p =
                (residue.positions[i] - residue.center) * BACKBONE_BOND_SCALE;
            output.setAtomPos(residue.atom_ids[i],
                              {target.x + p.x * cosine - p.y * sine,
                               target.y + p.x * sine + p.y * cosine, 0});
        }
    }
}

void place_backbone_hydrogens(const RDKit::ROMol& mol,
                              const ResidueMap& residues,
                              const ResidueKeys& placed_peptides,
                              RDKit::Conformer& output)
{
    // Correct explicit H atoms on CA and amide N after the complete backbone
    // has been placed, using the final positions of their heavy neighbors.
    const auto small_bond_length = RDDepict::BOND_LEN * SIDECHAIN_BOND_SCALE;
    for (const auto& key : placed_peptides) {
        const auto& residue = residues.at(key);
        const auto nitrogen = named_atom(mol, residue, "N");
        if (!nitrogen) {
            continue;
        }
        const auto alpha = named_atom(mol, residue, "CA");
        const auto carbonyl = named_atom(mol, residue, "C");
        if (alpha && carbonyl) {
            place_alpha_hydrogen(mol, residue, *alpha,
                                 output.getAtomPos(*nitrogen),
                                 output.getAtomPos(*alpha),
                                 output.getAtomPos(*carbonyl), output);
        }
        std::vector<unsigned int> heavy_neighbors;
        std::vector<unsigned int> hydrogens;
        for (const auto* neighbor :
             mol.atomNeighbors(mol.getAtomWithIdx(*nitrogen))) {
            if (neighbor->getAtomicNum() == 1) {
                hydrogens.push_back(neighbor->getIdx());
            } else {
                heavy_neighbors.push_back(neighbor->getIdx());
            }
        }
        // Leave terminal NH2/NH3 and substituted nitrogens to their fragment
        // depiction; only an amide N with two heavy bonds has one bisector.
        if (heavy_neighbors.size() != 2 || hydrogens.size() != 1) {
            continue;
        }
        const auto nitrogen_pos = output.getAtomPos(*nitrogen);
        RDGeom::Point3D bisector(0, 0, 0);
        for (auto atom_id : heavy_neighbors) {
            const auto bond = output.getAtomPos(atom_id) - nitrogen_pos;
            const auto length = std::hypot(bond.x, bond.y);
            if (length > 0) {
                bisector += bond / length;
            }
        }
        const auto length = std::hypot(bisector.x, bisector.y);
        if (length > 0) {
            output.setAtomPos(hydrogens.front(),
                              nitrogen_pos -
                                  bisector * (small_bond_length / length));
        }
    }
}

constexpr size_t SIDECHAIN_VARIANTS = 4;
constexpr double SIDECHAIN_TURN_ANGLE = std::numbers::pi / 12; // 15 degrees
constexpr double CLASH_DISTANCE_FRACTION = 0.8;
constexpr double ATOM_CLASH_WEIGHT = 20;
constexpr double BOND_CROSSING_WEIGHT = 100;
constexpr double DISULFIDE_LENGTH_WEIGHT = 2;

struct MovableSidechain {
    ResidueKey residue;
    unsigned int alpha;
    unsigned int beta;
    std::vector<unsigned int> atoms;
    std::vector<unsigned int> heavy_atoms;
    std::vector<const RDKit::Bond*> bonds;
    std::unordered_map<unsigned int, size_t> offsets;
    std::array<std::optional<std::vector<RDGeom::Point3D>>, SIDECHAIN_VARIANTS>
        coordinates;
    std::array<std::optional<double>, SIDECHAIN_VARIANTS> local_scores;
    bool allow_all_variants = false;
    size_t selected = 0;
};

std::optional<MovableSidechain>
make_movable_sidechain(const RDKit::ROMol& mol, const ResidueKey& key,
                       const ResidueGeometry& residue,
                       const RDKit::Conformer& output)
{
    // Collect the branch beyond CA for a residue whose side chain can move
    // without moving the backbone. Skip missing or ambiguous branches,
    // backbone-closing rings, and branches with only one heavy atom.
    const auto nitrogen = named_atom(mol, residue, "N");
    const auto alpha = named_atom(mol, residue, "CA");
    const auto carbonyl = named_atom(mol, residue, "C");
    if (!nitrogen || !alpha || !carbonyl) {
        return std::nullopt;
    }
    std::vector<unsigned int> roots;
    for (const auto* neighbor : mol.atomNeighbors(mol.getAtomWithIdx(*alpha))) {
        if (neighbor->getAtomicNum() != 1 && neighbor->getIdx() != *nitrogen &&
            neighbor->getIdx() != *carbonyl &&
            std::ranges::binary_search(residue.atom_ids, neighbor->getIdx())) {
            roots.push_back(neighbor->getIdx());
        }
    }
    if (roots.size() != 1) {
        return std::nullopt; // Gly, or an ambiguous alpha substituent.
    }

    MovableSidechain sidechain{key, *alpha, roots.front()};
    std::vector<unsigned int> pending{sidechain.beta};
    std::set<unsigned int> visited;
    while (!pending.empty()) {
        const auto id = pending.back();
        pending.pop_back();
        if (id == *nitrogen || id == *carbonyl) {
            return std::nullopt; // Proline or another backbone-closing ring.
        }
        if (id == *alpha || !std::ranges::binary_search(residue.atom_ids, id) ||
            !visited.insert(id).second) {
            continue;
        }
        for (const auto* neighbor : mol.atomNeighbors(mol.getAtomWithIdx(id))) {
            pending.push_back(neighbor->getIdx());
        }
    }
    sidechain.atoms.assign(visited.begin(), visited.end());
    for (size_t offset = 0; offset < sidechain.atoms.size(); ++offset) {
        const auto id = sidechain.atoms[offset];
        sidechain.offsets.emplace(id, offset);
        if (mol.getAtomWithIdx(id)->getAtomicNum() != 1) {
            sidechain.heavy_atoms.push_back(id);
        }
    }
    if (sidechain.heavy_atoms.size() <= 1) {
        return std::nullopt; // Ala's methyl cannot change its 2D shape.
    }
    for (const auto* bond : mol.bonds()) {
        const auto first = bond->getBeginAtomIdx();
        const auto second = bond->getEndAtomIdx();
        if ((first == *alpha || sidechain.offsets.contains(first)) &&
            (second == *alpha || sidechain.offsets.contains(second)) &&
            (sidechain.offsets.contains(first) ||
             sidechain.offsets.contains(second))) {
            sidechain.bonds.push_back(bond);
        }
    }

    // Keep the original branch. Other positions are built only if tried.
    auto& original = sidechain.coordinates[0].emplace();
    for (const auto id : sidechain.atoms) {
        original.push_back(output.getAtomPos(id));
    }
    return sidechain;
}

void apply_sidechain_variant(MovableSidechain& sidechain, size_t variant,
                             RDKit::Conformer& output)
{
    // Generate a side-chain variant only on first use: mirror the branch
    // across CA-CB for variant 1, or rotate it around CA for variants 2/3.
    // Then copy the cached coordinates into the output conformer.
    auto& positions = sidechain.coordinates[variant];
    if (!positions) {
        const auto& original = *sidechain.coordinates[0];
        const auto alpha_position = output.getAtomPos(sidechain.alpha);
        const auto beta_position =
            original[sidechain.offsets.at(sidechain.beta)];
        auto& generated = positions.emplace();
        generated.reserve(original.size());
        for (size_t i = 0; i < original.size(); ++i) {
            const auto& position = original[i];
            if (variant == 1) {
                generated.push_back(sidechain.atoms[i] == sidechain.beta
                                        ? position
                                        : reflect_across_bond(position,
                                                              alpha_position,
                                                              beta_position));
            } else {
                generated.push_back(rotate_around(position, alpha_position,
                                                  variant == 2
                                                      ? -SIDECHAIN_TURN_ANGLE
                                                      : SIDECHAIN_TURN_ANGLE));
            }
        }
    }
    for (size_t i = 0; i < sidechain.atoms.size(); ++i) {
        output.setAtomPos(sidechain.atoms[i], (*positions)[i]);
    }
}

double cross_2d(const RDGeom::Point3D& a, const RDGeom::Point3D& b)
{
    // Return the signed area of the parallelogram spanned by a and b.
    return a.x * b.y - a.y * b.x;
}

bool bonds_cross(const RDKit::Bond& first, const RDKit::Bond& second,
                 const RDKit::Conformer& output)
{
    // Count a crossing only when the bonds share no atom and each segment's
    // endpoints lie on opposite sides of the other segment. A bounding-box
    // check skips pairs that cannot meet; touching endpoints do not count.
    const auto a = first.getBeginAtomIdx();
    const auto b = first.getEndAtomIdx();
    const auto c = second.getBeginAtomIdx();
    const auto d = second.getEndAtomIdx();
    if (a == c || a == d || b == c || b == d) {
        return false;
    }
    const auto pa = output.getAtomPos(a);
    const auto pb = output.getAtomPos(b);
    const auto pc = output.getAtomPos(c);
    const auto pd = output.getAtomPos(d);
    if (std::max(pa.x, pb.x) <= std::min(pc.x, pd.x) ||
        std::max(pc.x, pd.x) <= std::min(pa.x, pb.x) ||
        std::max(pa.y, pb.y) <= std::min(pc.y, pd.y) ||
        std::max(pc.y, pd.y) <= std::min(pa.y, pb.y)) {
        return false;
    }
    const auto ab = pb - pa, cd = pd - pc;
    constexpr double epsilon = 1e-8;
    const auto first_side = cross_2d(ab, pc - pa);
    const auto second_side = cross_2d(ab, pd - pa);
    const auto third_side = cross_2d(cd, pa - pc);
    const auto fourth_side = cross_2d(cd, pb - pc);
    return ((first_side > epsilon && second_side < -epsilon) ||
            (first_side < -epsilon && second_side > epsilon)) &&
           ((third_side > epsilon && fourth_side < -epsilon) ||
            (third_side < -epsilon && fourth_side > epsilon));
}

bool is_backbone_atom(const RDKit::Atom* atom)
{
    // Recognize peptide backbone atoms by their PDB names, including OXT.
    return has_atom_name(atom, "N") || has_atom_name(atom, "CA") ||
           has_atom_name(atom, "C") || has_atom_name(atom, "O") ||
           has_atom_name(atom, "OXT");
}

double atom_clash_penalty(const RDGeom::Point3D& first,
                          const RDGeom::Point3D& second,
                          double minimum_distance)
{
    // Skip distant pairs with an axis check. Within minimum_distance,
    // penalize the square of their fractional overlap.
    const auto delta = first - second;
    if (std::abs(delta.x) >= minimum_distance ||
        std::abs(delta.y) >= minimum_distance) {
        return 0;
    }
    const auto distance = std::hypot(delta.x, delta.y);
    if (distance >= minimum_distance) {
        return 0;
    }
    const auto overlap = 1 - distance / minimum_distance;
    return ATOM_CLASH_WEIGHT * overlap * overlap;
}

struct Bounds2D {
    double min_x, max_x, min_y, max_y;
};

bool bounds_near(const Bounds2D& first, const Bounds2D& second,
                 double minimum_distance)
{
    // Two boxes can contain clashing atoms only if their gap is smaller
    // than the minimum allowed atom distance on both axes.
    return first.min_x <= second.max_x + minimum_distance &&
           second.min_x <= first.max_x + minimum_distance &&
           first.min_y <= second.max_y + minimum_distance &&
           second.min_y <= first.max_y + minimum_distance;
}

Bounds2D combined_bounds(const Bounds2D& first, const Bounds2D& second)
{
    // Return the smallest axis-aligned box containing both inputs.
    return {std::min(first.min_x, second.min_x),
            std::max(first.max_x, second.max_x),
            std::min(first.min_y, second.min_y),
            std::max(first.max_y, second.max_y)};
}

std::vector<Bounds2D>
possible_atom_bounds(const RDKit::ROMol& mol,
                     const std::vector<MovableSidechain>& sidechains,
                     const RDKit::Conformer& output)
{
    // Fixed atoms get point boxes. Each movable atom gets a square around
    // CA large enough to contain its rotations and reflection, without
    // generating any side-chain variants.
    std::vector<Bounds2D> bounds;
    bounds.reserve(mol.getNumAtoms());
    for (const auto* atom : mol.atoms()) {
        const auto point = output.getAtomPos(atom->getIdx());
        bounds.push_back({point.x, point.x, point.y, point.y});
    }
    for (const auto& sidechain : sidechains) {
        const auto center = output.getAtomPos(sidechain.alpha);
        for (const auto id : sidechain.atoms) {
            const auto point = output.getAtomPos(id);
            const auto radius =
                std::hypot(point.x - center.x, point.y - center.y) + 1e-8;
            bounds[id] = {center.x - radius, center.x + radius,
                          center.y - radius, center.y + radius};
        }
    }
    return bounds;
}

struct ClashPairs {
    std::vector<std::pair<unsigned int, unsigned int>> atoms;
    std::vector<std::pair<const RDKit::Bond*, const RDKit::Bond*>> bonds;
    std::vector<const RDKit::Bond*> disulfides;
    double minimum_distance;
    double sidechain_bond_length;
};

ClashPairs make_clash_pairs(const RDKit::ROMol& mol,
                            const std::vector<MovableSidechain>& sidechains,
                            const AtomResidueMap& atom_residue,
                            const RDKit::Conformer& output)
{
    // Prepare only atom and bond pairs that a side-chain variant could
    // change. Use bounds covering all variants to discard distant pairs;
    // keep disulfide bonds so their lengths can be scored separately.
    ClashPairs pairs;
    pairs.sidechain_bond_length = RDDepict::BOND_LEN * SIDECHAIN_BOND_SCALE;
    pairs.minimum_distance =
        pairs.sidechain_bond_length * CLASH_DISTANCE_FRACTION;
    const auto atom_bounds = possible_atom_bounds(mol, sidechains, output);
    std::vector<bool> movable(mol.getNumAtoms(), false);
    for (const auto& sidechain : sidechains) {
        for (const auto id : sidechain.atoms) {
            movable[id] = true;
        }
    }

    // Collect heavy atoms and bonds once, including bonds that may cross.
    std::vector<unsigned int> heavy_atoms;
    std::vector<const RDKit::Bond*> heavy_bonds;
    std::vector<Bounds2D> heavy_bond_bounds;
    for (const auto* atom : mol.atoms()) {
        if (atom->getAtomicNum() != 1) {
            heavy_atoms.push_back(atom->getIdx());
        }
    }
    for (const auto* bond : mol.bonds()) {
        if (bond->getBeginAtom()->getAtomicNum() == 1 ||
            bond->getEndAtom()->getAtomicNum() == 1) {
            continue;
        }
        heavy_bonds.push_back(bond);
        heavy_bond_bounds.push_back(
            combined_bounds(atom_bounds[bond->getBeginAtomIdx()],
                            atom_bounds[bond->getEndAtomIdx()]));
        if (bond->getBeginAtom()->getAtomicNum() == 16 &&
            bond->getEndAtom()->getAtomicNum() == 16 &&
            atom_residue.at(bond->getBeginAtomIdx()) !=
                atom_residue.at(bond->getEndAtomIdx())) {
            pairs.disulfides.push_back(bond);
        }
    }
    // Skip fixed pairs, bonded atoms, and atoms separated by one neighbor.
    for (size_t i = 0; i < heavy_atoms.size(); ++i) {
        for (size_t j = i + 1; j < heavy_atoms.size(); ++j) {
            const auto first = heavy_atoms[i];
            const auto second = heavy_atoms[j];
            if ((!movable[first] && !movable[second]) ||
                !bounds_near(atom_bounds[first], atom_bounds[second],
                             pairs.minimum_distance) ||
                mol.getBondBetweenAtoms(first, second)) {
                continue;
            }
            bool two_bonds_apart = false;
            for (const auto* neighbor :
                 mol.atomNeighbors(mol.getAtomWithIdx(first))) {
                if (mol.getBondBetweenAtoms(neighbor->getIdx(), second)) {
                    two_bonds_apart = true;
                    break;
                }
            }
            if (!two_bonds_apart) {
                pairs.atoms.emplace_back(first, second);
            }
        }
    }
    // Only bonds with a movable endpoint can gain or lose a crossing.
    for (size_t i = 0; i < heavy_bonds.size(); ++i) {
        for (size_t j = i + 1; j < heavy_bonds.size(); ++j) {
            const auto* first = heavy_bonds[i];
            const auto* second = heavy_bonds[j];
            if (!bounds_near(heavy_bond_bounds[i], heavy_bond_bounds[j], 0)) {
                continue;
            }
            if (movable[first->getBeginAtomIdx()] ||
                movable[first->getEndAtomIdx()] ||
                movable[second->getBeginAtomIdx()] ||
                movable[second->getEndAtomIdx()]) {
                pairs.bonds.emplace_back(first, second);
            }
        }
    }
    return pairs;
}

struct AffectedPairs {
    std::vector<size_t> atoms;
    std::vector<size_t> bonds;
    std::vector<size_t> disulfides;
};

void record_affected_pair(std::vector<AffectedPairs>& affected,
                          std::vector<size_t> AffectedPairs::*member,
                          size_t pair_index,
                          std::initializer_list<unsigned int> atom_ids,
                          const std::vector<size_t>& owner)
{
    // Add this pair to each distinct side chain that owns one of its atoms.
    // Ignore fixed atoms; member selects the atom, bond, or S-S pair list.
    std::array<size_t, 4> seen{};
    size_t count = 0;
    for (const auto id : atom_ids) {
        const auto sidechain = owner[id];
        if (sidechain == affected.size() ||
            std::find(seen.begin(), seen.begin() + count, sidechain) !=
                seen.begin() + count) {
            continue;
        }
        seen[count++] = sidechain;
        (affected[sidechain].*member).push_back(pair_index);
    }
}

std::vector<AffectedPairs>
index_affected_pairs(const RDKit::ROMol& mol,
                     const std::vector<MovableSidechain>& sidechains,
                     const ClashPairs& pairs)
{
    // Map each movable atom to its side chain, then index the atom clashes,
    // bond crossings, and disulfides affected by each side-chain move.
    // The sentinel owner value marks atoms that never move.
    std::vector<AffectedPairs> affected(sidechains.size());
    std::vector<size_t> owner(mol.getNumAtoms(), sidechains.size());
    for (size_t i = 0; i < sidechains.size(); ++i) {
        for (const auto atom : sidechains[i].atoms) {
            owner[atom] = i;
        }
    }
    for (size_t i = 0; i < pairs.atoms.size(); ++i) {
        const auto [first, second] = pairs.atoms[i];
        record_affected_pair(affected, &AffectedPairs::atoms, i,
                             {first, second}, owner);
    }
    for (size_t i = 0; i < pairs.bonds.size(); ++i) {
        const auto [first, second] = pairs.bonds[i];
        record_affected_pair(affected, &AffectedPairs::bonds, i,
                             {first->getBeginAtomIdx(), first->getEndAtomIdx(),
                              second->getBeginAtomIdx(),
                              second->getEndAtomIdx()},
                             owner);
    }
    for (size_t i = 0; i < pairs.disulfides.size(); ++i) {
        const auto* bond = pairs.disulfides[i];
        record_affected_pair(affected, &AffectedPairs::disulfides, i,
                             {bond->getBeginAtomIdx(), bond->getEndAtomIdx()},
                             owner);
    }
    return affected;
}

AffectedPairs combine_affected_pairs(const AffectedPairs& first,
                                     const AffectedPairs& second)
{
    // Merge the sorted pair indexes for two side chains, keeping each
    // shared atom, bond, or disulfide pair only once.
    AffectedPairs combined;
    const auto merge = [](const auto& a, const auto& b, auto& result) {
        std::set_union(a.begin(), a.end(), b.begin(), b.end(),
                       std::back_inserter(result));
    };
    merge(first.atoms, second.atoms, combined.atoms);
    merge(first.bonds, second.bonds, combined.bonds);
    merge(first.disulfides, second.disulfides, combined.disulfides);
    return combined;
}

double disulfide_length_penalty(const RDKit::Bond& bond,
                                const ClashPairs& pairs,
                                const RDKit::Conformer& output)
{
    // Penalize squared relative deviation from a side-chain bond length.
    // The small weight allows an S-S bridge to stretch when needed.
    const auto delta = output.getAtomPos(bond.getBeginAtomIdx()) -
                       output.getAtomPos(bond.getEndAtomIdx());
    const auto relative_error =
        std::hypot(delta.x, delta.y) / pairs.sidechain_bond_length - 1;
    return DISULFIDE_LENGTH_WEIGHT * relative_error * relative_error;
}

double clash_score(const AffectedPairs& affected, const ClashPairs& pairs,
                   const RDKit::Conformer& output)
{
    // Score only pairs touched by the moved side chains: atom overlap,
    // bond crossings, and disulfide bond-length error.
    double score = 0;
    for (const auto index : affected.atoms) {
        const auto [first, second] = pairs.atoms[index];
        score += atom_clash_penalty(output.getAtomPos(first),
                                    output.getAtomPos(second),
                                    pairs.minimum_distance);
    }
    for (const auto index : affected.bonds) {
        const auto [first, second] = pairs.bonds[index];
        if (bonds_cross(*first, *second, output)) {
            score += BOND_CROSSING_WEIGHT;
        }
    }
    for (const auto index : affected.disulfides) {
        score +=
            disulfide_length_penalty(*pairs.disulfides[index], pairs, output);
    }
    return score;
}

bool sidechain_needs_improvement(const AffectedPairs& affected,
                                 const ClashPairs& pairs,
                                 const RDKit::Conformer& output)
{
    // Try other positions if there is a clash, a bond crossing, or an S-S
    // bond longer than a normal small bond. This only starts a search;
    // longer S-S bonds remain allowed and receive a soft length penalty.
    for (const auto index : affected.atoms) {
        const auto [first, second] = pairs.atoms[index];
        if (atom_clash_penalty(output.getAtomPos(first),
                               output.getAtomPos(second),
                               pairs.minimum_distance) > 0) {
            return true;
        }
    }
    for (const auto index : affected.bonds) {
        const auto [first, second] = pairs.bonds[index];
        if (bonds_cross(*first, *second, output)) {
            return true;
        }
    }
    for (const auto index : affected.disulfides) {
        const auto* bond = pairs.disulfides[index];
        const auto delta = output.getAtomPos(bond->getBeginAtomIdx()) -
                           output.getAtomPos(bond->getEndAtomIdx());
        if (std::hypot(delta.x, delta.y) > pairs.sidechain_bond_length) {
            return true;
        }
    }
    return false;
}

struct AdjacentBackbone {
    std::vector<unsigned int> atoms;
    std::vector<const RDKit::Bond*> bonds;
};

std::map<ResidueKey, AdjacentBackbone>
index_adjacent_backbone(const RDKit::ROMol& mol, const ResidueMap& residues,
                        const AtomResidueMap& atom_residue,
                        const PeptideTopology& topology)
{
    // For each peptide residue, collect the heavy backbone atoms of its
    // immediate neighbors and the backbone bonds touching those neighbors.
    // Side-chain variants can then check nearby backbone clashes cheaply.
    std::map<ResidueKey, ResidueKeys> neighbors;
    for (const auto& [first, link] : topology.links) {
        neighbors[first].insert(link.next);
        neighbors[link.next].insert(first);
    }
    std::vector<const RDKit::Bond*> backbone_bonds;
    for (const auto* bond : mol.bonds()) {
        if (is_backbone_atom(bond->getBeginAtom()) &&
            is_backbone_atom(bond->getEndAtom())) {
            backbone_bonds.push_back(bond);
        }
    }

    std::map<ResidueKey, AdjacentBackbone> indexed;
    for (const auto& [residue_key, neighbor_keys] : neighbors) {
        auto& local = indexed[residue_key];
        for (const auto& neighbor_key : neighbor_keys) {
            for (const auto id : residues.at(neighbor_key).atom_ids) {
                const auto* atom = mol.getAtomWithIdx(id);
                if (atom->getAtomicNum() != 1 && is_backbone_atom(atom)) {
                    local.atoms.push_back(id);
                }
            }
        }
        for (const auto* bond : backbone_bonds) {
            if (neighbor_keys.contains(
                    atom_residue.at(bond->getBeginAtomIdx())) ||
                neighbor_keys.contains(
                    atom_residue.at(bond->getEndAtomIdx()))) {
                local.bonds.push_back(bond);
            }
        }
    }
    return indexed;
}

double adjacent_backbone_clash_score(const RDKit::ROMol& mol,
                                     const MovableSidechain& sidechain,
                                     const AdjacentBackbone& backbone,
                                     const ClashPairs& pairs,
                                     const RDKit::Conformer& output)
{
    // Score overlaps and bond crossings between this side chain and the
    // already indexed backbone of its peptide neighbors.
    double score = 0;
    for (const auto id : backbone.atoms) {
        for (const auto sidechain_id : sidechain.heavy_atoms) {
            if (mol.getBondBetweenAtoms(id, sidechain_id)) {
                continue;
            }
            score += atom_clash_penalty(output.getAtomPos(id),
                                        output.getAtomPos(sidechain_id),
                                        pairs.minimum_distance);
        }
    }
    for (const auto* moving_bond : sidechain.bonds) {
        for (const auto* fixed_bond : backbone.bonds) {
            if (bonds_cross(*moving_bond, *fixed_bond, output)) {
                score += BOND_CROSSING_WEIGHT;
            }
        }
    }
    return score;
}

Bounds2D sidechain_bounds(const MovableSidechain& sidechain,
                          const RDKit::Conformer& output)
{
    // Bound CA and every atom in the current side-chain variant. Include
    // hydrogens because pairwise clash checks also test their bonds.
    const auto alpha = output.getAtomPos(sidechain.alpha);
    Bounds2D bounds{alpha.x, alpha.x, alpha.y, alpha.y};
    for (const auto id : sidechain.atoms) {
        const auto point = output.getAtomPos(id);
        bounds.min_x = std::min(bounds.min_x, point.x);
        bounds.max_x = std::max(bounds.max_x, point.x);
        bounds.min_y = std::min(bounds.min_y, point.y);
        bounds.max_y = std::max(bounds.max_y, point.y);
    }
    return bounds;
}

double sidechain_pair_conflict(const RDKit::ROMol& mol,
                               const MovableSidechain& first,
                               const MovableSidechain& second,
                               const Bounds2D& first_bounds,
                               const Bounds2D& second_bounds,
                               const ClashPairs& pairs,
                               const RDKit::Conformer& output)
{
    // Score clashes between two current side-chain positions. Check atoms
    // and branch bonds only when their boxes are near, but always score an
    // S-S bond connecting the pair, including its crossings and length.
    double score = 0;
    if (bounds_near(first_bounds, second_bounds, pairs.minimum_distance)) {
        for (const auto first_id : first.heavy_atoms) {
            for (const auto second_id : second.heavy_atoms) {
                if (mol.getBondBetweenAtoms(first_id, second_id)) {
                    continue;
                }
                score += atom_clash_penalty(output.getAtomPos(first_id),
                                            output.getAtomPos(second_id),
                                            pairs.minimum_distance);
            }
        }
        for (const auto* first_bond : first.bonds) {
            for (const auto* second_bond : second.bonds) {
                if (bonds_cross(*first_bond, *second_bond, output)) {
                    score += BOND_CROSSING_WEIGHT;
                }
            }
        }
    }
    for (const auto* bond : pairs.disulfides) {
        const auto first_sulfur = bond->getBeginAtomIdx();
        const auto second_sulfur = bond->getEndAtomIdx();
        if (!((first.offsets.contains(first_sulfur) &&
               second.offsets.contains(second_sulfur)) ||
              (first.offsets.contains(second_sulfur) &&
               second.offsets.contains(first_sulfur)))) {
            continue;
        }
        score += disulfide_length_penalty(*bond, pairs, output);
        for (const auto& [first_bond, second_bond] : pairs.bonds) {
            if ((first_bond == bond || second_bond == bond) &&
                bonds_cross(*first_bond, *second_bond, output)) {
                score += BOND_CROSSING_WEIGHT;
            }
        }
    }
    return score;
}

struct SidechainVariantContext {
    const RDKit::ROMol& mol;
    const ClashPairs& pairs;
    const std::vector<AffectedPairs>& affected;
    const std::map<ResidueKey, AdjacentBackbone>& adjacent_backbone;
    RDKit::Conformer& output;

    double local_score(MovableSidechain& sidechain, size_t variant)
    {
        // Compute the adjacent-backbone score once per variant. The first
        // call applies that variant to output; callers restore the selected
        // coordinates after trying candidates.
        auto& score = sidechain.local_scores[variant];
        if (!score) {
            apply_sidechain_variant(sidechain, variant, output);
            const auto neighbor = adjacent_backbone.find(sidechain.residue);
            score = neighbor == adjacent_backbone.end()
                        ? 0
                        : adjacent_backbone_clash_score(
                              mol, sidechain, neighbor->second, pairs, output);
        }
        return *score;
    }

    bool available(MovableSidechain& sidechain, size_t variant)
    {
        // Normally allow only variants clear of adjacent backbone clashes.
        // If none are clear, later passes may try every variant.
        return sidechain.allow_all_variants ||
               local_score(sidechain, variant) == 0;
    }
};

void choose_variants_clear_of_adjacent_backbone(
    std::vector<MovableSidechain>& sidechains, SidechainVariantContext& context)
{
    // Try the original and three alternatives in order, choosing the first
    // one clear of neighboring backbone atoms and bonds. If none are clear,
    // keep the lowest local score and let later passes try all variants.
    for (auto& sidechain : sidechains) {
        bool has_clear_variant = false;
        size_t selected = 0;
        for (size_t variant = 0; variant < SIDECHAIN_VARIANTS; ++variant) {
            if (context.available(sidechain, variant)) {
                selected = variant;
                has_clear_variant = true;
                break;
            }
        }
        if (!has_clear_variant) {
            sidechain.allow_all_variants = true;
            for (size_t variant = 1; variant < SIDECHAIN_VARIANTS; ++variant) {
                if (*sidechain.local_scores[variant] <
                    *sidechain.local_scores[selected]) {
                    selected = variant;
                }
            }
        }
        sidechain.selected = selected;
        apply_sidechain_variant(sidechain, selected, context.output);
    }
}

void improve_single_sidechains(std::vector<MovableSidechain>& sidechains,
                               SidechainVariantContext& context)
{
    // Before pair trials, try each side chain alone if it clashes with a
    // fixed group, crosses a bond, or has a long S-S bridge. Keep a variant
    // only when its combined clash score improves.
    constexpr double improvement_epsilon = 1e-6;
    for (size_t i = 0; i < sidechains.size(); ++i) {
        auto& sidechain = sidechains[i];
        const auto& affected = context.affected[i];
        if (!sidechain_needs_improvement(affected, context.pairs,
                                         context.output)) {
            continue;
        }
        auto best_variant = sidechain.selected;
        auto best_score = clash_score(affected, context.pairs, context.output);
        for (size_t variant = 0; variant < SIDECHAIN_VARIANTS; ++variant) {
            if (!context.available(sidechain, variant) ||
                variant == sidechain.selected) {
                continue;
            }
            apply_sidechain_variant(sidechain, variant, context.output);
            const auto trial_score =
                clash_score(affected, context.pairs, context.output);
            if (trial_score < best_score - improvement_epsilon) {
                best_score = trial_score;
                best_variant = variant;
            }
        }
        sidechain.selected = best_variant;
        apply_sidechain_variant(sidechain, best_variant, context.output);
    }
}

void improve_clashing_sidechain_pairs(const RDKit::ROMol& mol,
                                      std::vector<MovableSidechain>& sidechains,
                                      SidechainVariantContext& context)
{
    // Rank side-chain pairs by current clash and S-S severity. For each
    // pair, try up to 4 x 4 available positions and keep a lower total
    // score. Limit passes, pair trials, and accepted moves per residue.
    constexpr size_t max_passes = 2;
    constexpr size_t max_pair_trials_per_residue = 4;
    constexpr size_t max_moves_per_residue = 2;
    constexpr double improvement_epsilon = 1e-6;
    std::vector<size_t> pair_trials(sidechains.size(), 0);
    std::vector<size_t> moves(sidechains.size(), 0);
    for (size_t pass = 0; pass < max_passes; ++pass) {
        // Find pairs worth improving, including distant Cys side chains
        // connected by a stretched S-S bond.
        std::vector<Bounds2D> bounds;
        bounds.reserve(sidechains.size());
        for (const auto& sidechain : sidechains) {
            bounds.push_back(sidechain_bounds(sidechain, context.output));
        }
        std::vector<std::tuple<double, size_t, size_t>> conflicts;
        for (size_t i = 0; i < sidechains.size(); ++i) {
            for (size_t j = i + 1; j < sidechains.size(); ++j) {
                const auto severity = sidechain_pair_conflict(
                    mol, sidechains[i], sidechains[j], bounds[i], bounds[j],
                    context.pairs, context.output);
                if (severity > improvement_epsilon) {
                    conflicts.emplace_back(severity, i, j);
                }
            }
        }
        std::ranges::sort(conflicts, std::greater{});
        bool improved = false;
        for (const auto& [severity, i, j] : conflicts) {
            (void)severity;
            if (pair_trials[i] >= max_pair_trials_per_residue ||
                pair_trials[j] >= max_pair_trials_per_residue) {
                continue;
            }
            ++pair_trials[i];
            ++pair_trials[j];
            auto& first = sidechains[i];
            auto& second = sidechains[j];
            const auto affected = combine_affected_pairs(context.affected[i],
                                                         context.affected[j]);
            const auto original_first = first.selected;
            const auto original_second = second.selected;
            auto best_first = original_first;
            auto best_second = original_second;
            const auto original_score =
                clash_score(affected, context.pairs, context.output);
            auto best_score = original_score;
            // Score only the pairs affected by these two side chains.
            for (size_t first_variant = 0; first_variant < SIDECHAIN_VARIANTS;
                 ++first_variant) {
                if (!context.available(first, first_variant) ||
                    (moves[i] >= max_moves_per_residue &&
                     first_variant != original_first)) {
                    continue;
                }
                apply_sidechain_variant(first, first_variant, context.output);
                for (size_t second_variant = 0;
                     second_variant < SIDECHAIN_VARIANTS; ++second_variant) {
                    if (!context.available(second, second_variant) ||
                        (moves[j] >= max_moves_per_residue &&
                         second_variant != original_second)) {
                        continue;
                    }
                    apply_sidechain_variant(second, second_variant,
                                            context.output);
                    const auto trial_score =
                        clash_score(affected, context.pairs, context.output);
                    if (trial_score < best_score - improvement_epsilon) {
                        best_score = trial_score;
                        best_first = first_variant;
                        best_second = second_variant;
                    }
                }
            }
            apply_sidechain_variant(first, best_first, context.output);
            apply_sidechain_variant(second, best_second, context.output);
            if (best_score < original_score - improvement_epsilon) {
                moves[i] += best_first != original_first;
                moves[j] += best_second != original_second;
                first.selected = best_first;
                second.selected = best_second;
                improved = true;
            }
        }
        if (!improved) {
            break;
        }
    }
}

void resolve_sidechain_clashes(const RDKit::ROMol& mol,
                               const ResidueMap& residues,
                               const AtomResidueMap& atom_residue,
                               const PeptideTopology& topology,
                               const ResidueKeys& placed_peptides,
                               RDKit::Conformer& output)
{
    // Collect movable branches and precompute clash indexes, then improve
    // their positions against adjacent backbone, single-side-chain, and
    // pairwise conflicts in that order.
    std::vector<MovableSidechain> sidechains;
    for (const auto& key : placed_peptides) {
        if (auto sidechain =
                make_movable_sidechain(mol, key, residues.at(key), output)) {
            sidechains.push_back(std::move(*sidechain));
        }
    }
    if (sidechains.empty()) {
        return;
    }
    const auto pairs = make_clash_pairs(mol, sidechains, atom_residue, output);
    const auto affected = index_affected_pairs(mol, sidechains, pairs);
    const auto adjacent_backbone =
        index_adjacent_backbone(mol, residues, atom_residue, topology);
    SidechainVariantContext context{mol, pairs, affected, adjacent_backbone,
                                    output};
    choose_variants_clear_of_adjacent_backbone(sidechains, context);
    improve_single_sidechains(sidechains, context);
    improve_clashing_sidechain_pairs(mol, sidechains, context);
}

ResidueKeys place_peptide_chains_with_best_phase(
    const RDKit::ROMol& mol, const RDKit::Conformer& monomer_conf,
    const ResidueMap& residues, const PeptideTopology& topology,
    const DisulfidePartners& disulfide_partners, double layout_scale,
    std::unique_ptr<RDKit::Conformer>& output)
{
    // Choose each chain's starting side-chain phase to minimize non-proline
    // LINEAR_CIS_N_CA shapes. Draw both phases only when the default needs
    // such a shape, and redraw a mixed phase only if chains disagree.
    // First draw the default phase and count those shape changes.
    LinearCounts default_counts;
    auto placed_peptides = place_peptide_chains(
        mol, monomer_conf, residues, topology, disulfide_partners, layout_scale,
        {}, default_counts, *output);
    const auto has_linear_cis = std::ranges::any_of(
        default_counts, [](const auto& entry) { return entry.second > 0; });
    if (!has_linear_cis) {
        return placed_peptides;
    }

    // Draw the opposite starting phase for every chain. Each chain may
    // choose its own phase; a tie keeps the default.
    ResidueKeys all_chain_starts;
    for (const auto& entry : default_counts) {
        all_chain_starts.insert(entry.first);
    }
    auto flipped_output = std::make_unique<RDKit::Conformer>(mol.getNumAtoms());
    flipped_output->set3D(false);
    LinearCounts flipped_counts;
    auto flipped_peptides = place_peptide_chains(
        mol, monomer_conf, residues, topology, disulfide_partners, layout_scale,
        all_chain_starts, flipped_counts, *flipped_output);
    ResidueKeys selected_flips;
    for (const auto& [start, count] : default_counts) {
        if (flipped_counts.at(start) < count) {
            selected_flips.insert(start);
        }
    }

    // Reuse an already drawn result when all chains agree. Otherwise draw
    // the mixture of default and flipped chains once.
    if (selected_flips == all_chain_starts) {
        output = std::move(flipped_output);
        return flipped_peptides;
    }
    if (!selected_flips.empty()) {
        output = std::make_unique<RDKit::Conformer>(mol.getNumAtoms());
        output->set3D(false);
        LinearCounts selected_counts;
        return place_peptide_chains(mol, monomer_conf, residues, topology,
                                    disulfide_partners, layout_scale,
                                    selected_flips, selected_counts, *output);
    }
    return placed_peptides;
}

unsigned int compute_atomistic_coords_using_monomer_layout(RDKit::ROMol& mol)
{
    // note: this assumes the polymer is peptidic. support for non-peptidic
    // chains is not implemented yet.
    // Build the monomer layout that will guide the atomistic coordinates.
    auto residues = get_residue_geometry(mol);
    if (residues.empty()) {
        throw std::invalid_argument("Monomer layout requires residue atoms");
    }

    auto monomer_mol = toMonomeric(mol);
    // This also creates the conformer read below.
    const auto layout_scale = prepare_monomer_layout(*monomer_mol, residues);

    // Find the peptide path and any disulfide constraints, then place the
    // peptide backbone using the starting phase best suited to each chain.
    const auto& monomer_conf = monomer_mol->getConformer();
    const auto atom_residue = map_atoms_to_residues(residues);
    const auto topology = find_peptide_links(mol, *monomer_mol, monomer_conf,
                                             residues, atom_residue);
    const auto disulfide_partners = find_disulfide_partners(mol, atom_residue);
    auto output = std::make_unique<RDKit::Conformer>(mol.getNumAtoms());
    output->set3D(false);
    const auto placed_peptides = place_peptide_chains_with_best_phase(
        mol, monomer_conf, residues, topology, disulfide_partners, layout_scale,
        output);

    // Place everything outside the peptide backbone, correct explicit
    // backbone hydrogens, and resolve side-chain clashes last.
    place_remaining_residues(mol, *monomer_mol, monomer_conf, residues,
                             atom_residue, placed_peptides, layout_scale,
                             *output);
    place_backbone_hydrogens(mol, residues, placed_peptides, *output);
    resolve_sidechain_clashes(mol, residues, atom_residue, topology,
                              placed_peptides, *output);

    // Replace the old coordinates with the completed 2D conformer.
    clear_molfile_bond_stereo(mol);
    mol.clearConformers();
    return mol.addConformer(output.release(), true);
}

} // namespace

RDGeom::Point3D compute_centroid(const std::vector<RDGeom::Point3D>& positions)
{
    // Calculate the centroid by averaging the coordinates
    size_t num_atoms = positions.size();
    RDGeom::Point3D centroid;
    for (const auto& coord : positions) {
        centroid += coord;
    }
    // avoid division by zero
    if (num_atoms > 0) {
        centroid /= static_cast<double>(num_atoms);
    }
    return centroid;
}

unsigned int compute2DCoords(RDKit::ROMol& mol,
                             const std::vector<unsigned int>& frozen_ids,
                             bool use_monomer_layout)
{
    // Monomeric input: lay out the monomer beads using monomeric coordgen.
    if (isMonomeric(mol)) {
        return compute_monomer_mol_coords(mol);
    }
    // Atomistic input: first lay out its monomers using monomeric coordgen,
    // then lay out its atoms to follow that monomer layout.
    if (use_monomer_layout && frozen_ids.empty() &&
        has_complete_residue_info(mol)) {
        return compute_atomistic_coords_using_monomer_layout(mol);
    }

    // Atomistic input: let RDKit lay it out, honoring any frozen atoms.
    RDDepict::Compute2DCoordParameters params;
    params.forceRDKit = true;
    params.useRingTemplates = true;

    RDGeom::INT_POINT2D_MAP coord_map;
    if (!frozen_ids.empty()) {
        auto& conf = mol.getConformer();
        for (auto idx : frozen_ids) {
            coord_map[idx] = conf.getAtomPos(idx);
        }
        params.coordMap = &coord_map;
    }

    clear_molfile_bond_stereo(mol);

    return RDDepict::compute2DCoords(mol, params);
}

} // namespace rdkit_extensions
} // namespace schrodinger
