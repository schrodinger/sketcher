/* -------------------------------------------------------------------------
 * Implements schrodinger::rdkit_extensions:: text block <-> rdkit
 mol conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */
#include "schrodinger/rdkit_extensions/convert.h"

#include <functional>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/inchi.h>

#include "schrodinger/rdkit_extensions/capture_rdkit_log.h"

namespace schrodinger
{
namespace rdkit_extensions
{

namespace
{
const std::string attachment_point_label_prefix{"_AP"};

/**
 * @internal
 * Guess by attempting conversions, and returning the first valid parse
 */
template <typename T>
boost::shared_ptr<T> auto_detect(
    const std::string& text, const std::vector<Format> formats,
    std::function<boost::shared_ptr<T>(const std::string&, Format)> text_to_)
{
    for (const auto& fmt : formats) {
        try {
            return text_to_(text, fmt);
        } catch (const std::exception&) {
            continue;
        }
    }
    throw std::invalid_argument("Unable to determine format");
}

void attachment_point_dummies_to_molattachpt_property(RDKit::RWMol& rdk_mol)
{
    rdk_mol.updatePropertyCache(false);
    std::vector<RDKit::Atom*> dummies_to_remove;
    for (auto atom : rdk_mol.atoms()) {
        if (is_attachment_point_dummy(*atom)) {
            dummies_to_remove.push_back(atom);

            // This should find only one neighbor
            for (auto parent : rdk_mol.atomNeighbors(atom)) {
                int attachment_point_type{1};
                auto num_exp_hs = parent->getNumExplicitHs();
                if (parent->hasProp(RDKit::common_properties::molAttachPoint)) {
                    attachment_point_type = -1;
                    parent->setNumExplicitHs(num_exp_hs + 2);
                } else {
                    parent->setNumExplicitHs(num_exp_hs + 1);
                }
                parent->setProp(RDKit::common_properties::molAttachPoint,
                                attachment_point_type);
            }
        }
    }
    rdk_mol.beginBatchEdit();
    for (auto atom : dummies_to_remove) {
        rdk_mol.removeAtom(atom);
    }
    rdk_mol.commitBatchEdit();

    rdk_mol.updatePropertyCache(false);
}

bool molattachpt_property_to_attachment_point_dummies(RDKit::RWMol& rdk_mol)
{
    std::vector<std::pair<int, int>> new_attachment_point_bonds;
    for (auto& atom : rdk_mol.atoms()) {
        int attach_point_type{0};
        if (atom->getPropIfPresent(RDKit::common_properties::molAttachPoint,
                                   attach_point_type)) {
            auto num_dummies = (attach_point_type == -1 ? 2 : 1);
            auto explicit_h_count = atom->getNumExplicitHs();
            for (auto i = 1; i <= num_dummies; ++i) {

                auto dummy_atom = new RDKit::QueryAtom(0);
                dummy_atom->setQuery(RDKit::makeAtomNullQuery());

                dummy_atom->setProp(RDKit::common_properties::atomLabel,
                                    attachment_point_label_prefix +
                                        std::to_string(i));

                bool update_label = false;
                bool take_ownership = true;
                auto dummy_idx =
                    rdk_mol.addAtom(dummy_atom, update_label, take_ownership);
                rdk_mol.addBond(atom, dummy_atom, RDKit::Bond::SINGLE);

                new_attachment_point_bonds.push_back(
                    {dummy_idx, atom->getIdx()});
            }
            if (explicit_h_count != 0) {
                atom->setNumExplicitHs(explicit_h_count - num_dummies);
            }
            atom->clearProp(RDKit::common_properties::molAttachPoint);
        }
    }

    if (!new_attachment_point_bonds.empty()) {
        // set coordinates for new attachment point atoms after adding all of
        // them
        // so that all new atoms are considered in coordinate generation
        for (const auto& new_bond : new_attachment_point_bonds) {
            RDKit::MolOps::setTerminalAtomCoords(rdk_mol, new_bond.first,
                                                 new_bond.second);
        }
        reapply_molblock_wedging(rdk_mol);
        RDKit::MolOps::assignChiralTypesFromBondDirs(rdk_mol);
        rdk_mol.updatePropertyCache(false);
        return true;
    }
    return false;
}

void preserve_wiggly_bonds(RDKit::RWMol& mol)
{
    for (auto& bond : mol.bonds()) {
        int wiggly_bond_v2000{0};
        int wiggly_bond_v3000{0};
        bond->getPropIfPresent(RDKit::common_properties::_MolFileBondStereo,
                               wiggly_bond_v2000);
        bond->getPropIfPresent(RDKit::common_properties::_MolFileBondCfg,
                               wiggly_bond_v3000);
        if (wiggly_bond_v2000 == 4 || wiggly_bond_v3000 == 2) {
            bond->setBondDir(RDKit::Bond::BondDir::UNKNOWN);
        }
    }
}

void fix_cxsmiles_rgroups(RDKit::ROMol& mol)
{
    for (auto atom : mol.atoms()) {
        if (atom->getAtomicNum() == 0) {

            // In SMILES RGroups have no query, but SMARTS & MDL
            // they should have an empty query
            if (atom->hasQuery()) {
                auto queryAtom = static_cast<RDKit::QueryAtom*>(atom);
                auto query_label = queryAtom->getQuery()->getTypeLabel();
                if (!query_label.empty()) {
                    continue;
                }
            }

            // If this is a null query, check for a proper RGroup atom label
            std::string atom_label;
            if (atom->getPropIfPresent(RDKit::common_properties::atomLabel,
                                       atom_label) &&
                atom_label.find("_R") == 0) {
                int idx = -1;
                size_t end = 0;
                try {
                    idx = std::stoi(atom_label.substr(2), &end);
                } catch (const std::invalid_argument&) {
                    continue;
                }
                if (idx > -1 && end == atom_label.size() - 2) {
                    // RDKit::setAtomRLabel() wouldn't work for R0
                    // (this will be corrected lated on)
                    atom->setProp(RDKit::common_properties::_MolFileRLabel,
                                  static_cast<unsigned>(idx));
                }
            }
        }
    }
}

void fix_r0_rgroup(RDKit::ROMol& mol)
{
    int highest_rlabel{0};
    std::vector<RDKit::Atom*> atoms_with_r0;
    for (auto atom : mol.atoms()) {
        int rlabel{-1};
        if (atom->getPropIfPresent(RDKit::common_properties::_MolFileRLabel,
                                   rlabel)) {
            highest_rlabel = std::max(highest_rlabel, rlabel);
            if (rlabel == 0) {
                atoms_with_r0.push_back(atom);
            }
        }
    }
    if (!atoms_with_r0.empty()) {
        ++highest_rlabel;
        for (auto atom : atoms_with_r0) {
            RDKit::setAtomRLabel(atom, highest_rlabel);
        }
    }
}

// For mol files, RDKit reads the molTotValence property from the atom block,
// and adjusts Hydrogens on the atom to match this valence, but still leaves the
// property in place. We should remove it, since it will make no sense if we
// modify
// the mol.
void cleanup_moltotvalence_property(RDKit::ROMol& mol)
{
    for (auto atom : mol.atoms()) {
        atom->clearProp(RDKit::common_properties::molTotValence);
    }
}

boost::shared_ptr<RDKit::ROMol>
molblock_to_romol(const std::string& text,
                  const CaptureRDErrorLog& rd_error_log, bool sanitize,
                  bool removeHs)
{
    bool strictParsing = false;
    // use SDMolSupplier because MolFromMolBlock does not preserve
    // structure level properties
    RDKit::SDMolSupplier reader;
    reader.setData(text, sanitize, removeHs, strictParsing);
    RDKit::ROMol* romol;
    try {
        romol = reader.next();
    } catch (std::runtime_error&) { // EOF hit
        throw std::invalid_argument("Failed to parse text: " +
                                    rd_error_log.messages());
    }
    if (romol == nullptr) {
        throw std::invalid_argument("Failed to parse text: " +
                                    rd_error_log.messages());
    }
    if (!reader.atEnd()) {
        throw std::invalid_argument(
            "Single molblock required; multiple present" +
            rd_error_log.messages());
    }
    return boost::shared_ptr<RDKit::ROMol>(romol);
}

void adjust_for_mdl_v2k_format(RDKit::RWMol& mol)
{
    auto nsgroups = getSubstanceGroups(mol).size();
    if (mol.getNumAtoms() > 999u || mol.getNumBonds() > 999u ||
        nsgroups > 999u) {
        throw std::runtime_error("Mol incompatible with MDL V2000 format");
    }
    // MDL V2000 format does not support enhanced stereo groups; if RDKit
    // detects then on write, it will force V3000 format to prevent
    // information loss; here, we explicitly clear enhanced stereo so that
    // V2000 can be written as requested.
    mol.setStereoGroups({});
}

} // unnamed namespace

boost::shared_ptr<RDKit::RWMol> text_to_rdmol(const std::string& text,
                                              Format format)
{
    if (format == Format::AUTO_DETECT) {
        // NOTE: Attempt SMILES before SMARTS, given not all SMARTS are SMILES
        return auto_detect<RDKit::RWMol>(text, {Format::MDL_MOLV3000,
                                                Format::PDB, Format::INCHI,
                                                Format::SMILES, Format::SMARTS},
                                         &text_to_rdmol);
    }

    CaptureRDErrorLog rd_error_log;

    RDKit::RWMol* mol = nullptr;
    bool sanitize = false;
    bool removeHs = false;
    switch (format) {
        case Format::SMILES:
        case Format::EXTENDED_SMILES: {
            int debugParse = 0;
            mol = RDKit::SmilesToMol(text, debugParse, sanitize);
            if (mol != nullptr) {
                fix_cxsmiles_rgroups(*mol);
            }
            break;
        }
        case Format::SMARTS:
            mol = RDKit::SmartsToMol(text);
            if (mol != nullptr) {
                fix_cxsmiles_rgroups(*mol);
            }
            break;
        case Format::MDL_MOLV2000:
        case Format::MDL_MOLV3000: {
            auto romol =
                molblock_to_romol(text, rd_error_log, sanitize, removeHs);
            mol = new RDKit::RWMol(*romol);
            break;
        }
        case Format::INCHI: {
            RDKit::ExtraInchiReturnValues aux;
            mol = RDKit::InchiToMol(text, aux, sanitize, removeHs);
            break;
        }
        case Format::INCHI_KEY:
            throw std::invalid_argument("Cannot read from INCHI_KEY");
        case Format::PDB:
            mol = RDKit::PDBBlockToMol(text, sanitize, removeHs);
            break;
        default:
            throw std::invalid_argument("Unsupported import format");
    }

    if (mol == nullptr) {
        throw std::invalid_argument("Failed to parse text: " +
                                    rd_error_log.messages());
    }

    bool has_attchpt = molattachpt_property_to_attachment_point_dummies(*mol);
    if (!has_attchpt) {
        preserve_wiggly_bonds(*mol);
        fix_r0_rgroup(*mol);
        mol->updatePropertyCache(false);
    }

    cleanup_moltotvalence_property(*mol);

    return boost::shared_ptr<RDKit::RWMol>(mol);
}

boost::shared_ptr<RDKit::ChemicalReaction>
text_to_reaction(const std::string& text, Format format)
{
    if (format == Format::AUTO_DETECT) {
        // NOTE: Text is always read as SMARTS. Reading text as
        // SMILES should be explicitly requested
        return auto_detect<RDKit::ChemicalReaction>(
            text, {Format::MDL_MOLV2000, Format::SMARTS}, &text_to_reaction);
    }

    CaptureRDErrorLog rd_error_log;

    RDKit::ChemicalReaction* reaction = nullptr;
    switch (format) {
        case Format::SMILES:
        case Format::SMARTS: {
            auto useSMILES = (format == Format::SMILES);
            reaction =
                RDKit::RxnSmartsToChemicalReaction(text, nullptr, useSMILES);

            // Bootstrap coords to preserve double bond stereo (SKETCH-1254)
            double spacing = 2.0;
            bool updateProps = true;
            RDDepict::compute2DCoordsForReaction(*reaction, spacing,
                                                 updateProps);

            break;
        }
        case Format::MDL_MOLV2000:
        case Format::MDL_MOLV3000: {
            bool sanitize = false;
            bool removeHs = false;
            reaction =
                RDKit::RxnBlockToChemicalReaction(text, sanitize, removeHs);
            break;
        }
        default:
            throw std::invalid_argument("Unsupported reaction import format");
    }
    if (reaction == nullptr) {
        throw std::invalid_argument("Failed to parse text: " +
                                    rd_error_log.messages());
    }

    for (auto reactant : reaction->getReactants()) {
        fix_r0_rgroup(*reactant);
    }

    for (auto product : reaction->getProducts()) {
        fix_r0_rgroup(*product);
    }

    return boost::shared_ptr<RDKit::ChemicalReaction>(reaction);
}

std::string rdmol_to_text(const RDKit::ROMol& mol, Format format)
{
    RDKit::RWMol rwmol(mol);

    bool kekulize = false;
    bool include_stereo = true;

    switch (format) {
        case Format::SMILES:
            return RDKit::MolToSmiles(mol, include_stereo, kekulize);
        case Format::EXTENDED_SMILES:
            rwmol.clearConformers(); // Don't pollute extension with coordinates
            return RDKit::MolToCXSmiles(rwmol, include_stereo, kekulize);
        case Format::SMARTS:
            return RDKit::MolToSmarts(mol);
        case Format::MDL_MOLV2000:
        case Format::MDL_MOLV3000: {
            attachment_point_dummies_to_molattachpt_property(rwmol);
            if (format == Format::MDL_MOLV2000) {
                adjust_for_mdl_v2k_format(rwmol);
            }
            int confId = -1;
            bool forceV3000 = format == Format::MDL_MOLV3000;
            return RDKit::SDWriter::getText(rwmol, confId, kekulize,
                                            forceV3000);
        }
        case Format::INCHI: {
            RDKit::ExtraInchiReturnValues aux;
            return RDKit::MolToInchi(mol, aux);
        }
        case Format::INCHI_KEY:
            return RDKit::MolToInchiKey(mol);
        case Format::PDB:
            return RDKit::MolToPDBBlock(mol);
        default:
            throw std::invalid_argument("Unsupported export format");
    }
    throw std::invalid_argument("Invalid format specified");
}

std::string reaction_to_text(const RDKit::ChemicalReaction& reaction,
                             Format format)
{
    switch (format) {
        case Format::SMILES:
            return RDKit::ChemicalReactionToRxnSmiles(reaction);
        case Format::SMARTS:
            return RDKit::ChemicalReactionToRxnSmarts(reaction);
        case Format::MDL_MOLV2000:
        case Format::MDL_MOLV3000: {
            bool separateAgents = false;
            bool forceV3000 = format == Format::MDL_MOLV3000;
            return RDKit::ChemicalReactionToRxnBlock(reaction, separateAgents,
                                                     forceV3000);
        }
        default:
            throw std::invalid_argument("Unsupported reaction export format");
    }
    throw std::invalid_argument("Invalid format specified");
}

bool is_attachment_point_dummy(const RDKit::Atom& atom)
{
    std::string label;
    return atom.getAtomicNum() == 0 && atom.getTotalDegree() == 1 &&
           atom.getPropIfPresent(RDKit::common_properties::atomLabel, label) &&
           label.find(attachment_point_label_prefix) == 0;
}

void add_enhanced_stereo_to_chiral_atoms(RDKit::ROMol& mol)
{
    auto stereo_groups = mol.getStereoGroups();
    std::unordered_set<RDKit::Atom*> seen_chiral_atoms;
    for (auto sg : stereo_groups) {
        auto atoms = sg.getAtoms();
        seen_chiral_atoms.insert(atoms.begin(), atoms.end());
    }

    std::vector<RDKit::Atom*> ungrouped_atoms;
    for (auto atom : mol.atoms()) {
        if ((atom->getChiralTag() == RDKit::Atom::CHI_TETRAHEDRAL_CW ||
             atom->getChiralTag() == RDKit::Atom::CHI_TETRAHEDRAL_CCW) &&
            seen_chiral_atoms.count(atom) == 0) {
            ungrouped_atoms.push_back(atom);
        }
    }

    if (ungrouped_atoms.empty()) {
        return;
    }

    int chiral_flag{0};
    mol.getPropIfPresent(RDKit::common_properties::_MolFileChiralFlag,
                         chiral_flag);

    if (chiral_flag == 1) {
        // We have a chiral flag: append to or create the ABS group
        // (ABS is unique!)

        auto is_abs = [](const auto& sg) {
            return sg.getGroupType() == RDKit::StereoGroupType::STEREO_ABSOLUTE;
        };
        auto abs_group =
            std::find_if(stereo_groups.begin(), stereo_groups.end(), is_abs);

        if (abs_group != stereo_groups.end()) {
            auto abs_atoms = abs_group->getAtoms();
            ungrouped_atoms.insert(ungrouped_atoms.end(), abs_atoms.begin(),
                                   abs_atoms.end());
            stereo_groups.erase(abs_group);
        }

        stereo_groups.emplace_back(RDKit::StereoGroupType::STEREO_ABSOLUTE,
                                   std::move(ungrouped_atoms));
    } else {
        // No chiral flag: ungrouped atoms go into a new AND group
        stereo_groups.emplace_back(RDKit::StereoGroupType::STEREO_AND,
                                   std::move(ungrouped_atoms));
    }

    mol.setStereoGroups(stereo_groups);
}

void reapply_molblock_wedging(RDKit::ROMol& rdk_mol)
{
    for (auto& bond : rdk_mol.bonds()) {
        // only change the wedging if the bond was wedged in the input data (we
        // recognize that by looking for the properties the mol file parser
        // sets):
        unsigned int bond_dir_val{0};
        if (bond->getPropIfPresent(RDKit::common_properties::_MolFileBondStereo,
                                   bond_dir_val)) {
            // v2000
            if (bond_dir_val == 1) {
                bond->setBondDir(RDKit::Bond::BondDir::BEGINWEDGE);
            } else if (bond_dir_val == 6) {
                bond->setBondDir(RDKit::Bond::BondDir::BEGINDASH);
            } else if (bond_dir_val == 4) {
                bond->setBondDir(RDKit::Bond::BondDir::UNKNOWN);
            }
        } else if (bond->getPropIfPresent(
                       RDKit::common_properties::_MolFileBondCfg,
                       bond_dir_val)) {
            // v3000
            if (bond_dir_val == 1) {
                bond->setBondDir(RDKit::Bond::BondDir::BEGINWEDGE);
            } else if (bond_dir_val == 3) {
                bond->setBondDir(RDKit::Bond::BondDir::BEGINDASH);
            } else if (bond_dir_val == 2) {
                bond->setBondDir(RDKit::Bond::BondDir::UNKNOWN);
            }
        } else {
            auto bond_dir = bond->getBondDir();
            if (bond_dir == RDKit::Bond::BondDir::BEGINWEDGE ||
                bond_dir == RDKit::Bond::BondDir::BEGINDASH ||
                bond_dir == RDKit::Bond::BondDir::UNKNOWN) {
                bond->setBondDir(RDKit::Bond::BondDir::UNKNOWN);
            }
        }
    }
}

} // namespace rdkit_extensions
} // namespace schrodinger
