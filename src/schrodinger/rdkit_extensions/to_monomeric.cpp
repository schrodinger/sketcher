#include "schrodinger/rdkit_extensions/atomistic_conversions.h"

#include <chrono>
#include <queue>
#include <span>
#include <unordered_map>
#include <memory>

#include <fmt/format.h>

#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include <rdkit/GraphMol/SubstanceGroup.h>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/rdkit_extensions/sup_utils.h"

namespace schrodinger
{
namespace rdkit_extensions
{

namespace
{

static constexpr const char* MONOMER_IDX{"monomerIndex"};
static constexpr const char* MONOMER_MAP_NUM{"monomerMapNumber"};
static constexpr const char* REFERENCE_IDX{"referenceIndex"};
static constexpr const char* ATTACH_NUM{
    "attachNumber"}; // Whether an atom corresponds to R1, R2, or R3

static constexpr int MIN_ATTCHPTS = 2;
static constexpr int SIDECHAIN_ATTCHPT = 8;
static constexpr int TERMINAL_ATTCHPT = 9;
static constexpr auto NO_ATTACHMENT = std::numeric_limits<unsigned int>::max();

// std::map to allow sequential/ordered iteration
using ChainsAndResidues =
    std::map<std::string,
             std::map<std::pair<int, std::string>, std::vector<unsigned int>>>;

struct ResidueQuery {
    std::unique_ptr<RDKit::RWMol> mol;
    std::vector<unsigned int> attch_map;
};

ResidueQuery prepare_static_mol_query(const char* smarts_query)
{
    ResidueQuery query;
    query.mol.reset(RDKit::SmartsToMol(smarts_query));

    // Maps SMARTS query index to attachment point # to avoid checking the
    // property for every match
    query.attch_map.resize(query.mol->getNumAtoms(), NO_ATTACHMENT);

    for (const auto atom : query.mol->atoms()) {
        if (atom->hasProp(RDKit::common_properties::molAtomMapNumber)) {
            query.attch_map[atom->getIdx()] = atom->getProp<unsigned int>(
                RDKit::common_properties::molAtomMapNumber);
        }
    }
    return query;
};

// clang-format off
// attachment points 1 and 2 are backbone attachment points,
// 8 is the side chain attachment point, 3 is cysteine's sulfur
static const ResidueQuery CYSTEINE_QUERY{prepare_static_mol_query(
    "[NX3,NX4+:1][CX4H]([CX4H2][S:3])[CX3:2](=[OX1])[O,N:9]")}; // matches C, dC, meC
static const ResidueQuery GENERIC_AMINO_ACID_QUERY{prepare_static_mol_query(
    "[NX3,NX4+:1][CX4H]([*:8])[CX3:2](=[OX1])[O,N:9]")};
static const ResidueQuery GLYCINE_AMINO_ACID_QUERY{prepare_static_mol_query(
    "[NX3,NX4+:1][CX4H2][CX3:2](=[OX1])[O,N:9]")}; // no side chain
// clang-format on

static const std::unordered_map<std::string, ChainType> BIOVIA_CHAIN_TYPE_MAP =
    {{"AA", ChainType::PEPTIDE}};

struct MonomerMatch {
    std::vector<unsigned int> atom_indices = {}; // atom indices of the match
    unsigned int r1 = NO_ATTACHMENT;             // attachment point 1
    unsigned int r2 = NO_ATTACHMENT;             // attachment point 2
    unsigned int r3 = NO_ATTACHMENT;             // attachment point 3
};

struct Linkage {
    unsigned int monomer_idx1;
    unsigned int monomer_idx2;
    unsigned int attach_from;
    unsigned int attach_to;

    Linkage(unsigned int idx1, unsigned int idx2, unsigned int from,
            unsigned int to) :
        monomer_idx1(idx1),
        monomer_idx2(idx2),
        attach_from(from),
        attach_to(to)
    {
    }

    std::string to_string() const
    {
        return fmt::format("R{}-R{}", attach_from, attach_to);
    }
};

bool alreadyMatched(const RDKit::ROMol& mol, std::span<const unsigned int> ids)
{
    // Make sure this match hasn't already been accounted for by a previous
    // match
    for (auto id : ids) {
        auto at = mol.getAtomWithIdx(id);
        unsigned int attch;
        if (at->hasProp(MONOMER_IDX) &&
            at->getPropIfPresent(ATTACH_NUM, attch) && attch == NO_ATTACHMENT) {
            return false;
        }
    }
    return true;
}

/*
 * Function that takes the SMARTS query and atomistic molecule and adds the atom
 * indices of the matches to the monomers vector
 *
 */
void addMatchesToMonomers(
    const ResidueQuery& query, RDKit::ROMol& atomistic_mol,
    std::vector<MonomerMatch>& monomers,
    std::unordered_map<unsigned int, unsigned int>& sidechain_attch_pts)
{
    // Set a final function check that ensures the entire match has not
    // already been accounted for by a previous SMARTS search.
    RDKit::SubstructMatchParameters params;
    params.useChirality = false;
    params.extraFinalCheck = &alreadyMatched;
    auto matches = RDKit::SubstructMatch(atomistic_mol, *query.mol, params);

    for (const auto& match : matches) {
        MonomerMatch monomer;
        auto monomer_idx = monomers.size();
        for (const auto& [query_idx, atom_idx] : match) {
            auto atom = atomistic_mol.getAtomWithIdx(atom_idx);
            auto map_no = query.attch_map[query_idx];
            if (map_no != TERMINAL_ATTCHPT && atom->hasProp(MONOMER_IDX)) {
                // This shouldn't happen, sanity check
                throw std::runtime_error(fmt::format(
                    "Atom {} belongs to more than 1 monomer", atom->getIdx()));
            }

            // Map number TERMINAL_ATTCHPT corresponds to the R1 in the next
            // monomer, so do not add it if it is present *unless* it is an
            // oxygen which implies the end of the PEPTIDE chain without a
            // special capping group.
            if (map_no != TERMINAL_ATTCHPT ||
                (atom->getAtomicNum() == 8 && atom->getDegree() == 1)) {
                atom->setProp<unsigned int>(MONOMER_IDX, monomer_idx);
                monomer.atom_indices.push_back(atom->getIdx());

                // This helps us identify linkages later
                if (map_no != SIDECHAIN_ATTCHPT) {
                    atom->setProp<unsigned int>(ATTACH_NUM, map_no);
                }
            }

            if (map_no != NO_ATTACHMENT) {
                if (map_no == 1) {
                    // Add methyl group if present (monomers like
                    // N-Methyl-Alanine)
                    for (auto nbr : atomistic_mol.atomNeighbors(atom)) {
                        if (nbr->getAtomicNum() == 6 && nbr->getDegree() == 1) {
                            monomer.atom_indices.push_back(nbr->getIdx());
                            nbr->setProp<unsigned int>(MONOMER_IDX,
                                                       monomer_idx);
                        }
                    }
                    monomer.r1 = atom->getIdx();
                } else if (map_no == 2) {
                    monomer.r2 = atom->getIdx();
                } else if (map_no == 3) {
                    monomer.r3 = atom->getIdx();
                } else if (map_no == SIDECHAIN_ATTCHPT) {
                    // if there is a side chain, the attachment point will be at
                    // the SIDECHAIN_IDX and will be indicated by the presence
                    // of the atom map number.
                    sidechain_attch_pts[monomer_idx] = atom_idx;
                }
            }
        }
        monomers.push_back(monomer);
    }
}

void addSidechainToMonomer(const RDKit::ROMol& atomistic_mol,
                           std::vector<unsigned int>& monomer,
                           unsigned int monomer_idx, unsigned int attch_at_idx)
{
    // BFS but use MONOMER_IDX as visited marker
    std::queue<unsigned int> q;
    q.push(attch_at_idx);
    while (!q.empty()) {
        auto at_idx = q.front();
        q.pop();
        auto at = atomistic_mol.getAtomWithIdx(at_idx);
        if (!at->hasProp(MONOMER_IDX)) {
            at->setProp<unsigned int>(MONOMER_IDX, monomer_idx);
            monomer.push_back(at_idx);
        }
        for (const auto& nbr : atomistic_mol.atomNeighbors(at)) {
            if (!nbr->hasProp(MONOMER_IDX)) {
                q.push(nbr->getIdx());
            }
        }
    }
}

void groupRemainingAtoms(const RDKit::ROMol& atomistic_mol,
                         std::vector<MonomerMatch>& monomers)
{
    /*
     * Group all remaining atoms into their own monomers by connectivity.
     *
     * This is done by performing a BFS search starting from each atom that
     * does not belong to a monomer, determined by the presence of the
     * MONOMER_IDX property.
     */
    std::vector<bool> visited(atomistic_mol.getNumAtoms(), false);
    for (unsigned int i = 0; i < atomistic_mol.getNumAtoms(); ++i) {
        auto at = atomistic_mol.getAtomWithIdx(i);
        if (at->hasProp(MONOMER_IDX)) {
            continue;
        }

        MonomerMatch monomer;
        std::queue<unsigned int> q;
        q.push(i);
        while (!q.empty()) {
            auto at_idx = q.front();
            q.pop();
            auto at = atomistic_mol.getAtomWithIdx(at_idx);
            if (!visited[at_idx]) {
                visited[at_idx] = true;
                monomer.atom_indices.push_back(at_idx);
                for (const auto& nbr : atomistic_mol.atomNeighbors(at)) {
                    if (visited[nbr->getIdx()]) {
                        continue;
                    }
                    if (!nbr->hasProp(MONOMER_IDX)) {
                        q.push(nbr->getIdx());
                    } else {
                        // Deduce attachment point from attachment point on
                        // neighboring monomer, if there.
                        unsigned int neighbor_attach_num = NO_ATTACHMENT;
                        if (nbr->getPropIfPresent(ATTACH_NUM,
                                                  neighbor_attach_num)) {
                            if (neighbor_attach_num == 1) {
                                // Assume this is a backbone linkage; this could
                                // also be an instance where the neighbor is a
                                // sidechain, but that must be deduced after
                                // ordering monomers
                                monomer.r2 = at_idx;
                                at->setProp<unsigned int>(ATTACH_NUM, 2);
                            } else if (neighbor_attach_num == 2) {
                                // This is a backbone linkage
                                monomer.r1 = at_idx;
                                at->setProp<unsigned int>(ATTACH_NUM, 1);
                            } else if (neighbor_attach_num == 3) {
                                // Assume this is an R3-R3 linkage; this could
                                // also be an instance where this monomer is a
                                // sidechain, but that must be deduced after
                                // ordering monomers
                                monomer.r3 = at_idx;
                                at->setProp<unsigned int>(ATTACH_NUM, 3);
                            }
                        }
                    }
                }
            }
        }

        if (monomer.atom_indices.size() == 1) {
            // Unless this is attached to an actual attachment point, we should
            // just add this atom to the monomer it is attached to.
            auto at = atomistic_mol.getAtomWithIdx(monomer.atom_indices[0]);
            if (at->getDegree() == 1) {
                auto neigh = *atomistic_mol.atomNeighbors(at).begin();
                unsigned int neigh_attach_num = 0;
                if (neigh->getPropIfPresent(ATTACH_NUM, neigh_attach_num) &&
                    neigh_attach_num == NO_ATTACHMENT) {
                    auto neigh_monomer_idx =
                        neigh->getProp<unsigned int>(MONOMER_IDX);
                    monomers[neigh_monomer_idx].atom_indices.push_back(
                        monomer.atom_indices[0]);
                    at->setProp<unsigned int>(MONOMER_IDX, neigh_monomer_idx);
                    continue;
                }
            }
        }
        // set MONOMER_IDX for all atoms in the monomer
        for (auto idx : monomer.atom_indices) {
            atomistic_mol.getAtomWithIdx(idx)->setProp<unsigned int>(
                MONOMER_IDX, monomers.size());
        }

        if (!monomer.atom_indices.empty()) {
            monomers.push_back(monomer);
        }
    }
}

void detectLinkages(const RDKit::ROMol& atomistic_mol,
                    std::vector<MonomerMatch>& monomers,
                    std::vector<Linkage>& linkages)
{
    // go through bonds, find bonds that are between different monomers,
    // and add them as linkages
    for (const auto& bond : atomistic_mol.bonds()) {
        auto at1 = bond->getBeginAtom();
        auto at2 = bond->getEndAtom();
        unsigned int monomer_idx1 = 0, monomer_idx2 = 0;
        if (at1->getPropIfPresent(MONOMER_IDX, monomer_idx1) &&
            at2->getPropIfPresent(MONOMER_IDX, monomer_idx2) &&
            monomer_idx1 != monomer_idx2) {
            unsigned int attach_num1 = NO_ATTACHMENT;
            unsigned int attach_num2 = NO_ATTACHMENT;
            at1->getPropIfPresent(ATTACH_NUM, attach_num1);
            at2->getPropIfPresent(ATTACH_NUM, attach_num2);

            if (attach_num1 == NO_ATTACHMENT && attach_num2 == NO_ATTACHMENT) {
                // This means we have two monomers that are connected and
                // the attachment point location cannot be determined from
                // both monomers. Since we cannot determine where one
                // monomer starts and the other ends, throw an error.
                throw std::runtime_error(fmt::format(
                    "Atoms {} and {} are bonded and in different monomers but "
                    "attachment points cannot be determined",
                    at1->getIdx(), at2->getIdx()));
            }

            // If no attachment point is present, we assume it
            // is R3 and store it for future reference
            if (attach_num1 == NO_ATTACHMENT) {
                attach_num1 = 3;
                monomers[monomer_idx1].r3 = at1->getIdx();
                at1->setProp<unsigned int>(ATTACH_NUM, attach_num1);
            }
            if (attach_num2 == NO_ATTACHMENT) {
                attach_num2 = 3;
                monomers[monomer_idx2].r3 = at2->getIdx();
                at2->setProp<unsigned int>(ATTACH_NUM, attach_num2);
            }

            if (attach_num2 > attach_num1) {
                // Ensure directionality of the linkage, R2-R1 or R3-R1
                linkages.emplace_back(monomer_idx2, monomer_idx1, attach_num2,
                                      attach_num1);
            } else {
                if (attach_num1 == 1 && attach_num2 == 1) {
                    // Special case of having a monomer at the beginning of
                    // a chain, attach_num1 should be 2
                    attach_num1 = 2;
                }
                if (attach_num1 == attach_num2 && attach_num1 != 3) {
                    throw std::runtime_error(fmt::format(
                        "Monomers {} and {} are connected with "
                        "the same attachment point {}, something is wrong",
                        monomer_idx1, monomer_idx2, attach_num1));
                }
                linkages.emplace_back(monomer_idx1, monomer_idx2, attach_num1,
                                      attach_num2);
            }
        }
    }
}

/*
 * Break an atomistic molecule into monomers
 *
 * Every atom should belong to either 1 or 2 monomers. If an atom belongs to 2
 * monomers, it represents a connection between the two monomers.
 *
 * Input ROMol is labeled as follows:
 * - MONOMER_IDX: index of the monomer the atom belongs to
 * - ATTACH_NUM: attachment point number (1, 2, 3, or NO_ATTACHMENT) for the
 * attachment point the atom represents within the monomer.
 *
 * Returns a list monomer index sets, which represents the atom indices of each
 * monomer.
 */
void identifyMonomers(RDKit::RWMol& atomistic_mol,
                      std::vector<MonomerMatch>& monomers)
{
    // Approach for identifying monomers:
    // 1. Find all matches with SMARTS queries for amino acids
    // 2. Add side chains to generic matches based on attachment points
    // 3. Identify and group any remaining atoms into 'unclassified' monomers.
    // These are grouped by connectivity.

    // Remove hydrogens that would otherwise not be included in the mapping
    RDKit::MolOps::removeAllHs(atomistic_mol);

    // Maps the monomer index to the atomistic atom index where the sidechain is
    // attached
    std::unordered_map<unsigned int, unsigned int> sidechain_attch_pts;

    // These searches are done in a specific order, from most to least specific
    // query. For example, we want to search with GENERIC_AMINO_ACID_QUERY
    // before GLYCINE_AMINO_ACID_QUERY because the glycine is a substructure of
    // the generic query and would result in incorrect matches. Start with
    // queries that may include R3 attachments (i.e. CYS)
    addMatchesToMonomers(CYSTEINE_QUERY, atomistic_mol, monomers,
                         sidechain_attch_pts);
    addMatchesToMonomers(GENERIC_AMINO_ACID_QUERY, atomistic_mol, monomers,
                         sidechain_attch_pts);
    addMatchesToMonomers(GLYCINE_AMINO_ACID_QUERY, atomistic_mol, monomers,
                         sidechain_attch_pts);

    // Add sidechains to monomers (if not specified in the above queries)
    for (const auto& [monomer_idx, sidechain_attach] : sidechain_attch_pts) {
        addSidechainToMonomer(atomistic_mol, monomers[monomer_idx].atom_indices,
                              monomer_idx, sidechain_attach);
    }
    groupRemainingAtoms(atomistic_mol, monomers);

    // Every atom should belong to a monomer
    for (const auto& atom : atomistic_mol.atoms()) {
        if (!atom->hasProp(MONOMER_IDX)) {
            throw std::runtime_error(fmt::format(
                "Atom {} does not belong to any monomer", atom->getIdx()));
        }
    }
}

void neutralizeAtoms(RDKit::ROMol& mol)
{
    // Algorithm for neutralizing molecules from
    // https://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules by Noel
    // O’Boyle Will neutralize the molecule by adding or removing hydrogens as
    // needed. This will ensure SMILES can be used to match atomistic structures
    // to the correct monomer.
    static const std::unique_ptr<RDKit::RWMol> neutralize_query(
        RDKit::SmartsToMol(
            "[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]"));
    for (const auto& match : RDKit::SubstructMatch(mol, *neutralize_query)) {
        auto atom = mol.getAtomWithIdx(match[0].second);
        auto chg = atom->getFormalCharge();
        auto hcount = atom->getTotalNumHs();
        atom->setFormalCharge(0);
        atom->setNumExplicitHs(hcount - chg);
        atom->updatePropertyCache();
    }
}

std::vector<std::string> enumerate_smiles(const std::string& smiles)
{
    using namespace RDKit::v2::SmilesParse;

    const static SmilesParserParams p{.removeHs = false, .replacements = {}};

    static RDKit::Atom dummy_atom(0);

    std::vector<std::string> enumerated_smiles;

    std::queue<std::unique_ptr<RDKit::RWMol>> q;

    q.emplace(MolFromSmiles(smiles, p));

    while (!q.empty()) {
        auto mol = std::move(q.front());
        q.pop();

        bool found_mapnum = false;
        for (const auto& atom : mol->atoms()) {
            if (atom->hasProp(RDKit::common_properties::molAtomMapNumber)) {
                found_mapnum = true;

                // Clear the map number: We don't want them to be
                // in the final SMILES, and pushing a copy of the
                // modified molecule back onto the queue.
                atom->clearProp(RDKit::common_properties::molAtomMapNumber);
                q.emplace(new RDKit::RWMol(*mol));

                // Now, replace the atom with a dummy, and push
                // the mol back onto the queue.
                mol->replaceAtom(atom->getIdx(), &dummy_atom);
                q.push(std::move(mol));

                // Only handle one map number at a time
                break;
            }
        }

        // This mol has been fully enumerated, so push the SMILES
        // to the output vector.
        if (!found_mapnum) {
            // The mols we check against have all Hs removed.
            RDKit::MolOps::removeAllHs(*mol);
            enumerated_smiles.push_back(RDKit::MolToSmiles(*mol));
        }
    }

    return enumerated_smiles;
}

std::unordered_map<std::string, std::string> enumerate_all_core_smiles()
{
    std::unordered_map<std::string, std::string> core_smiles_to_monomer;

    const auto& db = MonomerDatabase::instance();

    for (auto& [smiles, symbol] : db.getAllSMILES()) {
        for (auto&& enumerated_smiles : enumerate_smiles(smiles)) {
            // TO DO: should we check for duplicates ? how do we handle them?
            core_smiles_to_monomer.emplace(std::move(enumerated_smiles),
                                           symbol);
        }
    }

    return core_smiles_to_monomer;
}

std::optional<std::string>
findHelmSymbol(const RDKit::ROMol& atomistic_mol,
               const std::vector<unsigned int>& atom_indices)
{
    static RDKit::Atom dummy_atom(0);

    // This may get big and slow. Currently, enumeration of the 62 default
    // monomer definitions in default_monomer_definitions.json into 292
    // possible substitution states takes about 20 ms on my linux desktop.
    // But it may get untractable if we need to handle thousands of monomer
    // definitions or more.
    static auto amino_acids = enumerate_all_core_smiles();

    auto mol_fragment = ExtractMolFragment(atomistic_mol, atom_indices, false);

    constexpr bool update_label = false;
    constexpr bool take_ownership = false; // make a copy
    mol_fragment->beginBatchEdit();
    for (auto* at : mol_fragment->atoms()) {
        if (auto attach_num = NO_ATTACHMENT;
            at->getPropIfPresent(ATTACH_NUM, attach_num) &&
            attach_num != NO_ATTACHMENT) {

            if (attach_num == TERMINAL_ATTCHPT) {
                // This is a chain terminating O/N. We need to remove it, since
                // in the enumeration we replace these with dummy atoms
                mol_fragment->removeAtom(at->getIdx());
            } else {
                auto atom_idx = mol_fragment->addAtom(&dummy_atom, update_label,
                                                      take_ownership);
                mol_fragment->addBond(atom_idx, at->getIdx(),
                                      RDKit::Bond::SINGLE);
            }
        }
    }
    mol_fragment->commitBatchEdit();
    // we only replaced attachment points with dummy atoms, I think
    // we shouldn't require a ring finding or a sanitization,
    // just updating the valence (RDKit now resets them when adding bonds).
    mol_fragment->updatePropertyCache(false);

    auto monomer_smiles = RDKit::MolToSmiles(*mol_fragment, true);

    // Using the enumerated amino acids
    std::string monomer_symbol;
    if (amino_acids.find(monomer_smiles) != amino_acids.end()) {
        monomer_symbol = amino_acids.at(monomer_smiles);
        return monomer_symbol;
    }

    return std::nullopt;
}

void buildMonomerMol(const RDKit::ROMol& atomistic_mol,
                     std::vector<MonomerMatch>& monomers,
                     boost::shared_ptr<RDKit::RWMol> monomer_mol,
                     std::vector<Linkage>& linkages)
{
    // Start with all atoms in a single peptide chain
    monomer_mol->setProp<bool>(HELM_MODEL, true);

    // In the unique case only one monomer is identified and it has no known
    // monomer definition (meaning it is a SMILES monomer), place the entire
    // structure into a single CHEM smiles monomer.
    if (monomers.size() == 1) {
        auto& monomer = monomers.front();
        auto helm_symbol = findHelmSymbol(atomistic_mol, monomer.atom_indices);
        if (!helm_symbol) {
            addMonomer(*monomer_mol, RDKit::MolToSmiles(atomistic_mol), 1,
                       "CHEM1", MonomerType::SMILES);
            return;
        }
    }

    static constexpr const char* CHAIN_ID = "PEPTIDE1";

    int residue_num = 1;
    for (const auto& monomer : monomers) {
        auto helm_symbol = findHelmSymbol(atomistic_mol, monomer.atom_indices);

        // If the monomer is a known amino acid, use the 1-letter code
        if (helm_symbol) {
            addMonomer(*monomer_mol, *helm_symbol, residue_num, CHAIN_ID);
        } else {
            // We need to add R1/R2 attachment points to the monomer
            RDKit::RWMol atomistic_copy(atomistic_mol);

            static constexpr const char* ATTCH_PROP = "attach";
            if (monomer.r1 != NO_ATTACHMENT) {
                atomistic_mol.getAtomWithIdx(monomer.r1)
                    ->setProp<unsigned int>(ATTCH_PROP, 1);
            }
            if (monomer.r2 != NO_ATTACHMENT) {
                atomistic_mol.getAtomWithIdx(monomer.r2)
                    ->setProp<unsigned int>(ATTCH_PROP, 2);
            }
            if (monomer.r3 != NO_ATTACHMENT) {
                // Check if this monomer actually has an R3 attachment point
                bool has_r3 = false;
                auto this_at = atomistic_mol.getAtomWithIdx(monomer.r3);
                auto monomer_idx = this_at->getProp<unsigned int>(MONOMER_IDX);
                for (auto& neigh : atomistic_mol.atomNeighbors(this_at)) {
                    unsigned int neigh_monomer_idx = NO_ATTACHMENT;
                    if (neigh->getPropIfPresent(MONOMER_IDX,
                                                neigh_monomer_idx) &&
                        neigh_monomer_idx != monomer_idx) {
                        has_r3 = true;
                        break;
                    }
                }
                if (has_r3) {
                    atomistic_mol.getAtomWithIdx(monomer.r3)
                        ->setProp<unsigned int>(ATTCH_PROP, 3);
                }
            }

            auto mol_fragment =
                ExtractMolFragment(atomistic_mol, monomer.atom_indices, false);

            static constexpr bool update_label = true;
            static constexpr bool take_ownership = true;
            for (auto at : mol_fragment->atoms()) {
                unsigned int map_no;
                if (at->getPropIfPresent(ATTCH_PROP, map_no)) {
                    // add dummy atom with this attachment point
                    auto new_at = new RDKit::Atom(0);
                    new_at->setProp(RDKit::common_properties::molAtomMapNumber,
                                    map_no);
                    auto new_at_idx = mol_fragment->addAtom(
                        new_at, update_label, take_ownership);
                    mol_fragment->addBond(new_at_idx, at->getIdx(),
                                          RDKit::Bond::SINGLE);
                }
            }
            auto monomer_smiles = RDKit::MolToSmiles(*mol_fragment);

            addMonomer(*monomer_mol, monomer_smiles, residue_num, CHAIN_ID,
                       MonomerType::SMILES);
        }
        ++residue_num;
    }

    for (const auto& link : linkages) {
        addConnection(*monomer_mol, link.monomer_idx1, link.monomer_idx2,
                      link.to_string());
    }
}

size_t findStartMonomer(const std::vector<int>& parents)
{
    // Find the start monomer by picking arbitrary monomer then following
    // edges until a cycle is detected or the source is reached.
    auto current_monomer = 0;
    std::vector<bool> visited(parents.size(), false);
    while (parents[current_monomer] != -1 &&
           !visited[parents[current_monomer]]) {
        visited[current_monomer] = true;
        current_monomer = parents[current_monomer];
    }
    return static_cast<size_t>(current_monomer);
}

void orderMonomers(RDKit::ROMol& atomistic_mol,
                   std::vector<MonomerMatch>& monomers,
                   std::vector<Linkage>& linkages)
{
    // Create adjacency list representation of the monomer connectivity graph
    // that is partially directed; R2->R1 and R3->R1 connections are directed,
    // while R3-R3 connections are undirected.
    std::vector<std::vector<std::pair<size_t, Linkage*>>> adj_list(
        monomers.size());
    std::vector<int> parents(monomers.size(), -1);
    for (auto& link : linkages) {
        if (link.attach_from != link.attach_to) {
            // This skips disulfide (or other R3-R3 bonds) because that would
            // indicate the need for multiple chains, see SHARED-10787 Ensure we
            // go R2->R1 and R3->R1 for directionality
            adj_list[link.monomer_idx1].push_back({link.monomer_idx2, &link});

            // Since the directed bonds always go into R1, each monomer can have
            // at most one parent. Disulfide linkages are ignored for finding
            // the start monomer.
            parents[link.monomer_idx2] = link.monomer_idx1;
        }
    }
    auto start_monomer = findStartMonomer(parents);
    auto start_monomer_atom_idx = monomers[start_monomer].atom_indices.front();

    // Build ordered list starting from the best candidate
    std::vector<size_t> new_monomer_order;
    std::vector<bool> visited(monomers.size(), false);

    // Helper function for DFS traversal that follows R2->R1 directionality
    std::function<void(size_t)> dfs = [&](size_t current) {
        visited[current] = true;
        new_monomer_order.push_back(current);

        // Sort neighbors to prioritize R2->R1 connections
        std::vector<std::pair<int, size_t>> sorted_neighbors;
        for (auto& [neighbor, link] : adj_list[current]) {
            if (!visited[neighbor]) {
                // Determine priority:
                // * R3->R1 gets highest priority as it denotes a branched
                // connection
                //       i.e., A.B(C)D.E should be ordered as A-B-C-D-E not
                //       A-B-D-E-C
                // * R2->R1 gets next priority as it denotes a backbone
                // connection
                // * R3-R3 connections are ignored, see SHARED-10787
                int priority = 0;
                if (link->attach_from == 3 && link->attach_to == 1) {
                    priority = 2;
                } else if (link->attach_from == 2 && link->attach_to == 1) {
                    priority = 1;
                } else {
                    // Go N -> C direction only, no R1->R3 or R2->R3
                    continue;
                }
                sorted_neighbors.push_back({priority, neighbor});
            }
        }

        // Visit neighbors in priority order
        std::sort(
            sorted_neighbors.begin(), sorted_neighbors.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });
        for (auto& [_priority, neighbor] : sorted_neighbors) {
            if (!visited[neighbor]) {
                dfs(neighbor);
            }
        }
    };

    dfs(start_monomer);

    // If not all monomers were visited (disjoint chains), we have a
    // disconnected graph -- for now, throw an error.
    for (size_t i = 0; i < monomers.size(); ++i) {
        if (!visited[i]) {
            auto this_at_idx = monomers[i].atom_indices.front();
            throw std::runtime_error(fmt::format(
                "Could not traverse bonds from atom {} to atom {} - chains "
                "with disconnections or missing residues are not supported",
                start_monomer_atom_idx, this_at_idx));
        }
    }

    // Determine index mapping from old to new order
    std::vector<size_t> index_mapping(monomers.size());
    for (size_t i = 0; i < new_monomer_order.size(); ++i) {
        index_mapping[new_monomer_order[i]] = static_cast<size_t>(i);
    }

    // sort the monomers according to the new order
    std::vector<MonomerMatch> ordered_monomers(monomers.size());
    for (size_t i = 0; i < new_monomer_order.size(); ++i) {
        ordered_monomers[i] = monomers[new_monomer_order[i]];
    }
    monomers = std::move(ordered_monomers);

    // Update linkages with new monomer indices
    for (auto& link : linkages) {
        link.monomer_idx1 = index_mapping[link.monomer_idx1];
        link.monomer_idx2 = index_mapping[link.monomer_idx2];
    }
}

void removeWaters(RDKit::RWMol& mol)
{
    mol.beginBatchEdit();
    auto is_water = [](const RDKit::Atom* atom) {
        const auto res_info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());

        // TODO: This seems like it shouldn't be the job of this function; there
        // should be some sort of separate preprocessing step that removes
        // waters and other unwanted residues
        if (!res_info) {
            return false;
        }
        // strip whitespace from residue name
        auto res_name = res_info->getResidueName();
        res_name.erase(std::remove(res_name.begin(), res_name.end(), ' '),
                       res_name.end());
        if (res_info && res_name == "HOH") {
            return true;
        }
        return false;
    };
    for (auto atom : mol.atoms()) {
        if (is_water(atom)) {
            mol.removeAtom(atom);
        }
    }
    mol.commitBatchEdit();
}

void findChainsAndResidues(const RDKit::ROMol& mol,
                           ChainsAndResidues& chains_and_residues)
{
    // Find all chains and residues in the molecule
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        const auto atom = mol.getAtomWithIdx(i);
        const auto res_info = static_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());

        if (!res_info) {
            if (atom->getAtomicNum() == 1) {
                // This is a hydrogen that is part of a chain, so add it to the
                // chain
                continue;
            }
            throw std::runtime_error(fmt::format(
                "Atom {} does not have residue info", atom->getIdx()));
        }
        const auto chain_id = res_info->getChainId();
        const auto res_num = res_info->getResidueNumber();
        const auto ins_code = res_info->getInsertionCode();
        chains_and_residues[chain_id][std::make_pair(res_num, ins_code)]
            .push_back(i);
    }
}

std::string getMonomerSmiles(RDKit::ROMol& mol,
                             const std::vector<unsigned int>& atom_idxs,
                             const std::string& chain_id,
                             const std::pair<int, std::string>& current_key,
                             int res_num, bool end_of_chain)
{
    // Determine the atoms in current_res that connect to adjacent residues
    std::vector<std::pair<int, int>> attch_idxs; // adjacent res num, atom_idx
    for (auto idx : atom_idxs) {
        const auto at = mol.getAtomWithIdx(idx);
        for (const auto neigh : mol.atomNeighbors(at)) {
            if (neigh->getAtomicNum() == 1 && !neigh->hasProp(CHAIN_HYDROGEN)) {
                // skip Hs
                continue;
            }
            const auto res_info =
                dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                    neigh->getMonomerInfo());
            const auto key = std::make_pair(res_info->getResidueNumber(),
                                            res_info->getInsertionCode());
            if (key != current_key || res_info->getChainId() != chain_id) {
                // neighbor is in different residue, this will be an attachment
                // point
                attch_idxs.push_back(
                    {res_info->getResidueNumber(), at->getIdx()});
            }
        }
    }

    // This part is difficult; we need to figure out the order of the attachment
    // points on this monomer, Ideally we'd find the backbone and ensure to
    // label them directionally correctly using R2-R1
    std::ranges::sort(
        attch_idxs.begin(), attch_idxs.end(),
        [&res_num](const std::pair<int, int>& a, const std::pair<int, int>& b) {
            // we want the most 'adjacent' residue first in terms of residue
            // ordering; this gets more complicated if we have different chains
            if (std::abs(a.first - res_num) == std::abs(b.first - res_num)) {
                return a.first < b.first;
            }
            return std::abs(a.first - res_num) < std::abs(b.first - res_num);
        });

    static constexpr bool sanitize = false;
    auto mol_fragment = ExtractMolFragment(mol, atom_idxs, sanitize);

    // Add dummy atoms with an attachment point #s
    // For now, all monomers have at least attachment points 1 and 2 so
    // that backbone bonds can be formed (except the beginning monomer)
    int current_attchpt = 1;

    // If this is the beginning of the chain or a branch, start with attachment
    // point 2
    if (!end_of_chain && res_num == 1) {
        current_attchpt = 2;
    }

    static constexpr bool update_label = true;
    static constexpr bool take_ownership = true;
    for (const auto& [_, ref_idx] : attch_idxs) {
        for (auto at : mol_fragment->atoms()) {
            int orig_idx;
            if (at->getPropIfPresent(REFERENCE_IDX, orig_idx) &&
                orig_idx == ref_idx) {
                auto new_at = new RDKit::Atom(0);
                new_at->setProp(RDKit::common_properties::molAtomMapNumber,
                                current_attchpt);
                auto new_at_idx =
                    mol_fragment->addAtom(new_at, update_label, take_ownership);
                mol_fragment->addBond(new_at_idx, at->getIdx(),
                                      RDKit::Bond::SINGLE);
                mol.getAtomWithIdx(orig_idx)->setProp<unsigned int>(
                    MONOMER_MAP_NUM, current_attchpt);
                ++current_attchpt;
            }
        }
    }
    // There should always be enough attachment points so that
    // backbone connections can be made (R1 and R2)
    // TODO: Should this indicate a new chain?
    while (current_attchpt <= MIN_ATTCHPTS) {
        if (end_of_chain && current_attchpt > 1) {
            break;
        }
        auto new_at = new RDKit::Atom(0);
        new_at->setProp(RDKit::common_properties::molAtomMapNumber,
                        current_attchpt);
        mol_fragment->addAtom(new_at, update_label, take_ownership);
        ++current_attchpt;
    }

    // removing hydrogens to keep HELM string readable
    removeHs(*mol_fragment);
    return RDKit::MolToSmiles(*mol_fragment);
}

bool sameMonomer(RDKit::RWMol& atomistic_mol,
                 const std::vector<unsigned int>& atom_idxs,
                 const std::string& db_smiles)
{
    // Occasionally SMILES that cannot be kekulized are extracted from atomistic
    // mol, so skip sanitization
    constexpr int debug = 0;
    constexpr bool sanitize = false;

    // Monomer from original atomistic mol and monomer defined by database
    auto monomer_frag = ExtractMolFragment(atomistic_mol, atom_idxs, sanitize);
    std::unique_ptr<RDKit::RWMol> db_mol(
        RDKit::SmilesToMol(db_smiles, debug, sanitize));

    // Remove stereochemistry, atom map numbers, and neutralize the atoms
    // leaving groups are removed since they won't always be included in the
    // residue extracted from the atomistic molecule
    auto clean_mol = [](RDKit::RWMol& mol) {
        RDKit::MolOps::removeStereochemistry(mol);
        mol.beginBatchEdit();
        for (auto at : mol.atoms()) {
            unsigned int map_no;
            if (at->getPropIfPresent(RDKit::common_properties::molAtomMapNumber,
                                     map_no)) {
                // Label parent with this map number so we know which map number
                // the linkage uses
                for (const auto& nbr : mol.atomNeighbors(at)) {
                    nbr->setProp(RDKit::common_properties::molAtomMapNumber,
                                 map_no);
                }

                mol.removeAtom(at);
            }
        }
        mol.commitBatchEdit();
        removeHs(mol);
        // set aromaticity
        RDKit::MolOps::setAromaticity(mol);
    };
    clean_mol(*monomer_frag);
    clean_mol(*db_mol);

    // The DB monomer has had the leaving groups removed, while the residue
    // extracted from the atomistic mol may still have them present if it is at
    // the beginning or end of a chain. As a result, we need to allow for the DB
    // monomer to have one less atom than the residue extracted from the
    // atomistic mol
    auto match = RDKit::SubstructMatch(*monomer_frag,
                                       *db_mol); // (queryAtomIdx, molAtomIdx)

    if (match.size()) {

        // Any unmapped atom in monomer_frag must map to a leaving group in the
        // db_mol create vector of atom indices in the monomer_frag that are not
        // in the match
        for (auto at : monomer_frag->atoms()) {
            int target_idx = at->getIdx();
            auto it = std::find_if(match[0].begin(), match[0].end(),
                                   [target_idx](const std::pair<int, int>& p) {
                                       return p.second == target_idx;
                                   });
            if (it == match[0].end()) {
                // This atom should be terminal, otherwise the match is
                // definitely wrong
                if (at->getDegree() > 1) {
                    return false;
                }

                // Unmatched atoms should be leaving groups, we can check this
                // by matching
                RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
                boost::tie(nbrIdx, endNbrs) =
                    monomer_frag->getAtomNeighbors(at);

                it = std::find_if(match[0].begin(), match[0].end(),
                                  [nbrIdx](const std::pair<int, int>& p) {
                                      return p.second ==
                                             static_cast<int>(*nbrIdx);
                                  });
                if (it == match[0].end()) {
                    return false;
                }

                auto db_atom = db_mol->getAtomWithIdx(it->first);
                if (!db_atom->hasProp(
                        RDKit::common_properties::molAtomMapNumber)) {
                    return false;
                }
            }
        }

        // we are mapping atomistic atom indices to attachment point numbers
        // (molAtomMapNumber)
        for (auto [db_idx, at_idx] : match[0]) {
            auto db_at = db_mol->getAtomWithIdx(db_idx);
            auto at = monomer_frag->getAtomWithIdx(at_idx);
            unsigned int map_no;
            if (db_at->getPropIfPresent(
                    RDKit::common_properties::molAtomMapNumber, map_no)) {
                auto ref_idx = at->getProp<unsigned int>(REFERENCE_IDX);
                atomistic_mol.getAtomWithIdx(ref_idx)->setProp<unsigned int>(
                    MONOMER_MAP_NUM, map_no);
            }
        }
        return true;
    }

    return false;
}

int getAttchpt(const RDKit::Atom& monomer, const RDKit::Bond& bond,
               const RDKit::RWMol& atomistic_mol)
{
    // This should be the index of the atom in this monomer
    auto target_atom_idx = bond.getBeginAtomIdx();
    if (atomistic_mol.getAtomWithIdx(target_atom_idx)
            ->getProp<unsigned int>(MONOMER_IDX) != monomer.getIdx()) {
        target_atom_idx = bond.getEndAtomIdx();
    }
    unsigned int attchpt;
    if (atomistic_mol.getAtomWithIdx(target_atom_idx)
            ->getPropIfPresent(MONOMER_MAP_NUM, attchpt)) {
        return attchpt;
    } else {
        return -1;
    }
}

void detectLinkages(RDKit::RWMol& monomer_mol,
                    const RDKit::RWMol& atomistic_mol)
{
    // Find all linkages between monomers (used when PDB residue info is
    // available)
    for (const auto& bond : atomistic_mol.bonds()) {
        const auto begin_atom = bond->getBeginAtom();
        const auto end_atom = bond->getEndAtom();

        unsigned int begin_monomer_idx =
            std::numeric_limits<unsigned int>::max();
        unsigned int end_monomer_idx = std::numeric_limits<unsigned int>::max();
        if (!begin_atom->getPropIfPresent(MONOMER_IDX, begin_monomer_idx) ||
            !end_atom->getPropIfPresent(MONOMER_IDX, end_monomer_idx)) {
            // One of these is a hydrogen or unmapped for some reason
            continue;
        }

        // Check if this is already present as a backbone linkage
        if (begin_monomer_idx == end_monomer_idx ||
            monomer_mol.getBondBetweenAtoms(begin_monomer_idx,
                                            end_monomer_idx) != nullptr) {
            continue;
        }

        const auto begin_monomer =
            monomer_mol.getAtomWithIdx(begin_monomer_idx);
        const auto end_monomer = monomer_mol.getAtomWithIdx(end_monomer_idx);

        auto begin_attchpt = getAttchpt(*begin_monomer, *bond, atomistic_mol);
        auto end_attchpt = getAttchpt(*end_monomer, *bond, atomistic_mol);

        if (begin_attchpt == -1 || end_attchpt == -1) {
            // This happens when the input atomistic mol has bonds between
            // residues on atoms that are not marked as attachment points in the
            // monomer DB. This is important to capture so that we can determine
            // where we may need to add attachment points.
            auto begin_res_info =
                dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                    begin_atom->getMonomerInfo());
            auto end_res_info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                end_atom->getMonomerInfo());
            if (begin_res_info && end_res_info) {
                std::cerr << fmt::format(
                    "Bond between residue {}{} in chain {} and residue {}{} "
                    "in chain {} do not correspond to attachment points in "
                    "their residues.\n",
                    begin_res_info->getResidueName(),
                    begin_res_info->getResidueNumber(),
                    begin_res_info->getChainId(),
                    end_res_info->getResidueName(),
                    end_res_info->getResidueNumber(),
                    end_res_info->getChainId());
            } else {
                std::cerr << fmt::format(
                    "Bond between atoms {} and {} in monomers {} and {} do not "
                    "correspond to attachment points in their residues.\n",
                    begin_atom->getIdx(), end_atom->getIdx(),
                    begin_monomer->getProp<std::string>(ATOM_LABEL),
                    end_monomer->getProp<std::string>(ATOM_LABEL));
            }
            continue;
        }

        // Backbone connections (R2-R1) should be added in that order when
        // possible.
        if (begin_attchpt == 2 && end_attchpt == 1) {
            addConnection(monomer_mol, begin_monomer_idx, end_monomer_idx,
                          "R2-R1");
        } else if (begin_attchpt == 1 && end_attchpt == 2) {
            addConnection(monomer_mol, end_monomer_idx, begin_monomer_idx,
                          "R2-R1");
        } else if (begin_monomer_idx < end_monomer_idx) {
            addConnection(monomer_mol, begin_monomer_idx, end_monomer_idx,
                          fmt::format("R{}-R{}", begin_attchpt, end_attchpt));
        } else {
            addConnection(monomer_mol, end_monomer_idx, begin_monomer_idx,
                          fmt::format("R{}-R{}", end_attchpt, begin_attchpt));
        }
    }
}

std::optional<std::tuple<std::string, std::string, ChainType>>
getHelmInfo(const MonomerDatabase& db, const RDKit::Atom* atom,
            bool has_pdb_codes)
{
    if (has_pdb_codes) {
        // This comes from something like a PDB or MAE file
        auto res_info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());
        auto res_name = res_info->getResidueName();
        res_name.erase(
            std::remove_if(res_name.begin(), res_name.end(), ::isspace),
            res_name.end());
        return db.getHelmInfo(res_name);
    } else {
        // This comes from a SD file with SUP groups
        auto res_info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());
        std::string sup_class;
        if (!res_info ||
            !atom->getPropIfPresent(SGROUP_PROP_CLASS, sup_class)) {
            return std::nullopt;
        }

        auto pos = BIOVIA_CHAIN_TYPE_MAP.find(sup_class);
        if (pos == BIOVIA_CHAIN_TYPE_MAP.end()) {
            return std::nullopt;
        }
        auto chain_type = pos->second;
        auto res_name =
            res_info->getResidueName(); // this is the label, not the PDB code
        auto smiles = db.getMonomerSmiles(res_name, chain_type);

        if (smiles) {
            return std::make_tuple(res_name, *smiles, chain_type);
        }
        return std::nullopt;
    }
}

boost::shared_ptr<RDKit::RWMol>
pdbInfoAtomisticToMM(const RDKit::ROMol& input_mol, bool has_pdb_codes)
{
    // Make RWMol and remove waters
    RDKit::RWMol mol(input_mol);
    removeWaters(mol);

    // Set reference index for SMILES fragments
    for (auto at : mol.atoms()) {
        at->setProp(REFERENCE_IDX, at->getIdx());
    }

    // Map chain_id -> {residue mols}
    ChainsAndResidues chains_and_residues;
    findChainsAndResidues(mol, chains_and_residues);

    auto& db = MonomerDatabase::instance();

    std::map<ChainType, unsigned int> chain_counts = {{ChainType::PEPTIDE, 0},
                                                      {ChainType::RNA, 0},
                                                      {ChainType::DNA, 0},
                                                      {ChainType::CHEM, 0}};
    auto monomer_mol = boost::make_shared<RDKit::RWMol>();
    for (const auto& [chain_id, residues] : chains_and_residues) {
        // Use first residue to determine chain type. We assume that PDB data
        // is correct and there aren't multiple chain types in a single chain.
        // Default chain type is PEPTIDE if not specified.
        auto helm_info = getHelmInfo(
            db, mol.getAtomWithIdx(residues.begin()->second[0]), has_pdb_codes);
        auto chain_type =
            helm_info ? std::get<2>(*helm_info) : ChainType::PEPTIDE;
        std::string helm_chain_id = fmt::format("{}{}", toString(chain_type),
                                                ++chain_counts[chain_type]);
        // Assuming residues are ordered correctly
        size_t res_num = 1;
        for (const auto& [key, atom_idxs] : residues) {
            helm_info = getHelmInfo(db, mol.getAtomWithIdx(atom_idxs[0]),
                                    has_pdb_codes);
            bool end_of_chain = res_num == residues.size();
            size_t this_monomer;
            if (helm_info &&
                sameMonomer(mol, atom_idxs, std::get<1>(*helm_info))) {
                // Standard residue in monomer DB, Verify that the fragment
                // labeled as the residue matches what is in the monomer
                // database
                this_monomer = addMonomer(*monomer_mol, std::get<0>(*helm_info),
                                          res_num, helm_chain_id);
            } else {
                auto smiles = getMonomerSmiles(mol, atom_idxs, chain_id, key,
                                               res_num, end_of_chain);
                this_monomer = addMonomer(*monomer_mol, smiles, res_num,
                                          helm_chain_id, MonomerType::SMILES);
            }

            // Track which atoms are in which monomer
            for (auto idx : atom_idxs) {
                mol.getAtomWithIdx(idx)->setProp<unsigned int>(MONOMER_IDX,
                                                               this_monomer);
            }
            ++res_num;
        }
    }
    detectLinkages(*monomer_mol, mol);
    return monomer_mol;
}

bool hasPdbResidueInfo(const RDKit::ROMol& mol)
{
    return std::any_of(mol.atoms().begin(), mol.atoms().end(),
                       [](const auto& atom) {
                           return atom->getMonomerInfo() &&
                                  atom->getMonomerInfo()->getMonomerType() ==
                                      RDKit::AtomMonomerInfo::PDBRESIDUE;
                       });
}
} // unnamed namespace

boost::shared_ptr<RDKit::RWMol> toMonomeric(const RDKit::ROMol& mol,
                                            bool try_residue_info)
{
    if (mol.getNumAtoms() == 0) {
        throw std::runtime_error(
            "Input molecule has no atoms, cannot convert to "
            "monomeric representation.");
    }
    // First attempt to use residue information to build the monomeric molecule,
    // if the information isn't present fall back to SMARTS-based method
    RDKit::RWMol atomistic_mol(mol);
    neutralizeAtoms(atomistic_mol);
    if (try_residue_info) {
        if (hasPdbResidueInfo(atomistic_mol)) {
            auto monomer_mol = pdbInfoAtomisticToMM(atomistic_mol, true);
            assignChains(*monomer_mol);
            return monomer_mol;
        } else if (processSupGroups(atomistic_mol)) {
            auto monomer_mol = pdbInfoAtomisticToMM(atomistic_mol, false);
            assignChains(*monomer_mol);
            return monomer_mol;
        }
    }

    std::vector<MonomerMatch> monomers;
    identifyMonomers(atomistic_mol, monomers);

    std::vector<Linkage> linkages;
    detectLinkages(atomistic_mol, monomers, linkages);
    orderMonomers(atomistic_mol, monomers, linkages);

    boost::shared_ptr<RDKit::RWMol> monomer_mol =
        boost::make_shared<RDKit::RWMol>();
    buildMonomerMol(atomistic_mol, monomers, monomer_mol, linkages);
    assignChains(*monomer_mol);
    return monomer_mol;
}

std::vector<std::vector<unsigned int>> getMonomers(const RDKit::ROMol& mol)
{
    RDKit::RWMol atomistic_mol(mol);
    std::vector<MonomerMatch> monomers;
    identifyMonomers(atomistic_mol, monomers);

    std::vector<std::vector<unsigned int>> monomers_indices;
    for (auto& monomer : monomers) {
        monomers_indices.push_back(std::move(monomer.atom_indices));
    }
    return monomers_indices;
}

} // namespace rdkit_extensions
} // namespace schrodinger
