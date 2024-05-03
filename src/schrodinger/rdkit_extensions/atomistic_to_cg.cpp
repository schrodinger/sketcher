#include "schrodinger/rdkit_extensions/atomistic_to_cg.h"

#include <queue>
#include <unordered_map>
#include <memory>

#include <fmt/format.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

#include "schrodinger/rdkit_extensions/coarse_grain.h"
#include "schrodinger/rdkit_extensions/helm.h"

namespace schrodinger
{
namespace rdkit_extensions
{

namespace
{

const std::string ATTACH_FROM{"attachFrom"};
const std::string MONOMER_IDX1{"monomerIndex1"};
const std::string MONOMER_IDX2{"monomerIndex2"};

constexpr int SIDECHAIN_IDX = 2;
constexpr auto NO_ATTACHMENT = std::numeric_limits<unsigned int>::max();

// attachment points 1 and 2 are backbone attachment points
// and 3 is the side chain attachment point
const std::string GENERIC_AMINO_ACID_QUERY =
    "[NX3,NX4+:1][CX4H]([*:3])[CX3](=[OX1])[O,N:2]";
const std::string GLYCINE_AMINO_ACID_QUERY =
    "[N:1][CX4H2][CX3](=[OX1])[O,N:2]"; // no side chain

// SMILES monomer to CG monomer abbreviation mapping
// temporary for now, for proof of concept
// got most of these from pubchem. include the version with N and O
const std::unordered_map<std::string, std::string> amino_acids = {
    {"CC(N)C(=O)O", "A"},                 // Alanine (Ala)
    {"NC(N)=NCCCC(N)C(=O)O", "R"},        // Arginine (Arg)
    {"NC(=O)CC(N)C(=O)O", "N"},           // Asparagine (Asn)
    {"NC(CC(=O)O)C(=O)O", "D"},           // Aspartic acid (Asp)
    {"NC(CS)C(=O)O", "C"},                // Cysteine (Cys)
    {"NC(=O)CCC(N)C(=O)O", "Q"},          // Glutamine (Gln)
    {"NC(CCC(=O)O)C(=O)O", "E"},          // Glutamic acid (Glu)
    {"NCC(=O)O", "G"},                    // Glycine (Gly)
    {"NC(Cc1cnc[nH]1)C(=O)O", "H"},       // Histidine (His)
    {"CCC(C)C(N)C(=O)O", "I"},            // Isoleucine (Ile)
    {"CC(C)CC(N)C(=O)O", "L"},            // Leucine (Leu)
    {"NCCCCC(N)C(=O)O", "K"},             // Lysine (Lys)
    {"CSCCC(N)C(=O)O", "M"},              // Methionine (Met)
    {"NC(Cc1ccccc1)C(=O)O", "F"},         // Phenylalanine (Phe)
    {"O=C(O)C1CCCN1", "P"},               // Proline (Pro)
    {"NC(CO)C(=O)O", "S"},                // Serine (Ser)
    {"CC(O)C(N)C(=O)O", "T"},             // Threonine (Thr)
    {"NC(Cc1c[nH]c2ccccc12)C(=O)O", "W"}, // Tryptophan (Trp)
    {"NC(Cc1ccc(O)cc1)C(=O)O", "Y"},      // Tyrosine (Tyr)
    {"CC(C)C(N)C(=O)O", "V"},             // Valine (Val)
    {"CC(N)C(N)=O", "A"},
    {"NC(=O)C(N)CCCN=C(N)N", "R"}, // arginine, pubchem version
    {"N=C(N)NCCCC(N)C(N)=O", "R"}, // arginine, from HELM paper examples
                                   // (different double bond placement)
    {"NC(=O)CC(N)C(N)=O", "N"},
    {"NC(=O)C(N)CC(=O)O", "D"},
    {"NC(=O)C(N)CS", "C"},
    {"NC(=O)CCC(N)C(N)=O", "Q"},
    {"NC(=O)C(N)CCC(=O)O", "E"},
    {"NCC(N)=O", "G"},
    {"NC(=O)C(N)Cc1cnc[nH]1", "H"},
    {"CCC(C)C(N)C(N)=O", "I"},
    {"CC(C)CC(N)C(N)=O", "L"},
    {"NCCCCC(N)C(N)=O", "K"},
    {"CSCCC(N)C(N)=O", "M"},
    {"NC(=O)C(N)Cc1ccccc1", "F"},
    {"NC(=O)C1CCCN1", "P"},
    {"NC(=O)C(N)CO", "S"},
    {"CC(O)C(N)C(N)=O", "T"},
    {"NC(=O)C(N)Cc1c[nH]c2ccccc12", "W"},
    {"NC(=O)C(N)Cc1ccc(O)cc1", "Y"},
    {"CC(C)C(N)C(N)=O", "V"}};

struct Linkage {
    unsigned int monomer_idx1;
    unsigned int monomer_idx2;
    unsigned int attach_from;
    unsigned int attach_to;
    std::string to_string() const
    {
        return fmt::format("R{}-R{}", attach_from, attach_to);
    }
};

bool already_matched(const RDKit::ROMol& mol,
                     const std::vector<unsigned int>& ids)
{
    // Make sure this match hasn't already been accounted for by a previous
    // match
    for (auto id : ids) {
        auto at = mol.getAtomWithIdx(id);
        unsigned int attch;
        if (at->hasProp(MONOMER_IDX1) &&
            at->getPropIfPresent(ATTACH_FROM, attch) &&
            attch == NO_ATTACHMENT) {
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
void add_matches_to_monomers(
    const std::string& smarts_query, RDKit::ROMol& atomistic_mol,
    std::vector<std::vector<int>>& monomers,
    std::unordered_map<unsigned int, unsigned int>& attch_pts,
    std::vector<Linkage>& linkages)
{
    std::unique_ptr<RDKit::RWMol> query(RDKit::SmartsToMol(smarts_query));
    // maps SMARTS query index to attachment point #
    std::vector<unsigned int> attch_map(query->getNumAtoms(), NO_ATTACHMENT);
    for (const auto atom : query->atoms()) {
        if (atom->hasProp(RDKit::common_properties::molAtomMapNumber)) {
            attch_map[atom->getIdx()] = atom->getProp<unsigned int>(
                RDKit::common_properties::molAtomMapNumber);
        }
    }

    // Set a final function check that ensures the entire match has not
    // already been accounted for by a previous SMARTS search.
    RDKit::SubstructMatchParameters params;
    params.useChirality = false;
    params.extraFinalCheck = &already_matched;
    auto matches = RDKit::SubstructMatch(atomistic_mol, *query, params);

    for (const auto& match : matches) {
        std::vector<int> monomer;
        auto monomer_idx = monomers.size();
        for (const auto& [query_idx, atom_idx] : match) {
            auto atom = atomistic_mol.getAtomWithIdx(atom_idx);
            if (atom->hasProp(MONOMER_IDX1) && atom->hasProp(MONOMER_IDX2)) {
                // This shouldn't happen, sanity check
                throw std::runtime_error(fmt::format(
                    "Atom {} belongs to more than 2 monomers", atom->getIdx()));
            }

            if (atom->hasProp(MONOMER_IDX1)) {
                atom->setProp<unsigned int>(MONOMER_IDX2, monomer_idx);
                Linkage link = {atom->getProp<unsigned int>(MONOMER_IDX1),
                                static_cast<unsigned int>(monomer_idx),
                                atom->getProp<unsigned int>(ATTACH_FROM),
                                attch_map[query_idx]};
                // this is a linkage between two monomers
                linkages.push_back(link);
            } else {
                atom->setProp<unsigned int>(MONOMER_IDX1, monomer_idx);
                if (attch_map[query_idx] != NO_ATTACHMENT) {
                    atom->setProp(ATTACH_FROM, attch_map[query_idx]);
                }
            }
            monomer.push_back(atom_idx);

            // if there is a side chain, the attachment point will be at the
            // SIDECHAIN_IDX and will be indicated by the presence of the atom
            // map number. For now, assume there is a single side chain per
            // monomer
            if (query_idx == SIDECHAIN_IDX &&
                query->getAtomWithIdx(query_idx)->hasProp(
                    RDKit::common_properties::molAtomMapNumber)) {
                attch_pts[monomers.size()] = atom_idx;
            }
        }
        monomers.push_back(monomer);
    }
}

void add_sidechain_to_monomer(const RDKit::ROMol& atomistic_mol,
                              std::vector<int>& monomer,
                              unsigned int monomer_idx,
                              unsigned int attch_at_idx)
{

    // BFS but use MONOMER_IDX as visited marker
    std::queue<unsigned int> q;
    q.push(attch_at_idx);
    while (!q.empty()) {
        auto at_idx = q.front();
        q.pop();
        auto at = atomistic_mol.getAtomWithIdx(at_idx);
        if (!at->hasProp(MONOMER_IDX1)) {
            at->setProp<unsigned int>(MONOMER_IDX1, monomer_idx);
            monomer.push_back(at_idx);
        }
        for (const auto& nbr : atomistic_mol.atomNeighbors(at)) {
            if (!nbr->hasProp(MONOMER_IDX1)) {
                q.push(nbr->getIdx());
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
 * - firstMonomerIdx: index of the first monomer given atom belongs to
 * - secondMonomerIdx: index of the second monomer given atom belongs to (this
 * is optional, means there is a bond between two monomers)
 *
 * Returns a list monomer index sets, which represents the atom indices of each
 * monomer.
 */
void identify_monomers(RDKit::ROMol& atomistic_mol,
                       std::vector<std::vector<int>>& monomers,
                       std::vector<Linkage>& linkages)
{
    // Approach for identifying monomers:
    // 1. Find all matches with SMARTS queries for amino acids (TODO: Nucleic
    // acids & CHEM)
    // 2. Add side chains to generic matches based on attachment points
    // 3. Identify and group any remaining atoms into 'unclassified' monomers.
    // These are grouped
    //    by connectivity.
    std::unordered_map<unsigned int, unsigned int> attch_pts;
    add_matches_to_monomers(GENERIC_AMINO_ACID_QUERY, atomistic_mol, monomers,
                            attch_pts, linkages);
    add_matches_to_monomers(GLYCINE_AMINO_ACID_QUERY, atomistic_mol, monomers,
                            attch_pts, linkages);
    // TODO: nucleic acids and CHEM monomers

    // now, add sidechains onto each monomer
    for (size_t monomer_idx = 0; monomer_idx < monomers.size(); ++monomer_idx) {
        if (attch_pts.find(monomer_idx) != attch_pts.end()) {
            // there is a sidechain to add!
            add_sidechain_to_monomer(atomistic_mol, monomers[monomer_idx],
                                     monomer_idx, attch_pts[monomer_idx]);
        }
    }
}

void build_cg_mol(const RDKit::ROMol& atomistic_mol,
                  std::vector<std::vector<int>>& monomers,
                  boost::shared_ptr<RDKit::RWMol> cg_mol,
                  std::vector<Linkage>& linkages)
{
    // Start with all atoms in a single peptide chain
    cg_mol->setProp<bool>(HELM_MODEL, true);

    constexpr bool isomeric_smiles = false;
    int residue_num = 1;
    for (const auto& monomer : monomers) {
        auto monomer_smiles = RDKit::MolFragmentToSmiles(
            atomistic_mol, monomer, nullptr, nullptr, nullptr, isomeric_smiles);
        // We have to roundtrip to canonicalize smiles -- see RDKit issue #7214
        std::unique_ptr<RDKit::RWMol> canon_mol(
            RDKit::SmilesToMol(monomer_smiles));
        monomer_smiles = RDKit::MolToSmiles(*canon_mol);

        // If the monomer is a known amino acid, use the 1-letter code
        // NOTE: Setting the smilesSymbol is temporary & for testing purposes --
        // I think we'd actually want to set the atomLabel here
        if (amino_acids.find(monomer_smiles) != amino_acids.end()) {
            add_monomer(*cg_mol, amino_acids.at(monomer_smiles), residue_num,
                        "PEPTIDE1");
        } else {
            add_monomer(*cg_mol, monomer_smiles, residue_num, "PEPTIDE1",
                        MonomerType::SMILES);
        }
        ++residue_num;

        // TODO: Check for known nucleic acids and CHEM monomers
    }

    for (const auto& link : linkages) {
        // TODO: Non forward linkages
        add_connection(*cg_mol, link.monomer_idx1, link.monomer_idx2,
                       link.to_string());
    }
}

} // unnamed namespace

boost::shared_ptr<RDKit::RWMol> atomistic_to_cg(const RDKit::ROMol& mol)
{
    // Work on copy, for now
    RDKit::ROMol atomistic_mol(mol);
    std::vector<std::vector<int>> monomers;
    std::vector<Linkage> linkages;
    identify_monomers(atomistic_mol, monomers, linkages);

    boost::shared_ptr<RDKit::RWMol> cg_mol = boost::make_shared<RDKit::RWMol>();
    build_cg_mol(atomistic_mol, monomers, cg_mol, linkages);
    assign_chains(*cg_mol);

    // TODO
    // Now that we have the CG mol, we need to set the properties needed by the
    // HELM writer and other functions that work with CG mols created by the
    // HELM parser. I think this will include a few steps
    // 1. Break the CG Mol into polymers -- this would be done by connectivity
    // and monomer type (peptide, rna, dna, chem)
    // 2. Insure that the linkage information is correct -- backbone vs not, etc
    // 3. Set the polymers as substance groups on the molecule, and set
    // monomer-specific properties
    // 4. Maybe: make sure CG monomer indices are in connectivity order

    return cg_mol;
}

std::vector<std::vector<int>> get_monomers(const RDKit::ROMol& mol)
{
    RDKit::ROMol atomistic_mol(mol);
    std::vector<std::vector<int>> monomers;
    std::vector<Linkage> linkages;
    identify_monomers(atomistic_mol, monomers, linkages);
    return monomers;
}

} // namespace rdkit_extensions
} // namespace schrodinger
