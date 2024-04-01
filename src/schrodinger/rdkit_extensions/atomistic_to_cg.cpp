#include "schrodinger/rdkit_extensions/atomistic_to_cg.h"

#include <queue>
#include <unordered_map>
#include <memory>

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

// Assume that all of these queries are completely disjoint
// On these SMARTS queries, map no 1 exists on side chain
constexpr int SIDECHAIN_IDX = 2;
const std::string GENERIC_AMINO_ACID_QUERY =
    "[NX3,NX4+][CX4H]([*:1])[CX3](=[OX1])[O,N]";
const std::string GLYCINE_AMINO_ACID_QUERY =
    "N[CX4H2][CX3](=[OX1])[O,N]"; // no side chain

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

/*
 * Function that takes the SMARTS query and atomistic molecule and adds the atom
 * indices of the matches to the monomers vector
 *
 */
void add_matches_to_monomers(
    const std::string& smarts_query, RDKit::ROMol& atomistic_mol,
    std::vector<std::vector<int>>& monomers,
    std::unordered_map<unsigned int, unsigned int>& attch_pts,
    std::vector<std::pair<unsigned int, unsigned int>>& linkages)
{
    // TODO: maybe the query mol itself should be const somewhere
    std::unique_ptr<RDKit::RWMol> query =
        std::make_unique<RDKit::RWMol>(*RDKit::SmartsToMol(smarts_query));
    auto matches = RDKit::SubstructMatch(atomistic_mol, *query);

    for (const auto& match : matches) {
        std::vector<int> monomer;
        auto monomer_idx = monomers.size();
        for (const auto& [query_idx, atom_idx] : match) {
            auto atom = atomistic_mol.getAtomWithIdx(atom_idx);
            if (atom->hasProp(MONOMER_IDX_PROP1) &&
                atom->hasProp(MONOMER_IDX_PROP2)) {
                throw std::runtime_error(
                    "Atom belongs to more than 2 monomers");
            } else if (atom->hasProp(MONOMER_IDX_PROP1)) {
                atom->setProp<unsigned int>(MONOMER_IDX_PROP2, monomer_idx);
                // this is a linkage between two monomers
                linkages.push_back(
                    {atom->getProp<unsigned int>(MONOMER_IDX_PROP1),
                     monomer_idx});
            } else {
                atom->setProp<unsigned int>(MONOMER_IDX_PROP1, monomer_idx);
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

    // BFS but use MONOMER_IDX_PROP as visited marker
    std::queue<unsigned int> q;
    q.push(attch_at_idx);
    while (!q.empty()) {
        auto at_idx = q.front();
        q.pop();
        auto at = atomistic_mol.getAtomWithIdx(at_idx);
        if (!at->hasProp(MONOMER_IDX_PROP1)) {
            at->setProp<unsigned int>(MONOMER_IDX_PROP1, monomer_idx);
            monomer.push_back(at_idx);
        }
        for (const auto& nbr : atomistic_mol.atomNeighbors(at)) {
            if (!nbr->hasProp(MONOMER_IDX_PROP1)) {
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
void identify_monomers(
    RDKit::ROMol& atomistic_mol, std::vector<std::vector<int>>& monomers,
    std::vector<std::pair<unsigned int, unsigned int>>& linkages)
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

    // TODO: gather extraneous atoms into unclassified monomers
}

void build_cg_mol(const RDKit::ROMol& atomistic_mol,
                  std::vector<std::vector<int>>& monomers,
                  boost::shared_ptr<RDKit::RWMol> cg_mol,
                  std::vector<std::pair<unsigned int, unsigned int>>& linkages)
{
    // Start with all atoms in a single peptide chain
    cg_mol->setProp<bool>(HELM_MODEL, true);
    auto& chain = add_chain(*cg_mol, ChainType::PEPTIDE);

    constexpr bool isomeric_smiles = false;
    for (const auto& monomer : monomers) {
        auto monomer_smiles = RDKit::MolFragmentToSmiles(
            atomistic_mol, monomer, nullptr, nullptr, nullptr, isomeric_smiles);
        // roundtrip to canonicalize smiles
        monomer_smiles =
            RDKit::MolToSmiles(*RDKit::SmilesToMol(monomer_smiles));

        // If the monomer is a known amino acid, use the 1-letter code
        // NOTE: Setting the smilesSymbol is temporary & for testing purposes --
        // I think we'd actually want to set the atomLabel here
        if (amino_acids.find(monomer_smiles) != amino_acids.end()) {
            add_monomer(chain, amino_acids.at(monomer_smiles));
        } else {
            add_monomer(chain, monomer_smiles, MonomerType::SMILES);
        }

        // TODO: Check for known nucleic acids and CHEM monomers
    }

    for (const auto& [monomer_idx1, monomer_idx2] : linkages) {
        // TODO: Non forward linkages
        add_connection(*cg_mol, monomer_idx1, monomer_idx2);
    }
}

} // unnamed namespace

boost::shared_ptr<RDKit::RWMol> atomistic_to_cg(const RDKit::ROMol& mol)
{
    // Work on copy, for now
    RDKit::ROMol atomistic_mol(mol);
    std::vector<std::vector<int>> monomers;
    std::vector<std::pair<unsigned int, unsigned int>> linkages;
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
    std::vector<std::pair<unsigned int, unsigned int>> linkages;
    identify_monomers(atomistic_mol, monomers, linkages);
    return monomers;
}

} // namespace rdkit_extensions
} // namespace schrodinger
