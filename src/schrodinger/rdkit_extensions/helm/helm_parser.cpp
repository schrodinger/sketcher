#include "schrodinger/rdkit_extensions/helm/helm_parser.h"

#include <sstream>
#include <stdexcept>
#include <memory> // shared_ptr, make_shared
#include <string>

#include <GraphMol/Atom.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/helm/constants.h"

helm::HelmParser::HelmParser() :
    m_current_chain_id(1),
    m_current_residue_number(0),
    m_previous_residue_index(std::numeric_limits<unsigned int>::max()),
    m_mol(std::make_shared<::RDKit::RWMol>()){};

std::shared_ptr<::RDKit::RWMol>
helm::HelmParser::parse(const std::string& input_helm)
{
    std::istringstream stream(input_helm);
    helm::TokenScanner scanner(&stream);
    helm::TokenParser parser(scanner, *this);
    // success return 0
    if (parser.parse() != 0) {
        std::stringstream ss;
        ss << "Parsing failed around position " << m_num_chars_processed
           << ":\n";
        ss << input_helm << '\n';
        for (unsigned int i = 1; i < m_num_chars_processed; ++i) {
            ss << '_';
        }
        ss << "^";
        throw std::invalid_argument(ss.str());
    }
    return m_mol;
}

/*
 * Each atom should encode residue-level properties as follows:
 *   i. the pdb residue properties should be stored on the MonomerInfo object.
 *  ii. The monomer id should be set as the atom symbol so RDKit image
 *      generation apis can easily access that information.
 */
[[nodiscard]] static std::unique_ptr<::RDKit::Atom>
get_coarse_grain_atom_from_monomer_id(const std::string_view monomer_id,
                                      const int& residue_number,
                                      const unsigned int chain_id)
{
    auto* residue_info = new RDKit::AtomPDBResidueInfo();
    residue_info->setResidueNumber(residue_number);
    residue_info->setChainId(std::to_string(chain_id));

    const std::string residue_name{monomer_id};
    residue_info->setResidueName(residue_name);

    auto atom = std::make_unique<::RDKit::Atom>();
    atom->setMonomerInfo(residue_info);
    atom->setProp(ATOM_LABEL_PROP_NAME, residue_name);
    return atom;
}

void helm::HelmParser::add_monomer(const std::string_view monomer_id)
{
    auto coarse_grain_atom = get_coarse_grain_atom_from_monomer_id(
        monomer_id, ++m_current_residue_number, m_current_chain_id);

    const auto& current_atom_index =
        m_mol->addAtom(coarse_grain_atom.release(), true, true);
    if (m_previous_residue_index != std::numeric_limits<unsigned int>::max()) {
        m_mol->addBond(m_previous_residue_index, current_atom_index);
    }
    m_previous_residue_index = current_atom_index;
}

void helm::HelmParser::add_polymer(const std::string_view polymer_id)
{
    ++m_current_chain_id;
    m_current_residue_number = 0;
    m_previous_residue_index = std::numeric_limits<unsigned int>::max();
}

void helm::HelmParser::saveErrorInformation(
    const unsigned int num_chars_processed)
{
    m_num_chars_processed = num_chars_processed;
}
