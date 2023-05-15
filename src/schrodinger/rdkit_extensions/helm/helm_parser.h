#pragma once

#include <string>
#include <string_view>
#include <cstddef>
#include <istream>

#include <GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/helm/token_scanner.h"
#include "schrodinger/rdkit_extensions/helm/generated/helm_parser.tab.hh"

namespace helm
{

/*
 * A helper class to interact with a yacc/flex generated parser (TokenParser).
 * The associated TokenParser instance will use an instance of this class to
 * construct the biomolecule while it parser the HELM string.
 */
class HelmParser
{
  public:
    HelmParser();

    std::shared_ptr<::RDKit::RWMol> parse(const std::string&);
    void saveErrorInformation(const unsigned int num_chars_processed);

    void add_monomer(const std::string_view monomer_id);
    void add_polymer(const std::string_view polymer_id);

  private:
    unsigned int m_current_chain_id;
    unsigned int m_current_residue_number;
    unsigned int m_previous_residue_index;
    std::shared_ptr<::RDKit::RWMol> m_mol;

    unsigned int m_num_chars_processed;
};

} /* end namespace helm */
