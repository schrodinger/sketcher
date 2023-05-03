#include "schrodinger/sketcher/rdkit/periodic_table.h"

#include <GraphMol/PeriodicTable.h>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/trim.hpp>

namespace schrodinger
{
namespace sketcher
{

std::string atomic_number_to_symbol(unsigned int atomic_number)
{
    const auto* table = RDKit::PeriodicTable::getTable();
    return table->getElementSymbol(atomic_number);
}

std::string atomic_number_to_name(unsigned int atomic_number)
{
    const auto* table = RDKit::PeriodicTable::getTable();
    return table->getElementName(atomic_number);
}

int symbol_to_atomic_number(std::string element_symbol)
{
    // Temporarily silence rdkit logging
    RDLog::LogStateSetter silence_rdkit_logging;

    // Properly format as a chemical symbol
    boost::trim(element_symbol);
    boost::algorithm::to_lower(element_symbol);
    if (!element_symbol.empty()) {
        element_symbol[0] = std::toupper(element_symbol[0]);
    }

    const auto* table = RDKit::PeriodicTable::getTable();
    return table->getAtomicNumber(element_symbol);
}

bool is_atomic_number(unsigned int atomic_number)
{
    const unsigned int HEAVIEST_ELEMENT_NUMBER = 118;
    return (atomic_number > 0 && atomic_number <= HEAVIEST_ELEMENT_NUMBER);
}

} // namespace sketcher
} // namespace schrodinger