#include "schrodinger/sketcher/rdkit/periodic_table.h"

#include <GraphMol/GraphMol.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/SanitException.h>

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/trim.hpp>

namespace schrodinger
{
namespace sketcher
{

const int HEAVIEST_ELEMENT_NUMBER = 118;

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

bool is_atomic_number(int atomic_number)
{
    return (atomic_number > 0 && atomic_number <= HEAVIEST_ELEMENT_NUMBER);
}

bool has_valence_violation(const RDKit::Atom* atom)
{
    // Ignore dummy atoms, query atoms, or atoms attached to query bonds
    auto bonds = atom->getOwningMol().atomBonds(atom);
    auto is_query = [](auto b) { return b->hasQuery(); };
    if (atom->getAtomicNum() == 0 || atom->hasQuery() ||
        std::any_of(bonds.begin(), bonds.end(), is_query)) {
        return false;
    }

    try {
        bool strict = true;
        const_cast<RDKit::Atom*>(atom)->updatePropertyCache(strict);
    } catch (const RDKit::AtomValenceException&) {
        return true;
    }
    return false;
}

} // namespace sketcher
} // namespace schrodinger