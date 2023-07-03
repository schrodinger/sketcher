#include "schrodinger/sketcher/rdkit/periodic_table.h"

#include <GraphMol/GraphMol.h>
#include <GraphMol/PeriodicTable.h>
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
    auto atomic_number = atom->getAtomicNum();
    auto bonds = atom->getOwningMol().atomBonds(atom);
    auto is_query = [](auto b) { return b->hasQuery(); };
    if (atomic_number == 0 || atom->hasQuery() ||
        std::any_of(bonds.begin(), bonds.end(), is_query)) {
        return false;
    }

    // Account for charge via the isolobal analogy by modifying the atomic
    // number used to look up allowed valences
    atomic_number -= atom->getFormalCharge();
    // Because the lookup accounted for charge, we can ignore that contribution
    // to the current valence; however, still account for unpaired electrons
    int current_valence =
        atom->getTotalValence() + atom->getNumRadicalElectrons();

    // Special case proton, which is permissible
    if (atomic_number == 0 && current_valence == 0) {
        return false;
    }

    // If the isolobal analogy breaks down, assume an error
    if (atomic_number < 1 || atomic_number > HEAVIEST_ELEMENT_NUMBER) {
        return true;
    }

    const auto* table = RDKit::PeriodicTable::getTable();
    auto allowed_valence_list = table->getValenceList(atomic_number);
    if (allowed_valence_list == std::vector<int>({-1})) {
        return false; // only -1 indicates that any valence is permitted
    }

    return std::find(allowed_valence_list.begin(), allowed_valence_list.end(),
                     current_valence) == allowed_valence_list.end();
}

} // namespace sketcher
} // namespace schrodinger