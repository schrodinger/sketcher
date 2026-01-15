#include "schrodinger/sketcher/rdkit/monomer_connectors.h"

#include <string>

#include <rdkit/GraphMol/Bond.h>

#include "schrodinger/rdkit_extensions/helm.h"

#include <iostream>
namespace schrodinger
{
namespace sketcher
{

bool contains_two_monomer_linkages(const RDKit::Bond* bond)
{
    std::string linkage, custom_linkage;
    bond->getPropIfPresent(LINKAGE, linkage);
    bool custom_linkage_exists =
        bond->getPropIfPresent(CUSTOM_BOND, custom_linkage);
    return custom_linkage_exists && custom_linkage != linkage;
}

} // namespace sketcher
} // namespace schrodinger
