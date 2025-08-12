#include "schrodinger/rdkit_extensions/dummy_atom.h"

#include <rdkit/GraphMol/QueryAtom.h>

#include "schrodinger/rdkit_extensions/constants.h"

namespace schrodinger
{
namespace rdkit_extensions
{

[[nodiscard]] std::shared_ptr<RDKit::Atom> create_dummy_atom()
{
    auto atom = std::make_shared<RDKit::QueryAtom>();
    atom->setAtomicNum(DUMMY_ATOMIC_NUMBER);
    // This prevents the atom from being shown in SMARTS as [#0]
    atom->setQuery(RDKit::makeAtomNullQuery());
    // This prevents the atom from having a VAL in MDL MOL file output
    atom->setNoImplicit(false);
    return atom;
}

} // namespace rdkit_extensions
} // namespace schrodinger
