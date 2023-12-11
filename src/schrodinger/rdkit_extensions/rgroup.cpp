/* -------------------------------------------------------------------------
 * Implements schrodinger::rdkit_extensions:: rgroup querying and creation apis
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#include "schrodinger/rdkit_extensions/rgroup.h"

#include <memory>
#include <optional>
#include <stdexcept>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/QueryAtom.h>

#include "schrodinger/rdkit_extensions/constants.h"

namespace schrodinger
{
namespace rdkit_extensions
{

bool is_attachment_point_dummy(const RDKit::Atom& atom)
{
    std::string label;
    return atom.getAtomicNum() == 0 && atom.getTotalDegree() == 1 &&
           atom.getPropIfPresent(RDKit::common_properties::atomLabel, label) &&
           label.find(ATTACHMENT_POINT_LABEL_PREFIX) == 0;
}

[[nodiscard]] std::shared_ptr<RDKit::Atom>
make_new_r_group(const unsigned int r_group_num)
{
    if (r_group_num == 0) {
        throw std::invalid_argument("Rgroups can't have an index of '0'.");
    }

    auto atom = std::make_shared<RDKit::QueryAtom>();
    atom->setAtomicNum(DUMMY_ATOMIC_NUMBER);
    // This prevents the atom from being shown in SMARTS as [#0]
    atom->setQuery(RDKit::makeAtomNullQuery());
    // This prevents the RGroup from having a VAL in MOL file output
    atom->setNoImplicit(false);

    // we need to set the atom and dummy labels in the formats _RX and RX
    // respectively
    auto rlabel = R_GROUP_LABEL_PREFIX + std::to_string(r_group_num);
    auto dlabel = rlabel.substr(1);
    atom->setProp(RDKit::common_properties::atomLabel, rlabel);
    atom->setProp(RDKit::common_properties::dummyLabel, dlabel);

    atom->setProp(RDKit::common_properties::_MolFileRLabel, r_group_num);
    atom->setIsotope(r_group_num);
    return atom;
}

[[nodiscard]] std::optional<unsigned int>
get_r_group_number(const RDKit::Atom* const atom)
{
    unsigned int r_group_num = 0;
    if (!(atom->getAtomicNum() == DUMMY_ATOMIC_NUMBER &&
          atom->getPropIfPresent(RDKit::common_properties::_MolFileRLabel,
                                 r_group_num))) {
        return std::nullopt;
    }

    if (r_group_num == 0) {
        throw std::invalid_argument("Rgroups can't have an index of '0'.");
    }

    return std::make_optional<unsigned int>(r_group_num);
}
} // namespace rdkit_extensions
} // namespace schrodinger
