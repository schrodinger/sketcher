/* -------------------------------------------------------------------------
 * Implements schrodinger::rdkit_extensions:: miscellaneous mol operations
 mol conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */
#include "schrodinger/rdkit_extensions/molops.h"

#include <rdkit/GraphMol/MolOps.h>

namespace schrodinger
{
namespace rdkit_extensions
{

const std::string ATTACHMENT_POINT_LABEL_PREFIX{"_AP"};
const std::string R_GROUP_LABEL_PREFIX{"_R"};
const int DUMMY_ATOMIC_NUMBER{0};

} // namespace rdkit_extensions
} // namespace schrodinger
