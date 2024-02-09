/* -------------------------------------------------------------------------
 * Declares schrodinger::rdkit_extensions:: atomistic ROMol -> CG mol conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */
#pragma once

#include <boost/shared_ptr.hpp>

#include "schrodinger/rdkit_extensions/definitions.h"

const std::string MONOMER_IDX_PROP1{"monomerIndex1"};
const std::string MONOMER_IDX_PROP2{"monomerIndex2"};

namespace RDKit
{
class ROMol;
class RWMol;
} // namespace RDKit

namespace schrodinger
{
namespace rdkit_extensions
{

/**
 * Converts an atomistic ROMol into a CG ROMol that can
 * be written to HELM format using the HELM writer.
 *
 * @param atomistic_mol Atomistic molecule to convert to CG
 * @return CG molecule
 */
RDKIT_EXTENSIONS_API boost::shared_ptr<RDKit::RWMol>
atomistic_to_cg(const RDKit::ROMol& atomistic_mol);

/**
 * Identify monomers within an atomistic molecule
 *
 * For testing purposes.
 */
RDKIT_EXTENSIONS_API std::vector<std::vector<int>>
get_monomers(const RDKit::ROMol& mol);

} // namespace rdkit_extensions
} // namespace schrodinger
