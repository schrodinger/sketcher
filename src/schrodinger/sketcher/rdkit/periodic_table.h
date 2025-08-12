
#pragma once

#include <string>

#include "schrodinger/sketcher/definitions.h"

namespace RDKit
{
class Atom;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * @param atomic_number an atomic number
 * @return the element symbol associated with the specified atomic number
 * @throw Invariant if an unrecognized atomic number is specified
 */
SKETCHER_API std::string atomic_number_to_symbol(unsigned int atomic_number);

/**
 * @param atomic_number an atomic number
 * @return the full element name associated with the specified atomic number
 * @throw Invariant if an unrecognized atomic number is specified
 */
SKETCHER_API std::string atomic_number_to_name(unsigned int atomic_number);

/**
 * @param element_symbol a element symbol
 * @return the atomic number associated with the specified element symbol.
 * @throw Invariant if an unrecognized element symbol is specified
 */
SKETCHER_API int symbol_to_atomic_number(std::string element_symbol);

/**
 * @param atomic_number an atomic number
 * @return whether or not the given number corresponds to a valid element
 */
SKETCHER_API bool is_atomic_number(int atomic_number);

} // namespace sketcher
} // namespace schrodinger