#pragma once

#include <string>
#include <vector>

#include "schrodinger/rdkit_extensions/helm/helm_parser.h"

namespace helm
{

/*
 * Api to validate the parsed HELM information. We should try to catch
 * invalid information and unsupported features towards providing more
 * meaningful error messages to users.
 *
 * Some examples include:
 *      * duplicate polymer ids
 *      * duplicate polymer group ids
 *      * BLOB or unknown/wildcard residues having defined rgroups
 *
 * See implementation for a more exhaustive list of things we don't support.
 */
bool validate_parsed_info(const helm_info& parsed_info,
                          std::vector<std::string>& errors,
                          const std::string_view& input);
} // namespace helm
