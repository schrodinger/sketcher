#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "schrodinger/sketcher/definitions.h"

class QButtonGroup;

namespace schrodinger
{
namespace rdkit_extensions
{
struct MonomerInfo;
}

namespace sketcher
{

class ModularPopup;

/**
 * Populate `popup` with a horizontal row of QToolButtons: ID 0 is the
 * standard monomer (text = standard_symbol, tooltip = standard_name); IDs
 * 1+ are the analogs. Each button's object name is set to
 * "<object_name_prefix>_<symbol>_btn", and the id->symbol mapping is
 * written into `id_to_symbol`.
 *
 * @return The QButtonGroup containing the buttons. The caller must pass
 * it to ModularPopup::setButtonGroup() to finish initialization.
 */
SKETCHER_API QButtonGroup* build_monomer_symbol_buttons(
    ModularPopup* popup, const std::string& object_name_prefix,
    const std::string& standard_symbol, const std::string& standard_name,
    const std::vector<rdkit_extensions::MonomerInfo>& analogs,
    std::unordered_map<int, std::string>& id_to_symbol);

} // namespace sketcher
} // namespace schrodinger
