#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"

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
class SketcherModel;
enum class MonomerToolType;

/**
 * Find the popup entry matching the active monomer tool on this page.
 */
SKETCHER_API int get_monomer_symbol_button_id(
    const SketcherModel* model, MonomerToolType tool_type,
    const std::unordered_map<int, rdkit_extensions::MonomerID>& id_to_monomer);

/**
 * Populate `popup` with a row-major grid of QToolButtons: ID 0 is the
 * standard monomer (text = standard_symbol, tooltip = standard_name); IDs
 * 1+ are the analogs. Each button's object name is set to
 * "<object_name_prefix>_<symbol>_btn", and the id->monomer mapping is
 * written into `id_to_monomer`.
 * An empty standard_symbol omits the standard monomer; analog IDs start at 0.
 *
 * @return The QButtonGroup containing the buttons. The caller must pass
 * it to ModularPopup::setButtonGroup() to finish initialization.
 */
SKETCHER_API QButtonGroup* build_monomer_symbol_buttons(
    ModularPopup* popup, const std::string& object_name_prefix,
    const std::string& standard_symbol, const std::string& standard_name,
    const std::vector<rdkit_extensions::MonomerInfo>& analogs,
    std::unordered_map<int, rdkit_extensions::MonomerID>& id_to_monomer,
    rdkit_extensions::ChainType default_chain_type);

} // namespace sketcher
} // namespace schrodinger
