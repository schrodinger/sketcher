#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/modular_popup.h"

namespace schrodinger
{
namespace rdkit_extensions
{
struct MonomerInfo;
}

namespace sketcher
{

/**
 * Popup showing non-natural analogs for a standard nucleic acid base button.
 *
 * The first entry (ID 0) is the standard base; subsequent entries are
 * non-natural analogs (e.g. 5-methylcytosine) sourced from the monomer
 * database.
 */
class SKETCHER_API NucleicAcidSymbolPopup : public ModularPopup
{
  public:
    NucleicAcidSymbolPopup(
        const std::string& standard_symbol, const std::string& standard_name,
        const std::vector<rdkit_extensions::MonomerInfo>& analogs,
        QWidget* parent = nullptr);

    /**
     * @return The HELM symbol for the given button ID
     */
    QString getSymbolForId(int id) const;

  protected:
    void generateButtonPackets() override;
    int getButtonIDToCheck() override;

  private:
    std::unordered_map<int, std::string> m_id_to_symbol;
};

} // namespace sketcher
} // namespace schrodinger
