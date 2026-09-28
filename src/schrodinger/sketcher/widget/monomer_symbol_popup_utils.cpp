#include "schrodinger/sketcher/widget/monomer_symbol_popup_utils.h"

#include <QButtonGroup>
#include <QGridLayout>
#include <QToolButton>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/widget/modular_popup.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * Determine the number of columns we should use for a monomer popup containing
 * the specified number of monomers. This function tries to keep the popup
 * roughly square-ish while using the same number of columns for popups with
 * similar numbers of monomers.
 */
static size_t get_num_columns(const size_t num_monomers)
{
    const size_t MIN_NUMBER_OF_COLUMNS = 4;
    const float MAX_ROWS_TO_COLUMNS_RATIO = 1.25;

    size_t num_columns = MIN_NUMBER_OF_COLUMNS;
    while (true) {
        if (MAX_ROWS_TO_COLUMNS_RATIO * num_columns * num_columns >=
            num_monomers) {
            // the number of monomers will fit within num_columns using at most
            // `MAX_ROWS_TO_COLUMNS_RATIO * num_columns` rows, so use this
            // number of columns
            return num_columns;
        }
        num_columns *= 2;
    }
}

QButtonGroup* build_monomer_symbol_buttons(
    ModularPopup* popup, const std::string& object_name_prefix,
    const std::string& standard_symbol, const std::string& standard_name,
    const std::vector<rdkit_extensions::MonomerInfo>& analogs,
    std::unordered_map<int, std::string>& id_to_symbol)
{
    const auto num_monomers = analogs.size() + 1;
    const auto num_columns = get_num_columns(num_monomers);
    auto* layout = new QGridLayout(popup);
    layout->setContentsMargins(2, 2, 2, 2);
    layout->setSpacing(0);

    auto* group = new QButtonGroup(popup);

    int id = 0;
    auto make_button = [&](const std::string& symbol, const std::string& name) {
        auto* btn = new QToolButton(popup);
        btn->setText(QString::fromStdString(symbol));
        btn->setToolTip(QString::fromStdString(name));
        btn->setCheckable(true);
        btn->setMinimumSize(32, 30);
        btn->setMaximumSize(32, 32);
        btn->setObjectName(
            QString::fromStdString(object_name_prefix + "_" + symbol + "_btn"));
        btn->setStyleSheet(symbol.size() >= COMPACT_STYLE_MIN_LENGTH
                               ? ATOM_ELEMENT_OR_MONOMER_COMPACT_STYLE
                               : ATOM_ELEMENT_OR_MONOMER_STYLE);
        layout->addWidget(btn, id / num_columns, id % num_columns);
        group->addButton(btn);
        id_to_symbol[id] = symbol;
        ++id;
    };

    make_button(standard_symbol, standard_name);
    for (const auto& analog : analogs) {
        make_button(analog.symbol.value_or(""), analog.name.value_or(""));
    }
    return group;
}

} // namespace sketcher
} // namespace schrodinger
