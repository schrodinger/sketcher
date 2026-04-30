#include "schrodinger/sketcher/widget/monomer_symbol_popup_utils.h"

#include <QButtonGroup>
#include <QHBoxLayout>
#include <QToolButton>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/widget/modular_popup.h"

namespace schrodinger
{
namespace sketcher
{

QButtonGroup* build_monomer_symbol_buttons(
    ModularPopup* popup, const std::string& object_name_prefix,
    const std::string& standard_symbol, const std::string& standard_name,
    const std::vector<rdkit_extensions::MonomerInfo>& analogs,
    std::unordered_map<int, std::string>& id_to_symbol)
{
    auto* layout = new QHBoxLayout(popup);
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
        layout->addWidget(btn);
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
