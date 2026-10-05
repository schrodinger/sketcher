#include "schrodinger/sketcher/widget/monomer_symbol_popup_utils.h"

#include <QButtonGroup>
#include <QGridLayout>
#include <QScrollArea>
#include <QScrollBar>
#include <QToolButton>
#include <QVBoxLayout>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/widget/modular_popup.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

static constexpr int MAX_BUTTON_HEIGHT = 32;
/// The maximum number of rows visible at once before the user has to scroll
static constexpr int MAX_ROWS = 10;

int get_monomer_symbol_button_id(
    const SketcherModel* model, MonomerToolType tool_type,
    const std::unordered_map<int, rdkit_extensions::MonomerID>& id_to_monomer)
{
    if (model == nullptr || model->getMonomerToolType() != tool_type) {
        return -1;
    }
    rdkit_extensions::MonomerID selection;
    if (model->getDrawTool() == DrawTool::MONOMER_DB_MONOMER) {
        selection = model->getValue(ModelKey::MONOMER_DB_MONOMER)
                        .value<rdkit_extensions::MonomerID>();
    } else if (model->getDrawTool() == DrawTool::MONOMER) {
        if (tool_type == MonomerToolType::AMINO_ACID) {
            selection = {model->getValueString(ModelKey::AMINO_ACID_SYMBOL)
                             .toStdString(),
                         rdkit_extensions::ChainType::PEPTIDE};
        } else {
            selection = {model->getValue(ModelKey::NUCLEIC_ACID_SYMBOL)
                             .value<NucleicAcidMutation>()
                             .symbol.toStdString(),
                         rdkit_extensions::ChainType::RNA};
        }
    } else {
        return -1;
    }
    for (const auto& [id, monomer] : id_to_monomer) {
        if (monomer == selection) {
            return id;
        }
    }
    return -1;
}

QButtonGroup* build_monomer_symbol_buttons(
    ModularPopup* popup, const std::string& object_name_prefix,
    const std::string& standard_symbol, const std::string& standard_name,
    const std::vector<rdkit_extensions::MonomerInfo>& analogs,
    std::unordered_map<int, rdkit_extensions::MonomerID>& id_to_monomer,
    rdkit_extensions::ChainType default_chain_type)
{
    const auto num_monomers = analogs.size() + !standard_symbol.empty();
    const auto num_columns = num_monomers <= 20 ? 4 : 8;
    const bool needs_scroll = num_monomers > 80;
    auto* button_widget = needs_scroll ? new QWidget(popup) : popup;
    auto* layout = new QGridLayout(button_widget);
    layout->setContentsMargins(2, 2, 2, 2);
    layout->setSpacing(0);

    auto* group = new QButtonGroup(popup);

    int id = 0;
    auto make_button = [&](const std::string& symbol, const std::string& name,
                           rdkit_extensions::ChainType chain_type) {
        auto* btn = new QToolButton(button_widget);
        btn->setText(QString::fromStdString(symbol));
        btn->setToolTip(QString::fromStdString(name));
        btn->setCheckable(true);
        btn->setMinimumSize(32, 30);
        btn->setMaximumSize(32, MAX_BUTTON_HEIGHT);
        btn->setObjectName(
            QString::fromStdString(object_name_prefix + "_" + symbol + "_btn"));
        btn->setStyleSheet(symbol.size() >= COMPACT_STYLE_MIN_LENGTH
                               ? ATOM_ELEMENT_OR_MONOMER_COMPACT_STYLE
                               : ATOM_ELEMENT_OR_MONOMER_STYLE);
        layout->addWidget(btn, id / num_columns, id % num_columns);
        group->addButton(btn);
        id_to_monomer[id] = {symbol, chain_type};
        ++id;
    };

    if (!standard_symbol.empty()) {
        make_button(standard_symbol, standard_name, default_chain_type);
    }
    for (const auto& analog : analogs) {
        make_button(analog.symbol.value_or(""), analog.name.value_or(""),
                    analog.polymer_type == "CHEM"
                        ? rdkit_extensions::ChainType::CHEM
                        : default_chain_type);
    }
    if (needs_scroll) {
        // Keep the grid at its natural size and show at most ten full rows.
        // Reserve space for the scrollbar so no columns are clipped.
        auto* scroll_area = new QScrollArea(popup);
        scroll_area->setFrameShape(QFrame::NoFrame);
        scroll_area->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        scroll_area->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
        layout->setSizeConstraint(QLayout::SetFixedSize);
        scroll_area->setWidget(button_widget);
        const auto margins = layout->contentsMargins();
        scroll_area->setFixedSize(
            layout->sizeHint().width() +
                scroll_area->verticalScrollBar()->sizeHint().width(),
            MAX_ROWS * MAX_BUTTON_HEIGHT + margins.top() + margins.bottom());
        auto* popup_layout = new QVBoxLayout(popup);
        popup_layout->setContentsMargins(0, 0, 0, 0);
        popup_layout->addWidget(scroll_area);
    }
    return group;
}

} // namespace sketcher
} // namespace schrodinger
