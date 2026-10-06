#include "schrodinger/sketcher/dialog/custom_monomer_dialog.h"

#include <algorithm>
#include <map>
#include <unordered_set>

#include <QDialogButtonBox>
#include <QPalette>
#include <QPushButton>
#include <QResizeEvent>
#include <QScopedValueRollback>
#include <QStringList>

#include <stdexcept>

#include <rdkit/GraphMol/MolOps.h>

#include "schrodinger/rdkit_extensions/file_format.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/rdkit_extensions/rgroup.h"
#include "schrodinger/sketcher/dialog/message_box_dialog.h"
#include "schrodinger/sketcher/public_constants.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/sketcher_widget.h"
#include "schrodinger/sketcher/ui/ui_custom_monomer_dialog.h"
#include "schrodinger/sketcher/ui/ui_sketcher_widget.h"
#include "schrodinger/sketcher/widget/sketcher_side_bar.h"

using schrodinger::rdkit_extensions::ChainType;
using schrodinger::rdkit_extensions::Format;
using schrodinger::rdkit_extensions::get_r_group_number;

namespace schrodinger
{
namespace sketcher
{

/**
 * The width of the QDialog border added in schrodinger_livedesign.qss
 */
static constexpr int DIALOG_BORDER_WIDTH = 1;

/**
 * @return a list of all R-groups that appear more than once in the specified
 * molecule
 */
static std::vector<int> get_duplicate_r_groups(const RDKit::ROMol& mol)
{
    std::map<unsigned int, unsigned int> r_group_counts;
    for (const auto* atom : mol.atoms()) {
        if (const auto r_group_num = get_r_group_number(atom)) {
            ++r_group_counts[*r_group_num];
        }
    }

    std::vector<int> duplicate_r_groups;
    for (const auto& [r_group_num, count] : r_group_counts) {
        if (count > 1) {
            duplicate_r_groups.push_back(static_cast<int>(r_group_num));
        }
    }
    return duplicate_r_groups;
}

/**
 * @return formatted text listing all R-groups in the given list of R-groups
 */
static QString format_r_group_list(const std::vector<int>& r_group_numbers)
{
    if (r_group_numbers.empty()) {
        return {};
    }

    QStringList r_groups;
    for (const auto r_group_num : r_group_numbers) {
        r_groups.append("R" + QString::number(r_group_num));
    }

    if (r_groups.size() == 1) {
        return r_groups.front();
    }

    const auto last_r_group = r_groups.takeLast();
    const auto separator = r_groups.size() == 1 ? " " : ", ";
    return r_groups.join(", ") + separator + "and " + last_r_group;
}

static QString chain_type_display_name(const ChainType chain_type)
{
    switch (chain_type) {
        case ChainType::PEPTIDE:
            return "Peptide";
        case ChainType::RNA:
            return "Nucleic Acid";
        case ChainType::CHEM:
            return "Chem";
        default:
            throw std::invalid_argument(
                "Custom monomer dialog does not support this chain type");
    }
}

CustomMonomerDialog::CustomMonomerDialog(const ChainType chain_type,
                                         QWidget* parent) :
    ResizableModelDialog(parent),
    m_chain_type(chain_type)
{
    ui.reset(new Ui::CustomMonomerDialog());
    setupDialogUI(*ui);
    setStyleSheet(CUSTOM_MONOMER_DIALOG_STYLE);
    setWindowTitle("Sketch Custom " + chain_type_display_name(chain_type) +
                   " Monomer");
    ui->sketcher_widget->setInterfaceType(InterfaceType::ATOMISTIC);

#ifdef __EMSCRIPTEN__
    configureWasmTitleBar();
#endif
    qobject_cast<QVBoxLayout*>(layout())->setSpacing(0);
    ui->sketcher_widget->setInterfaceToggleVisible(false);

    // Remove the standard padding, but leave the dialog border exposed.
    m_dlg_layout->setContentsMargins(0, 0, 0, 0);
    layout()->setContentsMargins(DIALOG_BORDER_WIDTH, DIALOG_BORDER_WIDTH,
                                 DIALOG_BORDER_WIDTH, DIALOG_BORDER_WIDTH);

    // The full-width footer must not prevent shrinking into the compact layout.
    layout()->setSizeConstraint(QLayout::SetNoConstraint);
    ui->verticalLayout->setSizeConstraint(QLayout::SetNoConstraint);
    ui->verticalLayout->parentWidget()->setMinimumSize(0, 0);
    ensurePolished();

    // Preserve the footer's original background when it moves into the white
    // SketcherWidget, rather than inheriting its new parent's background.
    auto footer_palette = ui->button_bar->palette();
    footer_palette.setColor(QPalette::Window,
                            footer_palette.color(QPalette::Window));
    ui->button_bar->setPalette(footer_palette);
    ui->button_bar->setAutoFillBackground(true);

    m_layout_ready = true;
    updateButtonBarPlacement();

    connect(ui->sketcher_widget, &SketcherWidget::moleculeChanged, this,
            &CustomMonomerDialog::updateOkButton);
    connect(ui->sketcher_widget, &SketcherWidget::representationChanged, this,
            &CustomMonomerDialog::updateOkButton);
    updateOkButton();
}

CustomMonomerDialog::~CustomMonomerDialog() = default;

void CustomMonomerDialog::configureWasmTitleBar()
{
    if (m_title_bar == nullptr) {
        return;
    }
    // Retain the normal maximum height while allowing the title bar to shrink
    // into the space freed by hiding the interface toggle, accounting for the
    // top and bottom border margins.
    m_title_bar->setMinimumHeight(
        ui->sketcher_widget->getSideBar()->getInterfaceToggleHeight() -
        2 * DIALOG_BORDER_WIDTH);
    auto policy = m_title_bar->sizePolicy();
    policy.setVerticalPolicy(QSizePolicy::Preferred);
    m_title_bar->setSizePolicy(policy);
}

static QSize layout_item_size(const QLayoutItem* item, const bool minimum)
{
    return minimum ? item->minimumSize()
                   : item->sizeHint()
                         .expandedTo(item->minimumSize())
                         .boundedTo(item->maximumSize());
}

static QSize add_margins(QSize size, const QMargins& margins)
{
    return size + QSize(margins.left() + margins.right(),
                        margins.top() + margins.bottom());
}

QSize CustomMonomerDialog::dialogSizeHint(const bool minimum,
                                          const bool footer_below_view) const
{
    const auto& sketcher_ui = ui->sketcher_widget->m_ui;
    auto* view_layout = sketcher_ui->verticalLayout;
    QSize column_size(0, 0);
    int row_count = 0;
    for (int i = 0; i < view_layout->count(); ++i) {
        const auto* item = view_layout->itemAt(i);
        if (item->widget() == ui->button_bar || item->isEmpty()) {
            continue;
        }
        const auto item_size = layout_item_size(item, minimum);
        column_size.setWidth(std::max(column_size.width(), item_size.width()));
        column_size.rheight() += item_size.height();
        ++row_count;
    }
    QWidgetItem footer_item(ui->button_bar);
    const auto footer_size = layout_item_size(&footer_item, minimum);
    if (footer_below_view) {
        column_size.setWidth(
            std::max(column_size.width(), footer_size.width()));
        column_size.rheight() += footer_size.height();
        ++row_count;
    }
    column_size.rheight() +=
        std::max(0, row_count - 1) * view_layout->spacing();
    column_size = add_margins(column_size, view_layout->contentsMargins());

    const auto sidebar_size =
        layout_item_size(sketcher_ui->horizontalLayout->itemAt(0), minimum);
    QSize size(sidebar_size.width() + column_size.width() +
                   sketcher_ui->horizontalLayout->spacing(),
               std::max(sidebar_size.height(), column_size.height()));
    size = add_margins(size, sketcher_ui->horizontalLayout->contentsMargins() +
                                 ui->sketcher_widget->contentsMargins() +
                                 ui->verticalLayout_2->contentsMargins() +
                                 ui->sketcher_widget_holder->contentsMargins());
    if (!footer_below_view) {
        size.setWidth(std::max(size.width(), footer_size.width()));
        size.rheight() += footer_size.height() + ui->verticalLayout->spacing();
    }
    size = add_margins(
        size, ui->verticalLayout->contentsMargins() +
                  ui->verticalLayout->parentWidget()->contentsMargins() +
                  m_dlg_layout->contentsMargins());
    if (m_title_bar != nullptr) {
        QWidgetItem title_item(m_title_bar);
        const auto title_size = layout_item_size(&title_item, minimum);
        size.setWidth(std::max(size.width(), title_size.width()));
        size.rheight() += title_size.height() + layout()->spacing();
    }
    return add_margins(size, layout()->contentsMargins() + contentsMargins());
}

QSize CustomMonomerDialog::minimumSizeHint() const
{
    return m_layout_ready ? dialogSizeHint(true, true)
                          : ModalDialog::minimumSizeHint();
}

QSize CustomMonomerDialog::sizeHint() const
{
    return m_layout_ready
               ? dialogSizeHint(false, false).expandedTo(minimumSizeHint())
               : ModalDialog::sizeHint();
}

void CustomMonomerDialog::updateButtonBarPlacement()
{
    if (!m_layout_ready || m_updating_layout) {
        return;
    }
    QScopedValueRollback<bool> updating(m_updating_layout, true);
    const auto compact_minimum = minimumSizeHint();
    if (minimumSize() != compact_minimum) {
        setMinimumSize(compact_minimum);
    }
    const bool footer_below_view =
        height() < dialogSizeHint(true, false).height();
    if (footer_below_view == m_footer_below_view) {
        return;
    }
    if (footer_below_view) {
        ui->verticalLayout->removeWidget(ui->button_bar);
        ui->sketcher_widget->addWidgetBelowView(ui->button_bar);
    } else {
        ui->sketcher_widget->m_ui->verticalLayout->removeWidget(ui->button_bar);
        ui->verticalLayout->addWidget(ui->button_bar);
    }
    m_footer_below_view = footer_below_view;
    ui->button_bar->show();
    updateGeometry();
}

void CustomMonomerDialog::resizeEvent(QResizeEvent* event)
{
    updateButtonBarPlacement();
    ResizableModelDialog::resizeEvent(event);
}

bool CustomMonomerDialog::event(QEvent* event)
{
    const bool handled = ResizableModelDialog::event(event);
    if (event->type() == QEvent::LayoutRequest) {
        updateButtonBarPlacement();
    }
    return handled;
}

void CustomMonomerDialog::setRequiredAttachmentPoints(
    std::vector<int> required_attachment_points)
{
    std::ranges::sort(required_attachment_points);
    m_required_attachment_points.clear();
    std::ranges::unique_copy(required_attachment_points,
                             std::back_inserter(m_required_attachment_points));
}

void CustomMonomerDialog::addSMILES(const std::string& smiles)
{
    ui->sketcher_widget->addFromString(smiles, Format::EXTENDED_SMILES);
}

void CustomMonomerDialog::updateOkButton()
{
    bool valid = false;
    try {
        if (!ui->sketcher_widget->isEmpty()) {
            auto mol = ui->sketcher_widget->getRDKitMolecule();
            valid = mol->getNumAtoms() > 0 &&
                    RDKit::MolOps::getMolFrags(*mol, false).size() == 1;
            if (valid) {
                ui->sketcher_widget->getString(Format::EXTENDED_SMILES);
            }
        }
    } catch (const std::exception&) {
        valid = false;
    }
    ui->button_box->button(QDialogButtonBox::Ok)->setEnabled(valid);
}

/**
 * @return the test of the warning to use when the user has deleted the
 * specified required attachment points (i.e. attachment points that have a
 * bound connection)
 */
static QString
get_warning_text(const std::vector<int>& missing_attachment_points)
{
    const bool plural = missing_attachment_points.size() != 1;
    const auto attachment_points =
        format_r_group_list(missing_attachment_points);
    return attachment_points + (plural ? " have" : " has") +
           " been removed from this monomer but " + (plural ? "are" : "is") +
           " currently bound. Continuing will remove " +
           (plural ? "these connections." : "this connection.");
}

void CustomMonomerDialog::accept()
{
    const auto mol = ui->sketcher_widget->getRDKitMolecule();
    const auto duplicate_r_groups = get_duplicate_r_groups(*mol);
    if (!duplicate_r_groups.empty()) {
        show_error_dialog("Invalid Attachment Points",
                          "Multiple " +
                              format_r_group_list(duplicate_r_groups) +
                              " attachment points found. All attachment points "
                              "must be unique.",
                          this);
        return;
    }

    // warn before accepting if the user has deleted any required attachment
    // points (i.e. attachment points that have a bound connection)
    const auto missing_attachment_points =
        get_missing_required_attachment_points(*mol,
                                               m_required_attachment_points);
    const auto smiles = ui->sketcher_widget->getString(Format::EXTENDED_SMILES);
    if (!missing_attachment_points.empty()) {
        auto warning_text = get_warning_text(missing_attachment_points);
        auto* warning_dialog = show_warning_dialog("Remove Bound Connections?",
                                                   warning_text, this);
        connect(warning_dialog, &MessageBoxDialog::accepted, this,
                [this, smiles]() {
                    emit customMonomerAccepted(smiles, m_chain_type);
                    ResizableModelDialog::accept();
                });
        return;
    }

    emit customMonomerAccepted(smiles, m_chain_type);
    ResizableModelDialog::accept();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/custom_monomer_dialog.moc"
