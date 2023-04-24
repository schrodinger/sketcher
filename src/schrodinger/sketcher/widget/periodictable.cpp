#include "schrodinger/sketcher/widget/periodictable.h"

#include <QHBoxLayout>
#include <QList>
#include <QPainter>
#include <QPushButton>
#include <QStyleOption>
#include <cctype>

#include "schrodinger/sketcher/ChemicalKnowledge.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_periodictable.h"

namespace schrodinger
{
namespace sketcher
{

PeriodicTableWidget::PeriodicTableWidget(QWidget* parent) : SketcherView(parent)
{
    ui.reset(new Ui::PeriodicTableForm());
    ui->setupUi(this);
    setWindowFlags(this->windowFlags() | Qt::Window);
    setWindowTitle("Periodic Table");
    setStyleSheet(PERIODIC_TABLE_STYLE);
    resize(395, 210);

    for (auto button : ui->group->buttons()) {
        auto symbol = button->text().toStdString();
        auto atomic_number = symbol_to_atomic_number(symbol);
        auto full_name = atomic_number_to_name(atomic_number);
        auto tooltip = std::to_string(atomic_number) + " " + full_name;
        ui->group->setId(button, atomic_number);
        button->setToolTip(QString::fromStdString(tooltip));
    }

    connect(ui->group, &QButtonGroup::idClicked, this,
            &PeriodicTableWidget::onButtonClicked);
}

PeriodicTableWidget::~PeriodicTableWidget() = default;

void PeriodicTableWidget::setCloseOnClick(bool close_on_click)
{
    m_close_on_click = close_on_click;
}

void PeriodicTableWidget::onButtonClicked(int button_id)
{
    auto element = Element(button_id);
    emit elementSelected(element);
    auto model = getModel();
    if (model != nullptr) {
        if (!model->hasActiveSelection()) {
            // update the model
            std::unordered_map<ModelKey, QVariant> kv_pairs = {
                {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::ATOM)},
                {ModelKey::ATOM_TOOL, QVariant::fromValue(AtomTool::ELEMENT)},
                {ModelKey::ELEMENT, QVariant::fromValue(element)},
            };
            model->setValues(kv_pairs);
        } else {
            // do not update the model
            model->pingValue(ModelKey::ELEMENT, element);
        }
    }
    if (m_close_on_click) {
        close();
    }
}

void PeriodicTableWidget::paintEvent(QPaintEvent*)
{
    // NOTE: Duplicated code; see ModularPopup and FileSaveImagePopup
    QStyleOption opt;
    opt.initFrom(this);
    QPainter p(this);
    style()->drawPrimitive(QStyle::PE_Widget, &opt, &p, this);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/periodictable.moc"
