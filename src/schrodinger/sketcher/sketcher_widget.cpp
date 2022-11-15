#include "schrodinger/sketcher/sketcher_widget.h"

#include <iostream>

#include <QComboBox>
#include <QDoubleSpinBox>
#include <QLineEdit>
#include <QPushButton>
#include <QVBoxLayout>
#include <QWidget>

#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/view.h"
#include "schrodinger/sketcher/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_sketcher_widget.h"

namespace schrodinger
{
namespace sketcher
{

SketcherWidget::SketcherWidget(QWidget* parent) : QWidget(parent)
{
    m_sketcher_model = new SketcherModel(this);

    // FIXME: ignore UI and use temporary setup for testing (below)
    // ui.reset(new Ui::SketcherWidgetForm());
    // ui->setupUi(this);
    // ui->top_bar_wdg->setModel(m_sketcher_model);
    // ui->side_bar_wdg->setModel(m_sketcher_model);

    const QString DEFAULT_SMILES{"c1nccc2n1ccc2"};

    auto layout = new QVBoxLayout(this);

    auto scene = new Scene();
    auto view = new View(scene);
    layout->addWidget(view);

    auto smiles_le = new QLineEdit(this);
    smiles_le->setText(DEFAULT_SMILES);
    auto smile_btn = new QPushButton("Load SMILES", this);

    auto load_smiles = [scene, smiles_le]() {
        auto smiles = smiles_le->text();
        try {
            scene->loadSmiles(smiles.toStdString());
        } catch (const std::exception& exc) {
            std::cerr << exc.what() << std::endl;
        }
    };

    connect(smiles_le, &QLineEdit::returnPressed, load_smiles);
    connect(smile_btn, &QPushButton::clicked, load_smiles);
    layout->addWidget(smiles_le);
    layout->addWidget(smile_btn);

    auto font_sb = new QDoubleSpinBox(this);
    font_sb->setValue(scene->fontSize());
    connect(font_sb, &QDoubleSpinBox::valueChanged, scene, &Scene::setFontSize);
    layout->addWidget(font_sb);

    auto carbon_labels_combo = new QComboBox(this);
    carbon_labels_combo->addItems(
        {"No carbon labels", "Terminal carbons only", "All atoms labeled"});
    connect(carbon_labels_combo, &QComboBox::currentIndexChanged, scene,
            [scene](auto i) { scene->setCarbonsLabeled(CarbonLabels(i)); });
    layout->addWidget(carbon_labels_combo);

    scene->loadSmiles(DEFAULT_SMILES.toStdString());
}

SketcherWidget::~SketcherWidget() = default;

} // namespace sketcher
} // namespace schrodinger
