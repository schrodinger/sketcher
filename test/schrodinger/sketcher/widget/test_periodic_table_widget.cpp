
#define BOOST_TEST_MODULE Test_Sketcher

#include <QPushButton>
#include <QSignalSpy>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/ui/ui_periodic_table_widget.h"
#include "schrodinger/sketcher/widget/periodic_table_widget.h"

Q_DECLARE_METATYPE(schrodinger::sketcher::ModelKey);
Q_DECLARE_METATYPE(std::unordered_set<schrodinger::sketcher::ModelKey>);

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

class TestPeriodicTableWidget : public PeriodicTableWidget
{
  public:
    using PeriodicTableWidget::ui;
};

BOOST_AUTO_TEST_CASE(testGetElementButtons)
{
    TestPeriodicTableWidget pt_widget;
    std::unordered_set<std::string> symbols;
    for (auto button : pt_widget.ui->group->buttons()) {
        std::string text = button->text().toStdString();
        symbols.insert(text);
        auto exp_atomic_number = pt_widget.ui->group->id(button);
        BOOST_TEST(symbol_to_atomic_number(text) == exp_atomic_number);
    }

    // Make sure all buttons and symbols are unique
    int max_element = 118;
    BOOST_TEST(pt_widget.ui->group->buttons().size() == max_element);
    BOOST_TEST(symbols.size() == max_element);
    for (int i = 1; i <= max_element; ++i) {
        // and double check that all elements are mapped to buttons
        BOOST_TEST(pt_widget.ui->group->button(i) != nullptr);
    }
}

BOOST_AUTO_TEST_CASE(model_controller)
{
    auto test_scene = TestScene::getScene();
    auto model = test_scene->m_sketcher_model;
    TestPeriodicTableWidget pt_widget;
    pt_widget.setModel(model);
    import_mol_text(test_scene->m_mol_model, "CC");

    qRegisterMetaType<ModelKey>("ModelKey");
    qRegisterMetaType<std::unordered_set<ModelKey>>(
        "std::unordered_set<ModelKey>");
    QSignalSpy pinged_spy(model, &SketcherModel::valuePinged);
    QSignalSpy changed_spy(model, &SketcherModel::valuesChanged);

    for (auto button : pt_widget.ui->group->buttons()) {
        pinged_spy.clear();
        changed_spy.clear();
        button->click();
        auto exp_atomic_number = pt_widget.ui->group->id(button);
        BOOST_TEST(model->getValueInt(ModelKey::ELEMENT) == exp_atomic_number);
        BOOST_TEST(pinged_spy.count() == 3);
        BOOST_TEST(changed_spy.count() == 1);
    }

    // If an active selection, this widget should not change the model
    test_scene->m_mol_model->selectAll();
    pinged_spy.clear();
    changed_spy.clear();

    int exp_key_int = static_cast<int>(ModelKey::ELEMENT);
    auto exp_model_value = model->getElement();
    for (auto button : pt_widget.ui->group->buttons()) {
        button->click();
        BOOST_TEST(pinged_spy.count() == 1);
        BOOST_TEST(changed_spy.count() == 0);
        auto args = pinged_spy.takeLast();
        auto key_int = static_cast<int>(args.at(0).value<ModelKey>());
        auto value_int = args.at(1).value<QVariant>().toInt();
        BOOST_TEST(key_int == exp_key_int);
        auto exp_atomic_number = pt_widget.ui->group->id(button);
        BOOST_TEST(value_int == exp_atomic_number);
        BOOST_TEST(model->getElement() == exp_model_value);
    }
}

} // namespace sketcher
} // namespace schrodinger
