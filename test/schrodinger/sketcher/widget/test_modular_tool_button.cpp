#define BOOST_TEST_MODULE Test_Sketcher
#include <iostream>

#include <QButtonGroup>
#include <QLayout>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/widget/modular_popup.h"
#include "schrodinger/sketcher/widget/modular_tool_button.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

class DummyPopup : public ModularPopup
{

  public:
    DummyPopup() : ModularPopup()
    {
        setModel(new SketcherModel(this));
        auto button_layout = new QHBoxLayout(this);
        auto group = new QButtonGroup(this);
        for (int idx = 0; idx < m_button_count; ++idx) {
            auto button = new QToolButton(this);
            button->setText(QString::number(idx));
            group->addButton(button, idx);
            button_layout->addWidget(button);
        }
        setButtonGroup(group);
    }

    int m_button_count = 3;

  protected:
    void generateButtonPackets() override
    {
        auto button_layout = layout();
        for (int idx = 0; idx < button_layout->count(); ++idx) {
            auto button =
                dynamic_cast<QToolButton*>(button_layout->itemAt(idx));
            if (button != nullptr) {
                int button_idx = button->text().toInt();
                m_button_packets.emplace_back(button, button_idx);
            }
        }
    }

    int getButtonIDToCheck() override
    {
        auto model = getModel();
        if (model == nullptr) {
            return -1;
        }
        int rgroup_number = model->getValueInt(ModelKey::RGROUP_NUMBER);
        return rgroup_number >= m_button_count ? -1 : rgroup_number;
    }
};

class TestModularToolButton : public ModularToolButton
{
};

/**
 * Verify that the text on the button matches the text of the selected button
 * in the popoup.
 */
BOOST_AUTO_TEST_CASE(update_text)
{

    DummyPopup popup;
    TestModularToolButton button;
    button.setPopupWidget(&popup);

    for (int idx = 0; idx < popup.m_button_count; ++idx) {
        button.setEnumItem(idx);
        BOOST_TEST(button.getEnumItem() == idx);
        bool conversion_worked = false;
        int text_as_int = button.text().toInt(&conversion_worked);
        BOOST_TEST(conversion_worked);
        BOOST_TEST(text_as_int == idx);
        BOOST_TEST((button.toolTip()).contains("press & hold") == true);
    }
}

} // namespace sketcher
} // namespace schrodinger
