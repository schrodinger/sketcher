#include "schrodinger/sketcher/widget/widget_utils.h"

#include <QAbstractButton>
#include <QButtonGroup>

namespace schrodinger
{
namespace sketcher
{

/**
 * Uncheck all buttons in a button group.
 *
 * @param group The button group to uncheckify
 */
static void uncheck_all(QButtonGroup* group)
{
    bool exclusive = group->exclusive();
    group->setExclusive(false);
    for (auto button : group->buttons()) {
        button->setChecked(false);
    }
    group->setExclusive(exclusive);
}

void check_button_or_uncheck_group(QAbstractButton* button, QButtonGroup* group)
{
    if (button == nullptr) {
        uncheck_all(group);
    } else {
        button->setChecked(true);
    }
}

} // namespace sketcher
} // namespace schrodinger
