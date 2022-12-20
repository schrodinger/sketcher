#include "schrodinger/sketcher/widget/modular_element_button.h"
#include "schrodinger/sketcher/ChemicalKnowledge.h"
#include "schrodinger/sketcher/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

ModularElementButton::ModularElementButton(QWidget* parent) :
    QToolButton(parent)
{
}

Element ModularElementButton::getElement() const
{
    return m_element;
}

void ModularElementButton::setElement(Element element)
{
    m_element = element;
    updateButtonText();
}

void ModularElementButton::updateButtonText()
{
    setText(QString::fromStdString(
        atomic_number_to_symbol(static_cast<int>(m_element))));
}

} // namespace sketcher
} // namespace schrodinger

#include "modular_element_button.moc"
