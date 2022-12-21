#pragma once
#include "schrodinger/sketcher/definitions.h"
#include <QToolButton>

namespace schrodinger
{
namespace sketcher
{

enum class Element;

class SKETCHER_API ModularElementButton : public QToolButton
{
    Q_OBJECT
  public:
    ModularElementButton(QWidget* parent = nullptr);

    /**
     * @return the element currently assigned to this button
     */
    Element getElement() const;

    /**
     * @param element The element to assign to this button
     */
    void setElement(Element element);

  protected:
    void updateButtonText();
    Element m_element;
};

} // namespace sketcher
} // namespace schrodinger
