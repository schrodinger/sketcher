#include "schrodinger/sketcher/molviewer/atom_item.h"

#include <GraphMol/ROMol.h>

#include <Qt>
#include <QPainter>
#include <QPointF>
#include <QString>

#include "schrodinger/sketcher/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

AtomItem::AtomItem(RDKit::Atom* atom, Fonts& fonts, AtomItemSettings& settings,
                   QGraphicsItem* parent) :
    AbstractGraphicsItem(parent),
    m_atom(atom),
    m_fonts(fonts),
    m_settings(settings)
{
    setZValue(static_cast<qreal>(ZOrder::ATOM));
    m_selection_highlighting_path.addEllipse(
        QPointF(0, 0), ATOM_SELECTION_HIGHLIGHTING_RADIUS,
        ATOM_SELECTION_HIGHLIGHTING_RADIUS);
    m_predictive_highlighting_path.addEllipse(
        QPointF(0, 0), ATOM_PREDICTIVE_HIGHLIGHTING_RADIUS,
        ATOM_PREDICTIVE_HIGHLIGHTING_RADIUS);
    updateCachedData();
}

int AtomItem::type() const
{
    return Type;
}

void AtomItem::updateCachedData()
{
    // prepareGeometryChange notifies the scene to schedule a repaint and to
    // schedule a recheck of our bounding rect.  It must be called *before* the
    // bounding rect changes.
    prepareGeometryChange();
    m_label_is_visible = determineLabelIsVisible();
    m_subrects.clear();
    if (m_label_is_visible) {
        m_pen.setColor(m_settings.getAtomColor(m_atom->getAtomicNum()));
        m_main_label_text = QString::fromStdString(m_atom->getSymbol());
        qreal label_width =
            m_fonts.m_main_label_fm.boundingRect(m_main_label_text).width();
        qreal label_height = m_fonts.m_main_label_fm.height() * 0.8;
        // center a label_width by label_height rectangle about (0, 0)
        m_main_label_rect = QRectF(-label_width * 0.5, -label_height * 0.5,
                                   label_width, label_height);
        m_subrects.push_back(m_main_label_rect);
    } else {
        m_main_label_text = QString();
        m_main_label_rect = QRectF();
    }

    // merge all of the subrects with the predictive highlighting path to create
    // the shape and bounding rect
    m_shape = QPainterPath(m_predictive_highlighting_path);
    for (QRectF rect : m_subrects) {
        QPainterPath rect_path;
        rect_path.addRect(rect);
        m_shape |= rect_path;
    }
    m_bounding_rect = m_shape.boundingRect();
}

bool AtomItem::determineLabelIsVisible() const
{
    if (m_settings.m_carbon_labels == CarbonLabels::ALL) {
        return true;
    }
    if (m_atom->getAtomicNum() != static_cast<unsigned int>(Element::C)) {
        return true;
    }
    // TODO: visible if isotope
    // TODO: visible if valence error
    // if (m_settings.m_show_valence_errors && ...)
    if (m_atom->getFormalCharge() != 0) {
        return true;
    }
    int num_bonds = m_atom->getDegree();
    if (num_bonds == 0) {
        return true;
    } else if (num_bonds == 1 &&
               m_settings.m_carbon_labels == CarbonLabels::TERMINAL) {
        return true;
    }
    return false;
}

const std::vector<QRectF>& AtomItem::getSubrects() const
{
    return m_subrects;
}

bool AtomItem::labelIsVisible() const
{
    return m_label_is_visible;
}

void AtomItem::paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
                     QWidget* widget)
{
    if (m_label_is_visible) {
        painter->save();
        painter->setFont(m_fonts.m_main_label_font);
        painter->setPen(m_pen);
        painter->drawText(m_main_label_rect, Qt::AlignCenter,
                          m_main_label_text);
        painter->restore();
    }
}

} // namespace sketcher
} // namespace schrodinger
