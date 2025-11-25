#pragma once

#include <vector>

#include <QPen>
#include <QRectF>
#include <QString>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/abstract_atom_or_monomer_item.h"
#include "schrodinger/sketcher/molviewer/atom_display_settings.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/fonts.h"

class QGraphicsItem;

namespace RDKit
{
class Atom;
} // namespace RDKit

namespace RDGeom
{
class Point3D;
} // namespace RDGeom

namespace schrodinger
{
namespace sketcher
{

enum class HsDirection { RIGHT, LEFT, UP, DOWN };

/**
 * A Qt graphics item for representing atoms (but not monomers) in a molviewer
 * Scene.  This class is responsible for painting all aspects of the atom label,
 * which can include
 *
 * - the main text label for the atom.  This normally contains the atom's
 *   element, but can also contain the R group or a query string.
 * - implicit hydrogens
 * - lone pair indicator
 * - valence error indicator
 *
 * as well as a handful of additional aspects, depending on the atom.  Note that
 * labels are frequently not shown for certain carbons, so this class won't
 * paint anything in those scenarios.  Also note that this class is *not*
 * responsible for drawing bonds.  See BondItem for that.
 */
class SKETCHER_API AtomItem : public AbstractAtomOrMonomerItem
{
  public:
    /**
     * Note that this class does not take ownership of atom, fonts, or settings.
     * These objects must not be destroyed while the AtomItem is in use.
     *
     * @param atom The RDKit atom that this item should represent.
     *
     * @param fonts The fonts to be used when painting this item.  To change the
     * fonts, the calling code (i.e. the Scene) must update this object and then
     * call updateCachedData.
     *
     * @param settings The settings to be used for painting this item.  To
     * change settings, the calling code (i.e. the Scene) must update this
     * object and then call updateCachedData.
     *
     * @param parent The Qt parent for this item.  See the QGraphicsItem
     * documentation for additional information.
     *
     * @pre atom != nullptr
     * @pre atom->hasOwningMol()
     */
    AtomItem(const RDKit::Atom* atom, const Fonts& fonts,
             const AtomDisplaySettings& settings,
             QGraphicsItem* parent = nullptr);

    // Type and type() are required for qgraphicsitem_cast.  Note that this enum
    // implementation is based on the sample code from the QGraphicsItem::Type
    // documentation.
    enum { Type = static_cast<int>(ItemType::ATOM) };
    int type() const override;

    // Overridden AbstractGraphicsItem method
    void updateCachedData() override;

    // Overridden QGraphicsItem methods
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget = nullptr) override;

    /**
     * Return a list of bounding rectangles for all individual components of
     * this item.  These subrects are used to precisely trim the bond lines so
     * that they don't overlap with the atom labels.  Note that this method
     * returns subrects starting from the center of the atom and working
     * outwards.  This behavior is required for proper functioning of bond line
     * trimming due to assumptions made in trim_line_to_rect.
     */
    const std::vector<QRectF> getSubrects() const;

    /**
     * @param avoid_subrects whether or not to consider the atom's subrects
     * in the calculation: when picking a position for something close to the
     * atom (e.g. chiral labels), this make sure that the returned position
     * doesn't clash with the atom labels. When picking positions for something
     * far away from the label (e.g. another atom bound to this), this should be
     * false
     * @return a position on the canvas (in coordinates relative to this item)
     * that is near this atom but doesn't clash with it or its neighbors.
     */
    QPointF findPositionInEmptySpace(bool avoid_subrects) const;

    /**
     * Return whether the label for this atom is visible.  (Labels for some
     * carbons may be hidden depending on the atom item settings.)
     */
    bool labelIsVisible() const;

    QRectF getChiralityLabelRect() const;
    QString getChiralityLabelText() const;

  protected:
    /**
     * @return the direction Hs labels would be drawn on this atom if the label
     * is shown and it has hydrogens
     */
    HsDirection findHsDirection() const;

    QString m_main_label_text;
    QRectF m_main_label_rect;
    QRectF m_H_label_rect;
    QRectF m_H_count_label_rect;
    QString m_H_count_label_text;
    QRectF m_charge_and_radical_label_rect;
    QString m_charge_and_radical_label_text;
    QString m_isotope_label_text;
    QString m_chirality_label_text;
    QRectF m_chirality_label_rect;
    QString m_mapping_label_text;
    QRectF m_mapping_label_rect;
    QRectF m_isotope_label_rect;
    QRectF m_query_label_rect;
    QString m_query_label_text;
    QPainterPath m_squiggle_path; // used for attachment point squiggle

    const Fonts& m_fonts;
    const AtomDisplaySettings& m_settings;
    bool m_label_is_visible;
    bool m_valence_error_is_visible;
    QPen m_pen;
    QPen m_valence_error_pen;
    QPen m_chirality_pen;
    QPen m_squiggle_pen;
    QPen m_query_label_line_pen;
    QBrush m_valence_error_brush;

    /**
     * @return A painter path with the specified radius for use with either
     * selection highlighting or predictive highlighting.
     *
     */
    QPainterPath calcHighlightingPath(qreal radius);

    /**
     * @return a list of the rectangles that make up the atom label. This is
     * used to calculate subrects with the visible rectangles
     */
    std::vector<QRectF> getLabelRects() const;

    /**
     * Determine what type of label we should display for this atom
     * @return A tuple of
     *   - the text to display in the main label
     *   - the path for painting the attachment point squiggle.  This path will
     *     be empty if this atom is not an attachment point.
     *   - whether the main label should be visible
     *   - whether the valence error should be visible
     *   - whether this atom needs additional labels (e.g. charge label, isotope
     *     label, Hs label, etc.)
     *   - whether Hs labels should be displayed when the previous parameter is
     *     true (e.g. queries show additional labels but not Hs)
     *   - the text to display for queries. This is displayed on a separate
     *     label next to the atom. When this text is shown, the main label is
     *     hidden, unless the atom has no bonds, in which case we show an "*"
     */
    std::tuple<QString, QPainterPath, bool, bool, bool, bool, QString>
    determineLabelType() const;

    /**
     * @return the text to display in query label for this atom
     */
    QString getQueryLabel() const;

    /**
     * @return the smarts representation of the advanced properties of this atom
     */
    QString advancedPropertiesSmarts() const;

    /**
     * @return A path for painting the attachment point squiggle.  This method
     * should only be called for AtomItems that represent an attachment point.
     */
    QPainterPath getWavyLine() const;

    /**
     * @return whether the label for this atom is visible.  (Labels for some
     * carbons may be hidden depending on the atom item settings.)  Note that
     * this method is only called for "normal" atoms (i.e. atoms that represent
     * an element, as opposed to query atoms, R-groups, etc, all of which are
     * always visible) and that it does *not* take the valence error into
     * account.  If a valence error is visible, the label should always be
     * visible.
     */
    bool determineLabelIsVisible() const;

    /**
     * @Return whether the associated atom has a valence error that should be
     * displayed
     */
    bool determineValenceErrorIsVisible() const;

    /**
     * set all label rects to invalid rects and all label text to empty strings
     */
    void clearLabels();

    /**
     * create a label for the isotope number (if present)
     */
    void updateIsotopeLabel();

    /**
     * create a label for the mapping number (if present)
     */
    void updateMappingLabel();

    /**
     * Get tooltip text based on atom labels (chirality or query).
     * Returns tooltip text for labels related to stereochemistry or query
     * features. If an atom has both stereochemistry and query features, both
     * are shown on separate lines (e.g., "Stereo: (R)\nQuery: X2").
     * @return tooltip text, or empty string if no tooltip should be shown
     */
    QString getTooltip() const;

    /**
     * create a label for the charge and radical information  (if present)
     */
    void updateChargeAndRadicalLabel();

    /**
     * create a label for the hydrogens on this atom (if any)
     */
    void updateHsLabel();

    /**
     * position all the labels around the atom. Hs can be to the right, left,
     * top or bottom depending on the binding pattern
     */
    void positionLabels();

    /**
     * create and position labels (dots) for the unpaired electrons on this
     * atom (if any)
     */
    void updateUnpairedElectronsLabels();

    /**
     * create and position labels for stereo annotations (if they are present)
     * on this atom
     */
    void updateChiralityLabel();
};

} // namespace sketcher
} // namespace schrodinger
