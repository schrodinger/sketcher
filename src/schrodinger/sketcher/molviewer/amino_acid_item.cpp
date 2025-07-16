#include "schrodinger/sketcher/molviewer/amino_acid_item.h"

#include <cctype>

#include <QPainter>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/MonomerInfo.h>

#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/rdkit_extensions/helm.h"

namespace schrodinger
{
namespace sketcher
{

enum class AminoAcidType { STANDARD, D, OTHER };

AminoAcidItem::AminoAcidItem(const RDKit::Atom* monomer, const Fonts& fonts,
                             QGraphicsItem* parent) :
    AbstractAtomOrMonomerItem(monomer, parent),
    m_fonts(fonts)
{
    setZValue(static_cast<qreal>(ZOrder::MONOMER));
    updateCachedData();
}

int AminoAcidItem::type() const
{
    return Type;
}

static QString elide_text(const QString& text)
{
    if (text.length() > MAX_MONOMER_LABEL_LENGTH) {
        // shorten the string and add ellipses
        return text.left(MAX_MONOMER_LABEL_LENGTH - 1) +
               QString::fromUtf8("\u2026");
    } else {
        return text;
    }
}

static QString get_label_text(const RDKit::Atom* const monomer,
                              const RDKit::AtomPDBResidueInfo* const res_info)
{
    auto is_smiles = monomer->getProp<bool>(SMILES_MONOMER);
    if (is_smiles) {
        return SMILES_PLACEHOLDER_TEXT;
    }
    auto res_name = res_info->getResidueName();
    return elide_text(QString::fromStdString(res_name));
}

static AminoAcidType
get_amino_acid_type(const RDKit::AtomPDBResidueInfo* const res_info)
{
    auto res_name = res_info->getResidueName();
    if (res_name.length() == 1) {
        return AminoAcidType::STANDARD;
    } else if (std::tolower(res_name[0]) == 'd') {
        return AminoAcidType::D;
    } else {
        return AminoAcidType::OTHER;
    }
}

void AminoAcidItem::updateCachedData()
{
    prepareGeometryChange();

    const auto* monomer_info = m_atom->getMonomerInfo();
    const auto* res_info =
        dynamic_cast<const RDKit::AtomPDBResidueInfo*>(monomer_info);
    if (res_info == nullptr) {
        // this monomer is something other than an amino acid
        throw std::runtime_error("Atom has no residue info");
    }
    auto res_name = res_info->getResidueName();
    m_main_label_text = get_label_text(m_atom, res_info);
    qreal border_width, border_height, border_line_width;
    QColor border_color;
    QFont label_font;
    switch (get_amino_acid_type(res_info)) {
        case AminoAcidType::STANDARD:
            border_width = STANDARD_AA_BORDER_WIDTH;
            border_height = STANDARD_AA_BORDER_HEIGHT;
            border_line_width = STANDARD_AA_BORDER_LINE_WIDTH;
            border_color = STANDARD_AA_BORDER_COLOR;
            m_main_label_font = m_fonts.m_main_label_font;
            break;
        case AminoAcidType::D:
            border_width = D_AA_BORDER_WIDTH;
            border_height = D_AA_BORDER_HEIGHT;
            border_line_width = D_AA_BORDER_LINE_WIDTH;
            border_color = D_AA_BORDER_COLOR;
            m_main_label_font = m_fonts.m_d_amino_acid_font;
            break;
        case AminoAcidType::OTHER:
            border_width = OTHER_AA_BORDER_WIDTH;
            border_height = OTHER_AA_BORDER_HEIGHT;
            border_line_width = OTHER_AA_BORDER_LINE_WIDTH;
            border_color = OTHER_AA_BORDER_COLOR;
            m_main_label_font = m_fonts.m_other_amino_acid_font;
            break;
    }
    m_border_pen.setWidthF(border_line_width);
    m_border_pen.setColor(border_color);
    m_border_rect.setWidth(border_width);
    m_border_rect.setHeight(border_height);
    // TODO: add selection, predictive highlighting path - set bounding rect to
    // be largest of those
    m_bounding_rect = m_border_rect;
    m_shape.clear();
    m_shape.addRect(m_bounding_rect);

    auto color_find = AMINO_ACID_COLOR_BY_RES_NAME.find(res_name);
    if (color_find == AMINO_ACID_COLOR_BY_RES_NAME.end()) {
        m_border_brush.setColor(DEFAULT_AA_BACKGROUND_COLOR);
    } else {
        m_border_brush.setColor(color_find->second);
    }
}

void AminoAcidItem::paint(QPainter* painter,
                          const QStyleOptionGraphicsItem* option,
                          QWidget* widget)
{
    painter->save();
    painter->setPen(m_border_pen);
    painter->setBrush(m_border_brush);
    painter->drawRoundedRect(m_border_rect, AA_ROUNDING_RADIUS,
                             AA_ROUNDING_RADIUS);
    painter->restore();
}

} // namespace sketcher
} // namespace schrodinger
