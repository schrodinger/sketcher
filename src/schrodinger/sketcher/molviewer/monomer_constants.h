/**
 * Constants use in the molviewer Qt graphics view classes for monomeric models
 */

#pragma once

#include <unordered_map>

#include <QColor>
#include <QString>

namespace schrodinger
{
namespace sketcher
{

enum class MonomerColorType {
    // Amino Acid Monomers
    ACIDIC = 1,
    BASIC,
    ALIPHATIC,
    AROMATIC,
    HYDROPHILIC,
    THIOL,
    IMINO_ACID,
    VARIABLE_POLAR,

    // Nucleic Acid Monomers
    NA_BACKBONE,
    ADENINE,
    CYTOSINE,
    GUANINE,
    THYMINE,
    URACIL,

    // Unnatural AA or Nucleobase
    OTHER,

    // Chem -- everything else
    CHEM,

    // Text and Outlines
    BLACK,
    GRAY3,
    GRAY4,
    GRAY6,
    DARK_RED,

    // Default Background
    WHITE
};

const std::unordered_map<MonomerColorType, QColor> MONOMER_COLOR_MAP{
    // Amino Acid Monomers
    {MonomerColorType::ACIDIC, QColor("#FCB7D6")},
    {MonomerColorType::BASIC, QColor("#B8D1FC")},
    {MonomerColorType::ALIPHATIC, QColor("#B4FCB4")},
    {MonomerColorType::AROMATIC, QColor("#FBC6B7")},
    {MonomerColorType::HYDROPHILIC, QColor("#B9FCFD")},
    {MonomerColorType::THIOL, QColor("#FCFCB6")},
    {MonomerColorType::IMINO_ACID, QColor("#DADADA")},
    {MonomerColorType::VARIABLE_POLAR, QColor("#F3E0F9")},

    // Nucleic Acid Monomers
    {MonomerColorType::NA_BACKBONE, QColor("#E0E0E0")},
    {MonomerColorType::ADENINE, QColor("#B0A8F3")},
    {MonomerColorType::CYTOSINE, QColor("#FF9CB6")},
    {MonomerColorType::GUANINE, QColor("#6CF2BF")},
    {MonomerColorType::THYMINE, QColor("#F0F693")},
    {MonomerColorType::URACIL, QColor("#F1CC6C")},

    // Unnatural AA or Nucleobase
    {MonomerColorType::OTHER, QColor("#F0E7E1")},

    // Chem -- everything else
    {MonomerColorType::CHEM, QColor("#E5E4E2")},

    // Text and Outlines
    {MonomerColorType::BLACK, QColor("#000000")},
    {MonomerColorType::GRAY3, QColor("#333333")},
    {MonomerColorType::GRAY4, QColor("#444444")},
    {MonomerColorType::GRAY6, QColor("#666666")},
    {MonomerColorType::DARK_RED, QColor("#330000")},

    // Default Background
    {MonomerColorType::WHITE, QColor("#FFFFFF")},
};

const std::unordered_map<std::string, QColor> AMINO_ACID_COLOR_BY_RES_NAME{
    {"A", MONOMER_COLOR_MAP.at(MonomerColorType::ALIPHATIC)},
    {"B", MONOMER_COLOR_MAP.at(MonomerColorType::VARIABLE_POLAR)},
    {"C", MONOMER_COLOR_MAP.at(MonomerColorType::THIOL)},
    {"D", MONOMER_COLOR_MAP.at(MonomerColorType::ACIDIC)},
    {"E", MONOMER_COLOR_MAP.at(MonomerColorType::ACIDIC)},
    {"F", MONOMER_COLOR_MAP.at(MonomerColorType::AROMATIC)},
    {"G", MONOMER_COLOR_MAP.at(MonomerColorType::ALIPHATIC)},
    {"H", MONOMER_COLOR_MAP.at(MonomerColorType::BASIC)},
    {"I", MONOMER_COLOR_MAP.at(MonomerColorType::ALIPHATIC)},
    {"K", MONOMER_COLOR_MAP.at(MonomerColorType::BASIC)},
    {"L", MONOMER_COLOR_MAP.at(MonomerColorType::ALIPHATIC)},
    {"M", MONOMER_COLOR_MAP.at(MonomerColorType::ALIPHATIC)},
    {"N", MONOMER_COLOR_MAP.at(MonomerColorType::HYDROPHILIC)},
    {"P", MONOMER_COLOR_MAP.at(MonomerColorType::IMINO_ACID)},
    {"Q", MONOMER_COLOR_MAP.at(MonomerColorType::HYDROPHILIC)},
    {"R", MONOMER_COLOR_MAP.at(MonomerColorType::BASIC)},
    {"S", MONOMER_COLOR_MAP.at(MonomerColorType::HYDROPHILIC)},
    {"T", MONOMER_COLOR_MAP.at(MonomerColorType::HYDROPHILIC)},
    {"V", MONOMER_COLOR_MAP.at(MonomerColorType::ALIPHATIC)},
    {"W", MONOMER_COLOR_MAP.at(MonomerColorType::AROMATIC)},
    {"Y", MONOMER_COLOR_MAP.at(MonomerColorType::AROMATIC)},
    {"Z", MONOMER_COLOR_MAP.at(MonomerColorType::VARIABLE_POLAR)},
};

const std::unordered_map<std::string, QColor> NUCLEIC_ACID_COLOR_BY_RES_NAME{
    {"A", MONOMER_COLOR_MAP.at(MonomerColorType::ADENINE)},
    {"C", MONOMER_COLOR_MAP.at(MonomerColorType::CYTOSINE)},
    {"G", MONOMER_COLOR_MAP.at(MonomerColorType::GUANINE)},
    {"T", MONOMER_COLOR_MAP.at(MonomerColorType::THYMINE)},
    {"U", MONOMER_COLOR_MAP.at(MonomerColorType::URACIL)},
};

const std::string SMILES_PLACEHOLDER_TEXT = "***";
const qsizetype MAX_MONOMER_LABEL_LENGTH = 6;

const QColor DEFAULT_AA_BACKGROUND_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::OTHER);
const QColor MONOMER_LABEL_TEXT_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::GRAY3);
const qreal MONOMER_LABEL_TEXT_WIDTH = 0.5;

const qreal STANDARD_AA_BORDER_WIDTH = 25;
const qreal STANDARD_AA_BORDER_HEIGHT = 25;
const qreal STANDARD_AA_BORDER_LINE_WIDTH = 2;
const QColor STANDARD_AA_BORDER_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::GRAY4);

const qreal D_AA_BORDER_WIDTH = 32;
const qreal D_AA_BORDER_HEIGHT = 25;
const qreal D_AA_BORDER_LINE_WIDTH = 3;
// font size relative to the standard amino acid label font size
const qreal D_AA_FONT_RATIO = 1.0;
const QColor D_AA_BORDER_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::DARK_RED);

const qreal OTHER_AA_BORDER_WIDTH = 25;
const qreal OTHER_AA_BORDER_HEIGHT = 25;
const qreal OTHER_AA_BORDER_LINE_WIDTH = 3;
// font size relative to the standard amino acid label font size
const qreal OTHER_AA_FONT_RATIO = 0.875;
const QColor OTHER_AA_BORDER_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::DARK_RED);

const qreal AA_ROUNDING_RADIUS = 5;
const qreal MONOMER_RECT_MINIMUM_HORIZONTAL_TEXT_PADDING = 12.0;

const qreal NA_BASE_BORDER_WIDTH = 28.75;
const qreal NA_BASE_BORDER_HEIGHT = 28.75;
const QColor DEFAULT_NA_BACKGROUND_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::OTHER);
const qreal DIAMOND_AND_ELLIPSE_MIN_ASPECT_RATIO = 0.2;

const qreal STANDARD_NA_BORDER_LINE_WIDTH = 2;
const QColor STANDARD_NA_BORDER_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::GRAY4);

const qreal OTHER_NA_BORDER_LINE_WIDTH = 3;
const qreal OTHER_NA_SUGAR_AND_PHOSPHATE_FONT_RATIO = 1.0;
const qreal OTHER_NA_BASE_FONT_RATIO = 0.925;
const QColor OTHER_NA_BORDER_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::DARK_RED);

const qreal NA_SUGAR_BORDER_WIDTH = 25;
const qreal NA_SUGAR_BORDER_HEIGHT = 25;
const qreal NA_PHOSPHATE_BORDER_WIDTH = 26.25;
const qreal NA_PHOSPHATE_BORDER_HEIGHT = 26.25;
const QColor NA_BACKBONE_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::NA_BACKBONE);

const qreal CHEM_MONOMER_BORDER_WIDTH = 25;
const qreal CHEM_MONOMER_BORDER_HEIGHT = 25;
const qreal CHEM_MONOMER_BORDER_LINE_WIDTH = 2;
const QColor CHEM_MONOMER_BORDER_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::GRAY4);
const QColor CHEM_MONOMER_RECT_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::CHEM);

// the thickness of the "halo" around monomers for predictive and selection
// highlighting
const qreal MONOMER_PREDICTIVE_HIGHLIGHTING_THICKNESS = 7.5;
const qreal MONOMER_SELECTION_HIGHLIGHTING_THICKNESS = 8.5;

// The absolute minimum size for the colored area drawn behind a monomer. This
// floor is only relevant when running without a display, in which case this
// value can prevent, e.g., divide-by-zero errors.
const qreal MONOMER_BORDER_MIN_SIZE = 0.01;

const QColor AA_LINEAR_CONNECTOR_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::GRAY4);
const qreal AA_LINEAR_CONNECTOR_WIDTH = 5;

const QColor AA_BRANCHING_CONNECTOR_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::GRAY6);
const qreal AA_BRANCHING_CONNECTOR_WIDTH = 3;

const QColor DISULFIDE_CONNECTOR_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::GRAY6);
const qreal DISULFIDE_CONNECTOR_WIDTH = 3;

const QColor NA_BACKBONE_CONNECTOR_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::GRAY4);
const qreal NA_BACKBONE_CONNECTOR_WIDTH = 4;

const QColor NA_BASE_CONNECTOR_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::GRAY6);
const qreal NA_BASE_CONNECTOR_WIDTH = 4;

const QColor NA_BACKBONE_TO_BASE_CONNECTOR_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::GRAY6);
const qreal NA_BACKBONE_TO_BASE_CONNECTOR_WIDTH = 2;

const QColor CHEM_CONNECTOR_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::GRAY6);
const qreal CHEM_CONNECTOR_WIDTH = 2;

const qreal MONOMER_CONNECTOR_ARROWHEAD_RADIUS = 6;

/**
 * A scaling factor for generating the cursor hints for monomers. If this value
 * is greater than one, then the cursor hint for smaller monomers (as measured
 * by the size of the border outline around the label) will be zoomed out,
 * ensuring that the cursor hint is smaller than an actual monomer in the scene.
 * Increasing this value will zoom out further (i.e. will decrease the size of
 * the cursor hint). Note that this setting will have no effect on the cursor
 * hint of larger monomers, as those are already zoomed out to ensure that they
 * fit within the maximum cursor hint size (CURSOR_HINT_IMAGE_SIZE).
 */
const qreal MONOMER_CURSOR_HINT_MIN_SCENE_SIZE_SCALE = 1.4;

} // namespace sketcher
} // namespace schrodinger
