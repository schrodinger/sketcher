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

const QString SMILES_PLACEHOLDER_TEXT = "***";
const qsizetype MAX_MONOMER_LABEL_LENGTH = 6;

const QColor DEFAULT_AA_BACKGROUND_COLOR =
    MONOMER_COLOR_MAP.at(MonomerColorType::OTHER);

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

} // namespace sketcher
} // namespace schrodinger
