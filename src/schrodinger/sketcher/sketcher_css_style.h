#pragma once

#include <unordered_map>

#include <QString>

namespace schrodinger
{
namespace sketcher
{

const QString SKETCHER_WIDGET_STYLE{
    "schrodinger--sketcher--SketcherWidget { background-color: white; }"
    "schrodinger--sketcher--View { background: transparent; }"
    "schrodinger--sketcher--ModularPopup { background-color: white; }"
    "QToolButton { border: none; }"
    "QToolButton:disabled { color: #E4E4E4; }"
    "QToolButton:checked { background-color: #d4e6f1; }"
    "QToolButton:hover:!checked { background-color: #edf7fc; }"};

const QString SELECTION_ACTIVE_STYLE{
    "QWidget { background-color: #f3f6f0; }"
    "QToolButton:checked { background-color: #d4e6f1; }"
    "QToolButton:hover:!checked { background-color: #e3ebda; }"};

const QString TOOL_BUTTON_CORNER_ARROW_STYLE{
    "QToolButton {"
    " border : none;"
    " padding: 0px;"
    " margin: 0px;"
    " spacing: 0px; }"
    " QToolButton::menu-indicator {"
    " image: url(\":icons/menu_corner_arrow.svg\");"
    " subcontrol-position: right bottom;"
    " subcontrol-origin: padding;"
    " left: -2px; }"
    "QToolButton::menu-indicator:disabled { color: #E4E4E4; "
    " image: url(\":icons/menu_corner_arrow_dis.svg\"); }"};

const QString PALETTE_TITLE_STYLE{"QLabel { font: bold 9px; color: #666666; }"
                                  "QLabel:disabled { color: #E4E4E4; }"};

const QString TEXT_LINK_STYLE{
    "QToolButton { font: bold 10px; color: #3d5d71; }"
    "QToolButton:disabled { color: #E4E4E4; }"
    "QToolButton:hover { color: #5b8aa8; background-color: transparent; }"};

/**
 * A text link that uses a larger font (12 px instead of 10 px) and a brighter
 * blue
 */
const QString BRIGHTER_TEXT_LINK_STYLE{
    "QToolButton { font: bold 12px; color: #009ec3; }"
    "QToolButton:disabled { color: #E4E4E4; }"
    "QToolButton:hover { color: #00b6e0; background-color: transparent; }"};

const QString ATOM_ELEMENT_OR_MONOMER_STYLE{
    "QToolButton { font: bold 14px; color: #333333; }"
    "QToolButton:disabled { color: #E4E4E4; }"};

/// style for an unknown monomer (e.g. amino acid X or nucleotide N)
const QString UNKNOWN_MONOMER_STYLE{
    "QToolButton { font: bold italic 14px; color: #333333; }"
    "QToolButton:disabled { color: #E4E4E4; }"};

const QString ATOM_QUERY_STYLE{
    "QLabel { font: italic 7px; color: #666666; }"
    "QToolButton { font: bold italic 14px; color: #606060; }"
    "QToolButton:disabled { color: #E4E4E4; }"};

const QString PERIODIC_TABLE_STYLE{
    "#PeriodicTableForm { background-color:white }"
    "QToolButton { font: 10px; color: black; border-style: none; "
    "    min-width: 21px; max-width:21px; min-height:21px; max-height:21px; }"
    "QToolButton[class='hydrogen'] { background-color: #b2bcc2; }"
    "QToolButton[class='alkali_metals'] { background-color: #b7d9ec; }"
    "QToolButton[class='alkaline_earth_metals'] { background-color: #8fbed9; }"
    "QToolButton[class='transition_metals'] { background-color: #f2d2c6; }"
    "QToolButton[class='other_metals'] { background-color: #f2d2c6; }"
    "QToolButton[class='metalloids'] { background-color: #e1baad; }"
    "QToolButton[class='non_metals'] { background-color: #f2e8b7; }"
    "QToolButton[class='halogens'] { background-color: #f2e392; }"
    "QToolButton[class='noble_gases'] { background-color: #eccc75; }"
    "QToolButton[class='lanthanides'] { background-color: #cce5c3; }"
    "QToolButton[class='actinides'] { background-color: #afd1a2; }"
    "QLabel { font: 10px; color: black; border-style: none; padding: 0; }"};

const QString RENDER_OPTIONS_STYLE{
    "QLabel { font: italic 11px; color: #666666; }"};

const QString RENDER_OPTIONS_POPUP_STYLE{
    "#FileSaveImagePopup { background-color:white }"};

const QString BOND_QUERY_STYLE{
    "QToolButton { font: bold 12px; color: #333333; }"};

const QString ENUMERATION_STYLE{
    "QToolButton { font: bold 14px; color: #444444; }"
    "QToolButton:disabled { color: #E4E4E4; }"};

/// style for the popup that allows the user to specify a custom nucleotide
const QString CUSTOM_NUCLEOTIDE_STYLE{
    "#CustomNucleotidePopup { background-color:white }"
    "QLabel { font: 9px; color: #666666; }"};

} // namespace sketcher
} // namespace schrodinger
