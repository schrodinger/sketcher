#pragma once

#include <unordered_map>

#include <QString>

namespace schrodinger
{
namespace sketcher
{

const QString GENERAL_STYLE{
    "QMainWindow { background-color: white; }"
    "schrodinger--sketcher--ModularPopup { background-color: white; }"};

const QString TOOL_BUTTON_STYLE{
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
    " background-image:url(':icons/menu_corner_arrow.svg'); "
    " background-position: bottom right; background-repeat: no-repeat; }"
    "QToolButton:disabled { color: #E4E4E4; "
    " background-image:url(':icons/menu_corner_arrow_dis.svg'); }"};

const QString PALETTE_TITLE_STYLE{"QLabel { font: bold 9px; color: #666666; }"
                                  "QLabel:disabled { color: #E4E4E4; }"};

const QString TEXT_LINK_STYLE{
    "QToolButton { font: bold 10px; color: #3d5d71; }"
    "QToolButton:disabled { color: #E4E4E4; }"
    "QToolButton:hover { color: #5b8aa8; background-color: transparent; }"};

const QString ATOM_ELEMENT_STYLE{
    "QToolButton { font: bold 14px; color: #333333; }"
    "QToolButton:disabled { color: #E4E4E4; }"};

const QString ATOM_QUERY_STYLE{
    "QLabel { font: italic 7px; color: #666666; }"
    "QToolButton { font: bold italic 14px; color: #606060; }"};

const QString PERIODIC_TABLE_STYLE{
    "#PeriodicTableForm { background-color:white }"
    "QPushButton { font: 10px; color: black; border-style: none; "
    "    min-width: 21px; max-width:21px; min-height:21px; max-height:21px; }"
    "QPushButton[class='hydrogen'] { background-color: #b2bcc2; }"
    "QPushButton[class='alkali_metals'] { background-color: #b7d9ec; }"
    "QPushButton[class='alkaline_earth_metals'] { background-color: #8fbed9; }"
    "QPushButton[class='transition_metals'] { background-color: #f2d2c6; }"
    "QPushButton[class='other_metals'] { background-color: #f2d2c6; }"
    "QPushButton[class='metalloids'] { background-color: #e1baad; }"
    "QPushButton[class='non_metals'] { background-color: #f2e8b7; }"
    "QPushButton[class='halogens'] { background-color: #f2e392; }"
    "QPushButton[class='noble_gases'] { background-color: #eccc75; }"
    "QPushButton[class='lanthanides'] { background-color: #cce5c3; }"
    "QPushButton[class='actinides'] { background-color: #afd1a2; }"
    "QLabel { font: 10px; color: black; border-style: none; padding: 0; }"};

const QString RENDER_OPTIONS_POPUP_STYLE{
    "#FileSaveImagePopup { background-color:white }"};

const QString BOND_QUERY_STYLE{
    "QToolButton { font: bold 12px; color: #333333; }"};

const QString ENUMERATION_STYLE{
    "QToolButton { font: bold 14px; color: #444444; }"};

} // namespace sketcher
} // namespace schrodinger
