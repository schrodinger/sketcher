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
    "QToolButton:disabled { color: #EEEEEE; }"
    "QToolButton:checked { background-color: #d4e6f1; }"
    "QToolButton:hover:!checked { background-color: #edf7fc; }"};

const QString TOOL_BUTTON_WITH_POPUP_STYLE{
    "QToolButton { border: none; padding: 0px; margin: 0px; spacing: 0px; }" +
    TOOL_BUTTON_STYLE};

const QString TOOL_BUTTON_CORNER_ARROW_STYLE{
    "QToolButton {"
    " background-image:url(':icons/menu_corner_arrow.svg'); "
    " background-position: bottom right; background-repeat: no-repeat; }"
    "QToolButton:disabled {"
    " background-image:url(':icons/menu_corner_arrow_dis.svg'); }"};

const QString PALETTE_TITLE_STYLE{"QLabel { font: bold 9px; color: #666666; }"
                                  "QLabel:disabled { color: #EEEEEE; }"};

const QString SELECT_ACTIONS_STYLE{
    "QToolButton { font: bold 10px; color: #3d5d71; }"
    "QToolButton:hover { color: #5b8aa8; background-color: transparent; }" +
    TOOL_BUTTON_STYLE};

const QString ATOM_ELEMENT_STYLE{
    "QToolButton { font: bold 14px; color: #333333; }" + TOOL_BUTTON_STYLE};

const QString ATOM_QUERY_STYLE{
    "QLabel { font: italic 7px; color: #666666; }"
    "QToolButton { font: bold italic 14px; color: #606060; }" +
    TOOL_BUTTON_STYLE};

const QString PERIODIC_TABLE_STYLE{
    "#PeriodicTableForm { background-color:white }"
    "QPushButton { font: 10px; color: black; border-style: none; }"
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
    "QLabel { font: 10px; color: black; border-style: none; }"};

const QString BOND_QUERY_STYLE{
    "QToolButton { font: bold 12px; color: #333333; }" + TOOL_BUTTON_STYLE};

const QString ENUMERATION_STYLE{
    "QToolButton { font: bold 14px; color: #444444; }" + TOOL_BUTTON_STYLE};

} // namespace sketcher
} // namespace schrodinger
