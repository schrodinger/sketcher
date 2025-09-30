#pragma once
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <QGraphicsItem>
#include <QObject>
#include <QVariant>
#include <boost/functional/hash.hpp>

#include "schrodinger/sketcher/image_constants.h"
#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/atom_display_settings.h"
#include "schrodinger/sketcher/molviewer/bond_display_settings.h"

class QPointF;

namespace RDKit
{
class Atom;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

struct RenderOptions;

/**
 * Whether to the entire scene or a selection subset should be used
 */
enum class SceneSubset { ALL, SELECTION, HOVERED, SELECTED_OR_HOVERED };

/**
 * Possible selection tools to equip
 */
enum class SelectionTool {
    RECTANGLE,
    ELLIPSE,
    LASSO,
    FRAGMENT,
};

/**
 * All supported elements.
 */
enum class Element {
    H = 1,
    HE,
    LI,
    BE,
    B,
    C,
    N,
    O,
    F,
    NE,
    NA,
    MG,
    AL,
    SI,
    P,
    S,
    CL,
    AR,
    K,
    CA,
    SC,
    TI,
    V,
    CR,
    MN,
    FE,
    CO,
    NI,
    CU,
    ZN,
    GA,
    GE,
    AS,
    SE,
    BR,
    KR,
    RB,
    SR,
    Y,
    ZR,
    NB,
    MO,
    TC,
    RU,
    RH,
    PD,
    AG,
    CD,
    IN_, // This name must be mangled to avoid Windows compilation failure
    SN,
    SB,
    TE,
    I,
    XE,
    CS,
    BA,
    LA,
    CE,
    PR,
    ND,
    PM,
    SM,
    EU,
    GD,
    TB,
    DY,
    HO,
    ER,
    TM,
    YB,
    LU,
    HF,
    TA,
    W,
    RE,
    OS,
    IR,
    PT,
    AU,
    HG,
    TL,
    PB,
    BI,
    PO,
    AT,
    RN,
    FR,
    RA,
    AC,
    TH,
    PA,
    U,
    NP,
    PU,
    AM,
    CM,
    BK,
    CF,
    ES,
    FM,
    MD,
    NO,
    LR,
    RF,
    DB,
    SG,
    BH,
    HS,
    MT,
    DS,
    RG,
    CN,
    NH,
    FL,
    MC,
    LV,
    TS,
    OG,
};

/**
 * Atom queries.
 */
enum class AtomQuery {
    A,
    AH,
    Q,
    QH,
    M,
    MH,
    X,
    XH,
};

/**
 * Atom-related tools.
 *
 * For full specification of tool, see related enums `Element` and `AtomQuery`.
 */
enum class AtomTool {
    ELEMENT,
    QUERY,
};

/**
 * Possible bond tools.
 */
enum class BondTool {
    SINGLE,
    DOUBLE,
    TRIPLE,
    COORDINATE,
    ZERO,
    SINGLE_UP,
    SINGLE_DOWN,
    SINGLE_EITHER,
    DOUBLE_EITHER,
    AROMATIC,
    SINGLE_OR_DOUBLE,
    SINGLE_OR_AROMATIC,
    DOUBLE_OR_AROMATIC,
    ANY,
    ATOM_CHAIN
};

/**
 * Charge tools
 */
enum class ChargeTool {
    DECREASE,
    INCREASE,
};

/**
 * Ring fragment tools.
 */
enum class RingTool {
    CYCLOPROPANE,
    CYCLOBUTANE,
    CYCLOPENTANE,
    CYCLOPENTADIENE,
    CYCLOHEXANE,
    BENZENE,
    CYCLOHEPTANE,
    CYCLOOCTANE,
};

/**
 * Enumeration tools.
 */
enum class EnumerationTool : int {
    NEW_RGROUP,
    EXISTING_RGROUP,
    ATTACHMENT_POINT,
    ADD_MAPPING,
    REMOVE_MAPPING,
    RXN_ARROW,
    RXN_PLUS,
};

/**
 * Possible draw tools.
 *
 * The "atom" tool will require the atom type to be stored on the model
 * separately. This is separate to avoid enumerating each atom type in this enum
 * class. Likewise, other draw tools will have their own subclasses which will
 * be stored elsewhere in the model.
 */
enum class DrawTool {
    NO_INTERACTION,
    SELECT,
    MOVE_ROTATE,
    ERASE,
    ATOM,
    BOND,
    CHARGE,
    RING,
    ENUMERATION,
    RESIDUE,
    EXPLICIT_H,
};

/**
 * Different types of data that the model is responsible for.
 */
enum class ModelKey {
    NEW_STRUCTURES_REPLACE_CONTENT,
    SHOW_LID_LEGEND,
    ALLOW_MULTIPLE_RXNS,
    SHOW_VALENCE_ERRORS,
    COLOR_HETEROATOMS,
    SHOW_STEREO_LABELS,
    USE_IMPLICIT_HYDROGENS,
    SELECTION_TOOL,
    DRAW_TOOL,
    ATOM_TOOL,
    BOND_TOOL,
    CHARGE_TOOL,
    RING_TOOL,
    ENUMERATION_TOOL,
    ELEMENT,
    ATOM_QUERY,
    RGROUP_NUMBER,
    RESIDUE_TYPE,
};

/**
 * @return Every element of `ModelKey`.
 */
SKETCHER_API std::vector<ModelKey> get_model_keys();

/**
 * Model for shared state among the various sketcher views.
 *
 * In order to reduce coupling, any state that needs to be accessed by more than
 * one sketcher view should be stored on this model. By connecting this model to
 * all views that need access to this information, we can ensure that there is a
 * single source of truth for shared data and that all views can be immediately
 * synchronized if this data should be updated by any one of them.
 */
class SKETCHER_API SketcherModel : public QObject
{
    Q_OBJECT

  public:
    SketcherModel(QObject* parent = nullptr);

    /**
     * Retrieve arbitrary data from the model's state map.
     *
     * @param key A key corresponding to the desired state value
     * @return The value associated with the specified key
     */
    QVariant getValue(ModelKey key) const;

    /**
     * Convenience wrappers for model keys with an unique enum class value
     *
     * @return value for the key associated with the enum class
     */
    SelectionTool getSelectionTool() const;
    DrawTool getDrawTool() const;
    AtomTool getAtomTool() const;
    BondTool getBondTool() const;
    ChargeTool getChargeTool() const;
    RingTool getRingTool() const;
    EnumerationTool getEnumerationTool() const;
    Element getElement() const;
    AtomQuery getAtomQuery() const;

    /**
     * Retrieve data from the model's state map.
     *
     * @param key A key corresponding to the desired state value
     * @return The value associated with the specified key
     */
    bool getValueBool(ModelKey key) const;
    int getValueInt(ModelKey key) const;
    QString getValueString(ModelKey key) const;

    /**
     * @param key model key to ping without setting its value in the model
     * @param value value to ping the model with
     */
    template <class T> void pingValue(ModelKey key, T value)
    {
        emit valuePinged(key, QVariant::fromValue(value));
    }

    /**
     * @param key model key to set
     * @param value value to assign to the model
     */
    template <class T> void setValue(ModelKey key, T value)
    {
        setValues({{key, QVariant::fromValue(value)}});
    }

    /**
     * @param key_value_map (key, value) pairs to assign to the model
     */
    void setValues(const std::unordered_map<ModelKey, QVariant>& key_value_map);

    /**
     * @return whether a reaction is present
     */
    bool hasReaction() const;

    /**
     * When a reaction arrow is added, disable the reaction arrow tool unless
     * ALLOW_MULTIPLE_RXNS is enabled.
     */
    void onReactionArrowAdded();

    /**
     * Query the attached scene (if any) and return R group numbers.
     *
     * @return the R group numbers currently in use in the associated scene
     */
    std::set<unsigned int> getRGroupNumbers() const;

    /**
     * @return the smallest finite integer not currently used as an R group
     * number in the associated scene
     */
    unsigned int getNextRGroupNumber() const;

    /**
     * @return Whether the scene is empty
     */
    bool sceneIsEmpty() const;

    /**
     * @return Whether there are any items selected in the scene.
     */
    bool hasActiveSelection() const;
    virtual bool hasAtomSelection() const;
    virtual bool hasBondSelection() const;
    virtual bool hasNonMolecularObjectSelection() const;

    /**
     * @return Whether all items in the scene are selected.
     */
    bool allItemsSelected() const;

    /**
     * @return all interactive graphics items in the scene.
     */
    QList<QGraphicsItem*> getInteractiveItems() const;

    /**
     * @return all selected items in the scene.
     */
    QList<QGraphicsItem*> getSelection() const;

    /**
     * @return all objects associated with the currently-open context
     * menu (if any)
     */
    QList<QGraphicsItem*> getContextMenuObjects() const;

    /**
     * @return Whether an undo and redo can occur
     */
    std::pair<bool, bool> getUndoStackData() const;

    const AtomDisplaySettings* getAtomDisplaySettingsPtr() const;
    const BondDisplaySettings* getBondDisplaySettingsPtr() const;
    void setAtomDisplaySettings(const AtomDisplaySettings& settings);
    void setBondDisplaySettings(const BondDisplaySettings& settings);
    void setBondWidthScale(qreal scale);
    void setFontSize(int size);
    int getFontSize() const;

    /**
     * Update the color scheme, which controls colors for atoms and bonds, as
     * well as the background color
     */
    void setColorScheme(const ColorScheme scheme);

    /**
     * @return the currently active color scheme
     */
    ColorScheme getColorScheme() const;

    /**
     * @return a std::pair containing the color schemes for color and black &
     * white modes. Which one is actually in use depends on the state of
     * COLOR_HETEROATOMS
     */
    std::pair<ColorScheme, ColorScheme> getColorSchemes() const;

    QColor getBackgroundColor() const;

    /**
     * @return whether the current color scheme uses a dark background
     */
    bool hasDarkColorScheme() const;

    /**
     * @param opts render options to load into the model
     */
    void loadRenderOptions(const RenderOptions& opts);

    /**
     * By default, the select, move-rotate and delete tool will be switched to
     * the draw carbon tool when the scene is empty. This method allows that
     * behavior to be disabled, which is required when the Sketcher is in
     * select-only mode.
     */
    void setSelectToolAllowedWhenSceneEmpty(const bool allowed);

  signals:
    void valuePinged(ModelKey key, QVariant value);
    void valuesChanged(const std::unordered_set<ModelKey>& key);

    /**
     * Signal to notify the scene that the display settings have changed. This
     * triggers an update of the scene's atom and bond items so they can be
     * displayed with the new settings.
     */
    void displaySettingsChanged() const;

    /**
     * Signal to notify the scene that the background color has changed.
     */
    void backgroundColorChanged(const QColor& color);

    /**
     * Signal to notify views when the interactive items of the scene have
     * changed.
     *
     * Scene information is not stored on this model, so the responsibilty for
     * emitting this signal must go to the associated scene.
     */
    void interactiveItemsChanged() const;

    /**
     * Signal used to notify views that the scene selection has changed.
     *
     * Selection information is not stored on this model, so the responsibility
     * for emitting this signal must go to the associated scene.
     */
    void selectionChanged() const;

    /**
     * @return interactive items from the associated scene
     */
    QList<QGraphicsItem*> interactiveItemsRequested() const;

    /**
     * @return selection items from the associated scene
     */
    QList<QGraphicsItem*> selectionRequested() const;

    /**
     * Signal used to request context menu object information from the
     * associated scene. These items are the subject of the currently-open
     * context menu. If there is not a context menu open or if the context menu
     * was opened without being associated with any scene objects, the set
     * should remain empty.
     *
     * @return context menu information from the associated scene
     */
    QList<QGraphicsItem*> contextMenuObjectsRequested() const;

    /**
     * Signal used to request context menu object information from the
     * associated scene. These items are the subject of the currently-open
     * context menu. If there is not a context menu open or if the context menu
     * was opened without being associated with any scene objects, the set
     * should remain empty.
     *
     * @return context menu information from the associated scene
     */
    std::vector<const RDKit::Atom*> contextMenuAtomsRequested() const;

    /**
     * @return the number of reaction arrows in the scene
     */
    unsigned int reactionCountRequested() const;

    /**
     * @return Whether an undo operation can be performed
     */
    bool undoStackCanUndoRequested() const;

    /**
     * @return Whether a redo operation can be performed
     */
    bool undoStackCanRedoRequested() const;

    /**
     * Signal used to notify views that the scene undo stack has changed.
     *
     * The undo stack is not stored on this model, so the responsibility for
     * emitting this signal must go to the associated scene.
     */
    void undoStackDataChanged() const;

    /**
     * Signal used to request R group numbers in use in the associated scene.
     *
     * @param rgroup_numbers A empty set of numbers that should be modified in
     * place by the associated sketcher scene, if any.
     */
    void rGroupNumbersRequested(std::set<unsigned int>& rgroup_numbers) const;

    /**
     * Signal used to notify views that the extant R groups have changed.
     *
     * R groups are not stored on this model, so the responsibility for
     * emitting this signal must go to the associated scene.
     */
    void rGroupNumbersChanged() const;

  protected:
    /**
     * Update active tool based on change in selection
     */
    void onSelectionChanged();
    /**
     * Update active tool when items are deleted
     */
    void onInteractiveItemsChanged();

    /**
     * Update the atom display settings when the model is updated.
     * @param keys The keys that have been updated.
     */
    void updateAtomDisplaySettings(const std::unordered_set<ModelKey>& keys);

    void setBackgroundColor(QColor color);

    std::unordered_map<ModelKey, QVariant> m_model_map;

    AtomDisplaySettings m_atom_display_settings;
    BondDisplaySettings m_bond_display_settings;
    int m_font_size = DEFAULT_FONT_SIZE;
    QColor m_background_color = LIGHT_BACKGROUND_COLOR;
    bool m_allow_select_tool_when_scene_empty = false;

    ColorScheme m_color_scheme = ColorScheme::DEFAULT;
};

} // namespace sketcher
} // namespace schrodinger