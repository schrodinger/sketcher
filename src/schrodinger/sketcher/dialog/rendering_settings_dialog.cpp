#include "schrodinger/sketcher/dialog/rendering_settings_dialog.h"

#include <QPushButton>
#include <QTimer>

#include "mmshare-version.h"
#include "schrodinger/sketcher/ui/ui_rendering_settings_dialog.h"

Q_DECLARE_METATYPE(schrodinger::sketcher::ColorScheme);

namespace schrodinger
{
namespace sketcher
{

struct RenderingSettings {
    int m_atom_font_size = DEFAULT_FONT_SIZE;
    int m_bond_line_width = BOND_DEFAULT_PEN_WIDTH;
    CarbonLabels m_carbon_labels = CarbonLabels::NONE;
    bool m_color_heteroatoms = true;
    ColorScheme m_color_scheme = ColorScheme::DEFAULT;
    bool m_show_stereo_annotations = true;
    bool m_show_unknown_stereo_annotations = true;
    bool m_explicit_abs_labels_shown = false;
};

RenderingSettingsDialog::RenderingSettingsDialog(SketcherModel* model,
                                                 QWidget* parent) :
    ModalDialog(parent)
{
    // ModalDialog sets WA_DeleteOnClose to true, but we want to keep this class
    // when the dialog is closed, so we override it here
    setAttribute(Qt::WA_DeleteOnClose, false);
    m_ui.reset(new Ui::RenderingSettingsDialog());
    setupDialogUI(*m_ui);
    initColorModes();
    m_sketcher_model = model;
    m_update_timer = new QTimer(this);
    m_update_timer->setSingleShot(true);
    m_update_timer->setInterval(0);
    connect(m_update_timer, &QTimer::timeout, this,
            &RenderingSettingsDialog::onValuesChanged);
    connect_input_widgets_to_timer(this, m_update_timer);
    connect(m_ui->m_reset_to_default_pb, &QPushButton::clicked, this,
            &RenderingSettingsDialog::loadDefaults);
    connect(model, &SketcherModel::valuesChanged, this,
            &RenderingSettingsDialog::loadSettingsFromModel);

    // initialize the GUI
    loadSettingsFromModel();
}

RenderingSettingsDialog::~RenderingSettingsDialog() = default;

void RenderingSettingsDialog::onValuesChanged()
{
    updateWidgets();
    exportSettingsToModel();
}

void RenderingSettingsDialog::initColorModes()
{
    m_ui->m_color_mode_combo->addItem(
        "Default", QVariant::fromValue(ColorScheme::DEFAULT));
    m_ui->m_color_mode_combo->addItem("Avalon",
                                      QVariant::fromValue(ColorScheme::AVALON));
    m_ui->m_color_mode_combo->addItem("CDK",
                                      QVariant::fromValue(ColorScheme::CDK));
    m_ui->m_color_mode_combo->addItem(
        "Dark", QVariant::fromValue(ColorScheme::DARK_MODE));

    m_ui->m_bw_mode_combo->addItem(
        "Default", QVariant::fromValue(ColorScheme::BLACK_WHITE));
    m_ui->m_bw_mode_combo->addItem(
        "Dark", QVariant::fromValue(ColorScheme::WHITE_BLACK));
}

void RenderingSettingsDialog::loadSettingsFromModel()
{
    if (m_sketcher_model == nullptr) {
        return;
    }
    auto settings = getSettingsFromModel(m_sketcher_model);
    loadSettings(settings);
}

void RenderingSettingsDialog::updateWidgets()
{
    bool check_label_carbons = m_ui->m_label_carbons_cb->isChecked();
    m_ui->m_label_all_C_rb->setEnabled(check_label_carbons);
    m_ui->m_label_terminal_C_rb->setEnabled(check_label_carbons);

    bool show_stereo = m_ui->m_show_stereo_cb->isChecked();
    m_ui->m_with_ABS_rb->setEnabled(show_stereo);
    m_ui->m_without_ABS_rb->setEnabled(show_stereo);
    m_ui->m_undefined_centers_labels_cb->setEnabled(show_stereo);

    auto sync_comboboxes = [](QComboBox* visible_combobox,
                              QComboBox* hidden_combobox) {
        bool dark_mode = visible_combobox->currentText() == "Dark";
        hidden_combobox->setCurrentText(dark_mode ? "Dark" : "Default");
    };

    if (m_ui->m_color_heteroatoms_cb->isChecked()) {
        sync_comboboxes(m_ui->m_color_mode_combo, m_ui->m_bw_mode_combo);
    } else {
        sync_comboboxes(m_ui->m_bw_mode_combo, m_ui->m_color_mode_combo);
    }
}

void RenderingSettingsDialog::loadSettings(RenderingSettings& settings)
{
    m_ui->m_atom_font_size_sb->setValue(settings.m_atom_font_size);
    m_ui->m_bond_line_width_sb->setValue(settings.m_bond_line_width);
    m_ui->m_label_carbons_cb->setChecked(settings.m_carbon_labels !=
                                         CarbonLabels::NONE);
    m_ui->m_label_all_C_rb->setChecked(settings.m_carbon_labels ==
                                       CarbonLabels::ALL);
    m_ui->m_label_terminal_C_rb->setChecked(settings.m_carbon_labels ==
                                            CarbonLabels::TERMINAL);

    m_ui->m_color_heteroatoms_cb->setChecked(settings.m_color_heteroatoms);

    auto update_combo_boxes = [settings](QComboBox* combobox_to_show,
                                         QComboBox* combobox_to_hide) {
        combobox_to_hide->hide();
        combobox_to_show->show();
        auto index = combobox_to_show->findData(
            QVariant::fromValue(settings.m_color_scheme));
        if (index != -1) {
            combobox_to_show->setCurrentIndex(index);
        }
    };
    if (settings.m_color_heteroatoms) {
        update_combo_boxes(m_ui->m_color_mode_combo, m_ui->m_bw_mode_combo);

    } else {
        update_combo_boxes(m_ui->m_bw_mode_combo, m_ui->m_color_mode_combo);
    }
    m_ui->m_show_stereo_cb->setChecked(settings.m_show_stereo_annotations);

    m_ui->m_with_ABS_rb->setChecked(settings.m_explicit_abs_labels_shown);
    m_ui->m_without_ABS_rb->setChecked(!settings.m_explicit_abs_labels_shown);
    m_ui->m_undefined_centers_labels_cb->setChecked(
        settings.m_show_unknown_stereo_annotations);
    updateWidgets();
}

RenderingSettings RenderingSettingsDialog::getSettingsFromModel(
    const sketcher::SketcherModel* model)
{
    if (model == nullptr) {
        throw std::runtime_error(
            "RenderingSettingsDialog::getSettingsFromModel: model is null");
    }
    RenderingSettings settings;
    settings.m_atom_font_size = model->getFontSize();
    settings.m_bond_line_width =
        model->getBondDisplaySettingsPtr()->m_bond_width;
    settings.m_carbon_labels =
        model->getAtomDisplaySettingsPtr()->m_carbon_labels;
    settings.m_color_heteroatoms =
        model->getValueBool(ModelKey::COLOR_HETEROATOMS);
    settings.m_color_scheme = model->getColorScheme();
    settings.m_show_stereo_annotations =
        model->getValueBool(ModelKey::SHOW_STEREOCENTER_LABELS);
    settings.m_explicit_abs_labels_shown =
        model->getAtomDisplaySettingsPtr()->m_explicit_abs_labels_shown;

    auto stereo_labels_visibility =
        model->getAtomDisplaySettingsPtr()->m_stereo_labels_visibility;
    if (stereo_labels_visibility != StereoLabels::NONE) {
        settings.m_show_unknown_stereo_annotations =
            stereo_labels_visibility == StereoLabels::ALL;
    }
    return settings;
}

void RenderingSettingsDialog::loadDefaults()
{
    RenderingSettings settings;
    loadSettings(settings);
    exportSettingsToModel();
}

RenderingSettings RenderingSettingsDialog::getSettingsFromPanel()
{
    RenderingSettings settings;
    settings.m_atom_font_size = m_ui->m_atom_font_size_sb->value();
    settings.m_bond_line_width = m_ui->m_bond_line_width_sb->value();
    if (!m_ui->m_label_carbons_cb->isChecked()) {
        settings.m_carbon_labels = CarbonLabels::NONE;
    } else {
        if (m_ui->m_label_terminal_C_rb->isChecked()) {
            settings.m_carbon_labels = CarbonLabels::TERMINAL;
        } else if (m_ui->m_label_all_C_rb->isChecked()) {
            settings.m_carbon_labels = CarbonLabels::ALL;
        }
    }
    settings.m_color_heteroatoms = m_ui->m_color_heteroatoms_cb->isChecked();
    settings.m_color_scheme = qvariant_cast<schrodinger::sketcher::ColorScheme>(
        (settings.m_color_heteroatoms) ? m_ui->m_color_mode_combo->currentData()
                                       : m_ui->m_bw_mode_combo->currentData());

    settings.m_show_stereo_annotations = m_ui->m_show_stereo_cb->isChecked();

    settings.m_explicit_abs_labels_shown = m_ui->m_with_ABS_rb->isChecked();

    settings.m_show_unknown_stereo_annotations =
        m_ui->m_undefined_centers_labels_cb->isChecked();

    return settings;
}

void RenderingSettingsDialog::exportSettingsToModel()
{
    if (m_sketcher_model == nullptr) {
        return;
    }
    auto settings = getSettingsFromPanel();
    m_sketcher_model->setFontSize(settings.m_atom_font_size);
    bool dark_background = settings.m_color_scheme == ColorScheme::DARK_MODE ||
                           settings.m_color_scheme == ColorScheme::WHITE_BLACK;
    m_sketcher_model->setBackgroundColor(
        dark_background ? DARK_BACKGROUND_COLOR : LIGHT_BACKGROUND_COLOR);

    m_sketcher_model->setColorScheme(settings.m_color_scheme);

    std::unordered_map<ModelKey, QVariant> kv_pairs = {
        {ModelKey::COLOR_HETEROATOMS,
         QVariant::fromValue(settings.m_color_heteroatoms)},
        {ModelKey::SHOW_STEREOCENTER_LABELS,
         QVariant::fromValue(settings.m_show_stereo_annotations)}};
    m_sketcher_model->setValues(kv_pairs);

    AtomDisplaySettings atom_settings(
        *(m_sketcher_model->getAtomDisplaySettingsPtr()));
    atom_settings.m_carbon_labels = settings.m_carbon_labels;

    atom_settings.m_explicit_abs_labels_shown =
        settings.m_explicit_abs_labels_shown;

    atom_settings.m_stereo_labels_visibility =
        settings.m_show_stereo_annotations
            ? (settings.m_show_unknown_stereo_annotations ? StereoLabels::ALL
                                                          : StereoLabels::KNOWN)
            : StereoLabels::NONE;

    m_sketcher_model->setAtomDisplaySettings(atom_settings);

    BondDisplaySettings bond_settings(
        *(m_sketcher_model->getBondDisplaySettingsPtr()));
    bond_settings.m_bond_width = settings.m_bond_line_width;
    bond_settings.m_color =
        atom_settings.getAtomColor(static_cast<int>(Element::C));
    m_sketcher_model->setBondDisplaySettings(bond_settings);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/rendering_settings_dialog.moc"
