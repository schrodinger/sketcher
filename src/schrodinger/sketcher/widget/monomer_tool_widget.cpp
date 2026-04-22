#include "schrodinger/sketcher/widget/monomer_tool_widget.h"

#include <array>

#include <boost/assign.hpp>

#include <QButtonGroup>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h" // ChainType
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_monomer_tool_widget.h"
#include "schrodinger/sketcher/widget/amino_acid_symbol_popup.h"
#include "schrodinger/sketcher/widget/custom_nucleotide_popup.h"
#include "schrodinger/sketcher/widget/tool_button_with_popup.h"
#include "schrodinger/sketcher/widget/modular_tool_button.h"
#include "schrodinger/sketcher/widget/nucleotide_popup.h"
#include "schrodinger/sketcher/widget/widget_utils.h"

namespace schrodinger
{
namespace sketcher
{

MonomerToolWidget::MonomerToolWidget(QWidget* parent) :
    AbstractDrawToolWidget(parent)
{
    ui.reset(new Ui::MonomerToolWidget());
    ui->setupUi(this);

    for (auto* btn_group :
         {ui->amino_monomer_group, ui->nucleic_monomer_group}) {
        for (auto* btn : btn_group->buttons()) {
            if (auto* btn_with_popup =
                    dynamic_cast<ToolButtonWithPopup*>(btn)) {
                // make sure that we call the subclass's version of
                // setStyleSheet, since it's overriden but not virtual
                btn_with_popup->setStyleSheet(ATOM_ELEMENT_OR_MONOMER_STYLE);
            } else {
                btn->setStyleSheet(ATOM_ELEMENT_OR_MONOMER_STYLE);
            }
        }
    }
    ui->unk_btn->setStyleSheet(UNKNOWN_MONOMER_STYLE);
    ui->na_n_btn->setStyleSheet(UNKNOWN_MONOMER_STYLE);
    ui->amino_monomer_btn->setStyleSheet(TEXT_LINK_STYLE);
    ui->nucleic_monomer_btn->setStyleSheet(TEXT_LINK_STYLE);

    using ButtonAminoAcidBimapType =
        boost::bimap<QAbstractButton*, AminoAcidTool>;
    using ButtonNucleicAcidBimapType =
        boost::bimap<QAbstractButton*, NucleicAcidTool>;
    // clang-format off
    m_button_amino_acid_bimap =
        boost::assign::list_of<ButtonAminoAcidBimapType::relation>
            (ui->ala_btn, AminoAcidTool::ALA)
            (ui->arg_btn, AminoAcidTool::ARG)
            (ui->asn_btn, AminoAcidTool::ASN)
            (ui->asp_btn, AminoAcidTool::ASP)
            (ui->cys_btn, AminoAcidTool::CYS)
            (ui->gln_btn, AminoAcidTool::GLN)
            (ui->glu_btn, AminoAcidTool::GLU)
            (ui->gly_btn, AminoAcidTool::GLY)
            (ui->his_btn, AminoAcidTool::HIS)
            (ui->ile_btn, AminoAcidTool::ILE)
            (ui->leu_btn, AminoAcidTool::LEU)
            (ui->lys_btn, AminoAcidTool::LYS)
            (ui->met_btn, AminoAcidTool::MET)
            (ui->phe_btn, AminoAcidTool::PHE)
            (ui->pro_btn, AminoAcidTool::PRO)
            (ui->ser_btn, AminoAcidTool::SER)
            (ui->thr_btn, AminoAcidTool::THR)
            (ui->trp_btn, AminoAcidTool::TRP)
            (ui->tyr_btn, AminoAcidTool::TYR)
            (ui->val_btn, AminoAcidTool::VAL)
            (ui->unk_btn, AminoAcidTool::UNK);
    m_button_nucleic_acid_bimap =
        boost::assign::list_of<ButtonNucleicAcidBimapType::relation>
            (ui->na_a_btn, NucleicAcidTool::A)
            (ui->na_u_btn, NucleicAcidTool::U)
            (ui->na_g_btn, NucleicAcidTool::G)
            (ui->na_c_btn, NucleicAcidTool::C)
            (ui->na_t_btn, NucleicAcidTool::T)
            (ui->na_n_btn, NucleicAcidTool::N)
            (ui->na_r_btn, NucleicAcidTool::R)
            (ui->na_dr_btn, NucleicAcidTool::dR)
            (ui->na_p_btn, NucleicAcidTool::P)
            (ui->na_rna_btn, NucleicAcidTool::RNA_NUCLEOTIDE)
            (ui->na_dna_btn, NucleicAcidTool::DNA_NUCLEOTIDE)
            (ui->na_custom_nt_btn, NucleicAcidTool::CUSTOM_NUCLEOTIDE);
    // clang-format on

    connect(ui->amino_or_nucleic_group, &QButtonGroup::buttonClicked, this,
            &MonomerToolWidget::onAminoOrNucleicBtnClicked);
    connect(ui->amino_monomer_group, &QButtonGroup::buttonClicked, this,
            &MonomerToolWidget::onAminoAcidClicked);
    connect(ui->nucleic_monomer_group, &QButtonGroup::buttonClicked, this,
            &MonomerToolWidget::onNucleicAcidClicked);

    m_rna_popup = new NucleotidePopup(NucleicAcidTool::RNA_NUCLEOTIDE,
                                      ModelKey::RNA_NUCLEOBASE, "R", "U", this);
    ui->na_rna_btn->setPopupWidget(m_rna_popup);

    m_dna_popup =
        new NucleotidePopup(NucleicAcidTool::DNA_NUCLEOTIDE,
                            ModelKey::DNA_NUCLEOBASE, "dR", "T", this);
    ui->na_dna_btn->setPopupWidget(m_dna_popup);

    m_custom_nt_popup = new CustomNucleotidePopup(this);
    ui->na_custom_nt_btn->setPopupWidget(m_custom_nt_popup);

    // Set up amino acid analog popups from the monomer database
    static const std::unordered_map<std::string, std::string>
        STANDARD_AA_NAMES = {
            {"A", "Alanine"},    {"R", "Arginine"},      {"N", "Asparagine"},
            {"D", "Aspartate"},  {"C", "Cysteine"},      {"Q", "Glutamine"},
            {"E", "Glutamate"},  {"G", "Glycine"},       {"H", "Histidine"},
            {"I", "Isoleucine"}, {"L", "Leucine"},       {"K", "Lysine"},
            {"M", "Methionine"}, {"F", "Phenylalanine"}, {"P", "Proline"},
            {"S", "Serine"},     {"T", "Threonine"},     {"W", "Tryptophan"},
            {"Y", "Tyrosine"},   {"V", "Valine"},        {"X", "Unknown"},
        };

    auto analogs_by_aa =
        rdkit_extensions::MonomerDatabase::instance()
            .getMonomersByNaturalAnalog(rdkit_extensions::ChainType::PEPTIDE);

    for (auto& entry : m_button_amino_acid_bimap.left) {
        auto* button = entry.first;
        auto aa_tool = entry.second;
        auto symbol = AMINO_ACID_TOOL_TO_RES_NAME.at(aa_tool);
        auto it = analogs_by_aa.find(symbol);
        if (it == analogs_by_aa.end() || it->second.empty()) {
            continue;
        }
        auto name_it = STANDARD_AA_NAMES.find(symbol);
        auto name = name_it != STANDARD_AA_NAMES.end() ? name_it->second : "";
        auto* popup = new AminoAcidSymbolPopup(symbol, name, it->second, this);
        auto* modular_btn = qobject_cast<ModularToolButton*>(button);
        Q_ASSERT_X(modular_btn, "MonomerToolWidget",
                   "Expected ModularToolButton");
        modular_btn->setPopupWidget(popup);
        modular_btn->setEnumItem(0); // default to standard AA
        modular_btn->showPopupIndicator(false);
        m_amino_acid_symbol_popups[button] = popup;
    }
}

MonomerToolWidget::~MonomerToolWidget() = default;

void MonomerToolWidget::setModel(SketcherModel* model)
{
    AbstractDrawToolWidget::setModel(model);
    m_rna_popup->setModel(model);
    m_dna_popup->setModel(model);
    m_custom_nt_popup->setModel(model);
    for (auto& [button, popup] : m_amino_acid_symbol_popups) {
        popup->setModel(model);
    }
    updateCheckedButton();
}

void MonomerToolWidget::updateWidgetsEnabled()
{
    auto model = getModel();
    if (model == nullptr) {
        return;
    }
    auto has_peptide = model->hasPeptideSelection();
    auto has_na_base = model->hasNABaseSelection();
    auto has_na_sugar = model->hasNASugarSelection();
    auto has_na_phosphate = model->hasNAPhosphateSelection();
    bool has_any_na = has_na_base || has_na_sugar || has_na_phosphate;
    bool has_monomer_selection = has_peptide || has_any_na;

    // Apply selection active style to amino/nucleic pages
    ui->amino_page->setStyleSheet(has_peptide ? SELECTION_ACTIVE_STYLE : "");
    ui->nucleic_page->setStyleSheet(has_any_na ? SELECTION_ACTIVE_STYLE : "");

    // Enable amino acid buttons only if no monomer selection or peptides
    // are selected
    for (auto* btn : ui->amino_monomer_group->buttons()) {
        btn->setEnabled(!has_monomer_selection || has_peptide);
    }

    // Enable individual nucleic acid component buttons only if that
    // component type is in the selection (or there's no monomer selection)
    std::array base_btns{ui->na_a_btn, ui->na_c_btn, ui->na_g_btn,
                         ui->na_u_btn, ui->na_t_btn, ui->na_n_btn};
    for (auto* btn : base_btns) {
        btn->setEnabled(!has_monomer_selection || has_na_base);
    }
    ui->na_r_btn->setEnabled(!has_monomer_selection || has_na_sugar);
    ui->na_dr_btn->setEnabled(!has_monomer_selection || has_na_sugar);
    ui->na_p_btn->setEnabled(!has_monomer_selection || has_na_phosphate);

    // Full nucleotide buttons are never available for mutation
    ui->na_rna_btn->setEnabled(!has_monomer_selection);
    ui->na_dna_btn->setEnabled(!has_monomer_selection);
    ui->na_custom_nt_btn->setEnabled(!has_monomer_selection);
}

std::unordered_set<QAbstractButton*> MonomerToolWidget::getCheckableButtons()
{
    auto buttons = AbstractDrawToolWidget::getCheckableButtons();
    for (auto group : {ui->amino_monomer_group, ui->nucleic_monomer_group}) {
        for (auto button : group->buttons()) {
            buttons.insert(button);
        }
    }
    return buttons;
}

void MonomerToolWidget::updateCheckedButton()
{
    auto model = getModel();
    if (model == nullptr) {
        return;
    }
    QAbstractButton* amino_button = nullptr;
    QAbstractButton* nucleic_button = nullptr;
    bool has_sel = model->hasActiveSelection();
    if (model->getDrawTool() == DrawTool::MONOMER) {
        if (model->getMonomerToolType() == MonomerToolType::AMINO_ACID) {
            ui->amino_monomer_btn->setChecked(true);
            ui->amino_or_nucleic_stack->setCurrentWidget(ui->amino_page);
            if (!has_sel) {
                auto amino_acid = model->getAminoAcidTool();
                amino_button = m_button_amino_acid_bimap.right.at(amino_acid);
            }
        } else {
            ui->nucleic_monomer_btn->setChecked(true);
            ui->amino_or_nucleic_stack->setCurrentWidget(ui->nucleic_page);
            if (!has_sel) {
                auto nucleic_acid = model->getNucleicAcidTool();
                nucleic_button =
                    m_button_nucleic_acid_bimap.right.at(nucleic_acid);
            }
        }
    }
    check_button_or_uncheck_group(amino_button, ui->amino_monomer_group);
    check_button_or_uncheck_group(nucleic_button, ui->nucleic_monomer_group);

    // Only show the popup indicator arrow on the currently checked button
    for (auto& [btn, popup] : m_amino_acid_symbol_popups) {
        auto* modular_btn = qobject_cast<ModularToolButton*>(btn);
        Q_ASSERT_X(modular_btn, "MonomerToolWidget",
                   "Expected ModularToolButton");
        modular_btn->showPopupIndicator(btn == amino_button);
    }

    // sync amino acid analog popup state
    auto analog = model->getValueString(ModelKey::AMINO_ACID_SYMBOL);
    for (auto& [btn, popup] : m_amino_acid_symbol_popups) {
        auto* modular_btn = qobject_cast<ModularToolButton*>(btn);
        Q_ASSERT_X(modular_btn, "MonomerToolWidget",
                   "Expected ModularToolButton");
        // find the ID matching the current analog symbol
        for (auto& packet : popup->getButtonPackets()) {
            if (popup->getSymbolForId(packet.enum_int) == analog) {
                modular_btn->setEnumItem(packet.enum_int);
                break;
            }
        }
    }

    ui->na_rna_btn->setEnumItem(static_cast<int>(model->getRNANucleobase()));
    ui->na_dna_btn->setEnumItem(static_cast<int>(model->getDNANucleobase()));

    // update the text displayed on the custom nucleotide button
    auto [sugar, base, phosphate] = model->getCustomNucleotide();
    const QString custom_nt_name_fmt("%1(%2)%3");
    auto custom_nt_name = custom_nt_name_fmt.arg(sugar, base, phosphate);
    ui->na_custom_nt_btn->setText(custom_nt_name);
}

void MonomerToolWidget::onAminoOrNucleicBtnClicked(QAbstractButton* button)
{
    QWidget* page = nullptr;
    MonomerToolType tool_type;
    if (button == ui->amino_monomer_btn) {
        page = ui->amino_page;
        tool_type = MonomerToolType::AMINO_ACID;
    } else {
        page = ui->nucleic_page;
        tool_type = MonomerToolType::NUCLEIC_ACID;
    }
    ui->amino_or_nucleic_stack->setCurrentWidget(page);
    getModel()->setValues(
        {{ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::MONOMER)},
         {ModelKey::MONOMER_TOOL_TYPE, QVariant::fromValue(tool_type)}});
    updateCheckedButton();
}

/**
 * When a tool button is clicked, update the model with the appropriate value
 * for the button
 */
template <typename T> static void
on_tool_clicked(SketcherModel* model, const ModelKey key,
                const boost::bimap<QAbstractButton*, T>& button_tool_bimap,
                QAbstractButton* button)
{
    auto tool = button_tool_bimap.left.at(button);
    ping_or_set_model_value(model, key, tool);
}

/**
 * If the model has a selection, ping the specified key/value.  Otherwise, set
 * the key to the value and ensure that the draw tool is set to MONOMER.
 */
template <typename T> static void
ping_or_set_model_value(SketcherModel* model, const ModelKey key, const T value)
{
    if (model->hasActiveSelection()) {
        // ping the model to indicate that we want to replace the selection
        // without changing the tool
        model->pingValue(key, value);
    } else {
        // change the tool
        model->setValues(
            {{ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::MONOMER)},
             {key, QVariant::fromValue(value)}});
    }
}

void MonomerToolWidget::onAminoAcidClicked(QAbstractButton* button)
{
    on_tool_clicked<AminoAcidTool>(getModel(), ModelKey::AMINO_ACID_TOOL,
                                   m_button_amino_acid_bimap, button);

    // Also set which specific analog is selected
    QString symbol;
    auto it = m_amino_acid_symbol_popups.find(button);
    if (it != m_amino_acid_symbol_popups.end()) {
        auto* modular_btn = qobject_cast<ModularToolButton*>(button);
        Q_ASSERT_X(modular_btn, "MonomerToolWidget",
                   "Expected ModularToolButton");
        symbol = it->second->getSymbolForId(modular_btn->getEnumItem());
    } else {
        auto tool = m_button_amino_acid_bimap.left.at(button);
        symbol = QString::fromStdString(AMINO_ACID_TOOL_TO_RES_NAME.at(tool));
    }
    ping_or_set_model_value(getModel(), ModelKey::AMINO_ACID_SYMBOL, symbol);
}

void MonomerToolWidget::onNucleicAcidClicked(QAbstractButton* button)
{
    on_tool_clicked<NucleicAcidTool>(getModel(), ModelKey::NUCLEIC_ACID_TOOL,
                                     m_button_nucleic_acid_bimap, button);
    // if one of the full nucleotide buttons was clicked, we also need to inform
    // the model of the nucleotide type
    if (button == ui->na_rna_btn || button == ui->na_dna_btn) {
        auto* modular_button = static_cast<ModularToolButton*>(button);
        auto base = static_cast<StdNucleobase>(modular_button->getEnumItem());
        auto key = button == ui->na_rna_btn ? ModelKey::RNA_NUCLEOBASE
                                            : ModelKey::DNA_NUCLEOBASE;
        ping_or_set_model_value(getModel(), key, base);
    }
}

} // namespace sketcher
} // namespace schrodinger
