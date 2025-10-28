#pragma once
#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/sketcher_view_with_white_background.h"

namespace Ui
{
class CustomNucleotidePopup;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * Popup used to provide sugar, base, and phosphate entry for the custom
 * nucleotide
 */
class SKETCHER_API CustomNucleotidePopup
    : public SketcherViewWithWhiteBackground
{
  public:
    CustomNucleotidePopup(QWidget* parent = nullptr);
    ~CustomNucleotidePopup();

    void setModel(SketcherModel* model) override;

  protected:
    std::unique_ptr<Ui::CustomNucleotidePopup> ui;
    bool m_updating_model = false;

    void
    onModelValuesChanged(const std::unordered_set<ModelKey>& keys) override;
    void onTextEdited();
    void updateFromModel();
};

} // namespace sketcher
} // namespace schrodinger
