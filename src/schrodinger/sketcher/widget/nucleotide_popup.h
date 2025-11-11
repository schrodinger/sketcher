#pragma once
#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/modular_popup.h"

class QString;

namespace Ui
{
class NucleotidePopup;
}

namespace schrodinger
{
namespace sketcher
{

enum class NucleicAcidTool;
enum class ModelKey;
enum class StdNucleobase;

/**
 * Popup used to provide base choices for the nucleotide buttons
 */
class SKETCHER_API NucleotidePopup : public ModularPopup
{
  public:
    /**
     * @param tool The tool that this popup is setting the base for. Should be
     * either RNA_NUCLEOTIDE or DNA_NUCLEOTIDE.
     * @param ModelKey The key that the base should be stored in.  Should be
     * either RNA_NUCLEOBASE or DNA_NUCLEOBASE.
     * @param sugar The sugar to display on the buttons. Should be either "R" or
     * "dR".
     * @param u_or_t The name of the pyrimidine.  Should be either "U" or "T".
     * @param parent The Qt parent of this popup
     */
    NucleotidePopup(const NucleicAcidTool tool, const ModelKey model_key,
                    const QString& sugar, const QString& u_or_t,
                    QWidget* parent = nullptr);
    ~NucleotidePopup();

  protected:
    NucleicAcidTool m_tool;
    ModelKey m_model_key;
    QString m_sugar;
    QString m_u_or_t;

    std::unique_ptr<Ui::NucleotidePopup> ui;
    void generateButtonPackets() override;
    int getButtonIDToCheck() override;
};

} // namespace sketcher
} // namespace schrodinger
