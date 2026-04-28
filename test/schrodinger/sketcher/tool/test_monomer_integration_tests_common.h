#include <memory>
#include <string>

#include <QApplication>
#include <QEvent>
#include <QGraphicsItem>
#include <QGraphicsSceneMouseEvent>
#include <QPointF>
#include <QString>
#include <QVariant>

#include <boost/test/unit_test.hpp>

#include <rdkit/GraphMol/ROMol.h>

#include "../test_common.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/unbound_monomeric_attachment_point_item.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * Process all Qt events, includeing DeferredDelete events.
 */
void process_qt_events()
{
    // call processEvents multiple times in case an any current events put new
    // events on the queue (e.g. starting a timer with a timeout of 0)
    for (int i = 0; i < 3; ++i) {
        QApplication::processEvents();
        // processEvents will never process DeferredDelete events, but
        // sendPostedEvents will if we explicitly pass their event type.
        // (Despite what Qt's documentation claims, passing 0 as the event type
        // processes everything *other than* DeferredDeletes.)
        QApplication::sendPostedEvents(nullptr, QEvent::DeferredDelete);
    }
}

/**
 * Set both the scene and screen pos for an event
 */
void set_event_pos(QGraphicsSceneMouseEvent& event, const QPointF& pos)
{
    event.setScenePos(pos);
    event.setScreenPos(pos.toPoint());
}

/**
 * Fixture that provides a clean Scene with MolModel and SketcherModel for
 * testing monomer drawing tools.
 */
struct MonomerToolTestFixture {
    std::shared_ptr<TestScene> m_scene;
    MolModel* m_mol_model;
    SketcherModel* m_sketcher_model;

    MonomerToolTestFixture()
    {
        m_scene = TestScene::getScene();
        m_mol_model = m_scene->m_mol_model;
        m_sketcher_model = m_scene->m_sketcher_model;
        m_sketcher_model->setValue(
            ModelKey::INTERFACE_TYPE,
            static_cast<int>(InterfaceType::ATOMISTIC_OR_MONOMERIC));
        process_qt_events();
    }

    void setAminoAcidTool(AminoAcidTool tool)
    {
        m_sketcher_model->setValues(
            {{ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::MONOMER)},
             {ModelKey::TOOL_SET, QVariant::fromValue(ToolSet::MONOMERIC)},
             {ModelKey::MONOMER_TOOL_TYPE,
              QVariant::fromValue(MonomerToolType::AMINO_ACID)},
             {ModelKey::AMINO_ACID_TOOL, QVariant::fromValue(tool)},
             {ModelKey::AMINO_ACID_SYMBOL, QString("")}});
        process_qt_events();
    }

    void setNucleicAcidTool(NucleicAcidTool tool)
    {
        m_sketcher_model->setValues(
            {{ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::MONOMER)},
             {ModelKey::TOOL_SET, QVariant::fromValue(ToolSet::MONOMERIC)},
             {ModelKey::MONOMER_TOOL_TYPE,
              QVariant::fromValue(MonomerToolType::NUCLEIC_ACID)},
             {ModelKey::NUCLEIC_ACID_TOOL, QVariant::fromValue(tool)}});
        process_qt_events();
    }

    void setRNANucleotideTool(StdNucleobase base)
    {
        m_sketcher_model->setValues(
            {{ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::MONOMER)},
             {ModelKey::TOOL_SET, QVariant::fromValue(ToolSet::MONOMERIC)},
             {ModelKey::MONOMER_TOOL_TYPE,
              QVariant::fromValue(MonomerToolType::NUCLEIC_ACID)},
             {ModelKey::NUCLEIC_ACID_TOOL,
              QVariant::fromValue(NucleicAcidTool::RNA_NUCLEOTIDE)},
             {ModelKey::RNA_NUCLEOBASE, QVariant::fromValue(base)}});
        process_qt_events();
    }

    void setDNANucleotideTool(StdNucleobase base)
    {
        m_sketcher_model->setValues(
            {{ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::MONOMER)},
             {ModelKey::TOOL_SET, QVariant::fromValue(ToolSet::MONOMERIC)},
             {ModelKey::MONOMER_TOOL_TYPE,
              QVariant::fromValue(MonomerToolType::NUCLEIC_ACID)},
             {ModelKey::NUCLEIC_ACID_TOOL,
              QVariant::fromValue(NucleicAcidTool::DNA_NUCLEOTIDE)},
             {ModelKey::DNA_NUCLEOBASE, QVariant::fromValue(base)}});
        process_qt_events();
    }

    void setCustomNucleotideTool(const QString& sugar, const QString& base,
                                 const QString& phosphate)
    {
        m_sketcher_model->setValues(
            {{ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::MONOMER)},
             {ModelKey::TOOL_SET, QVariant::fromValue(ToolSet::MONOMERIC)},
             {ModelKey::MONOMER_TOOL_TYPE,
              QVariant::fromValue(MonomerToolType::NUCLEIC_ACID)},
             {ModelKey::NUCLEIC_ACID_TOOL,
              QVariant::fromValue(NucleicAcidTool::CUSTOM_NUCLEOTIDE)},
             {ModelKey::CUSTOM_NUCLEOTIDE,
              QVariant::fromValue<std::tuple<QString, QString, QString>>(
                  {sugar, base, phosphate})}});
        process_qt_events();
    }

    void importMolText(const std::string& text)
    {
        import_mol_text(m_mol_model, text);
        process_qt_events();
    }

    QPointF getMonomerPos(unsigned int monomer_idx)
    {
        auto mol = m_mol_model->getMol();
        BOOST_REQUIRE(mol);
        BOOST_REQUIRE(monomer_idx < mol->getNumAtoms());
        return to_scene_xy(mol->getConformer().getAtomPos(monomer_idx));
    }

    QPointF getAttachmentPointPos(unsigned int monomer_idx,
                                  const std::string& ap_display_name)
    {
        auto monomer_pos = getMonomerPos(monomer_idx);
        auto* monomer_item = m_scene->getTopInteractiveItemAt(
            monomer_pos, InteractiveItemFlag::MONOMER);
        BOOST_REQUIRE(monomer_item != nullptr);

        // Search through child items to find the attachment point
        for (auto* child : monomer_item->childItems()) {
            auto* ap_item =
                qgraphicsitem_cast<UnboundMonomericAttachmentPointItem*>(child);
            if (ap_item) {
                auto ap = ap_item->getAttachmentPoint();
                if (ap.display_name == ap_display_name) {
                    auto offset_vec = direction_to_qt_vector(ap.direction);
                    auto offset_dist = UNBOUND_AP_LINE_LENGTH - 1 +
                                       STANDARD_AA_BORDER_WIDTH / 2;
                    return monomer_pos + offset_dist * offset_vec;
                }
            }
        }
        throw std::runtime_error("Attachment point " + ap_display_name +
                                 " not found");
    }

    void verifyHELM(const std::string& expected)
    {
        auto actual = get_mol_text(m_mol_model, rdkit_extensions::Format::HELM);
        BOOST_TEST(actual == expected);
    }

    void mouseMove(const QPointF& pos,
                   const Qt::MouseButtons btns = Qt::NoButton)
    {
        QGraphicsSceneMouseEvent event(QEvent::GraphicsSceneMouseMove);
        set_event_pos(event, pos);
        event.setButton(Qt::NoButton);
        event.setButtons(btns);
        m_scene->mouseMoveEvent(&event);
        process_qt_events();
    }

    void mousePress(const QPointF& pos)
    {
        QGraphicsSceneMouseEvent event(QEvent::GraphicsSceneMousePress);
        set_event_pos(event, pos);
        event.setButton(Qt::LeftButton);
        event.setButtons(Qt::LeftButton);
        m_scene->mousePressEvent(&event);
        process_qt_events();
    }

    void mouseRelease(const QPointF& pos)
    {
        QGraphicsSceneMouseEvent event(QEvent::GraphicsSceneMouseRelease);
        set_event_pos(event, pos);
        event.setButton(Qt::LeftButton);
        event.setButtons(Qt::NoButton);
        m_scene->mouseReleaseEvent(&event);
        process_qt_events();
    }

    void mouseClick(const QPointF& pos)
    {
        mouseMove(pos);
        mousePress(pos);
        mouseRelease(pos);
    }

    void mouseDrag(const QPointF& start, const QPointF& end)
    {
        mouseMove(start);
        mousePress(start);
        mouseMove(end, Qt::LeftButton);
        mouseRelease(end);
    }

    void confirmIsEmpty()
    {
        auto mol = m_mol_model->getMol();
        BOOST_TEST(mol->getNumAtoms() == 0);
    }
};

} // namespace sketcher
} // namespace schrodinger
