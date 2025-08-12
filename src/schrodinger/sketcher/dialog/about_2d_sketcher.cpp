#include "schrodinger/sketcher/dialog/about_2d_sketcher.h"

#include <ctime>
#include <iomanip>
#include <sstream>

#include <QDesktopServices>
#include <QUrl> // required with Qt5
#include <rdkit/RDGeneral/versions.h>

#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_about_2d_sketcher.h"
#include "schrodinger/sketcher/version.h"

namespace schrodinger
{
namespace sketcher
{

About2DSketcher::About2DSketcher(QWidget* parent) : ModalDialog(parent)
{
    m_ui.reset(new Ui::About2DSketcher());
    setupDialogUI(*m_ui);

    auto schrodinger_version = std::string("Release ") +
                               std::string(SKETCHER_RELEASE) + " (build " +
                               std::to_string(SKETCHER_BUILD) + ")";
    m_ui->schrodinger_version_lbl->setText(schrodinger_version.c_str());

    auto rdkit_version = std::string("Version ") + RDKit::rdkitVersion;
    m_ui->rdkit_version_lbl->setText(rdkit_version.c_str());

    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::stringstream ss;
    ss << std::put_time(&tm, "%Y");
    auto copyright_text = std::string("© ") + ss.str() + " Schrödinger, Inc.";
    m_ui->copyright_lbl->setText(copyright_text.c_str());

    m_ui->eula_btn->setStyleSheet(TEXT_LINK_STYLE);
    connect(m_ui->eula_btn, &QToolButton::clicked, this, []() {
        QUrl url("https://www.schrodinger.com/salesagreements");
        QDesktopServices::openUrl(url);
    });
}

About2DSketcher::~About2DSketcher() = default;

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/about_2d_sketcher.moc"
