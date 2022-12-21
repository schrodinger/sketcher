#include "schrodinger/sketcher/widget/sketcher_bottom_bar.h"
#include "schrodinger/sketcher/ui/ui_sketcher_bottom_bar.h"
#include <QPushButton>
#include <QInputDialog>
#include <QLineEdit>
#include <QString>

namespace schrodinger
{
namespace sketcher
{

SketcherBottomBar::SketcherBottomBar(QWidget* parent) : QWidget(parent)
{
    ui.reset(new Ui::SketcherBottomBarForm());
    ui->setupUi(this);

    connect(ui->reset_btn, &QPushButton::pressed, this,
            &SketcherBottomBar::resetRequested);
    connect(ui->reload_btn, &QPushButton::pressed, this,
            &SketcherBottomBar::reloadRequested);
    connect(ui->save_as_new_btn, &QPushButton::pressed, this,
            &SketcherBottomBar::onSaveAsNewClicked);
    connect(ui->save_changes_btn, &QPushButton::pressed, this,
            &SketcherBottomBar::saveChangesRequested);
}

SketcherBottomBar::~SketcherBottomBar() = default;

void SketcherBottomBar::onSaveAsNewClicked()
{
    bool dialog_accepted = false;
    QString entry_title = QInputDialog::getText(
        this, tr("Save As New Entry"), tr("Input entry title:"),
        QLineEdit::Normal, tr(""), &dialog_accepted);
    if (dialog_accepted) {
        emit saveAsNewEntryRequested(entry_title.toStdString());
    }
}

}; // namespace sketcher
}; // namespace schrodinger

#include "schrodinger/sketcher/widget/sketcher_bottom_bar.moc"
