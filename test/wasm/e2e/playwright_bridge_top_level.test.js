import { expect, test } from '@playwright/test';
import {
  clickMenuButtonRow,
  clickPopupTool,
  getExportedSmiles,
  isWidgetVisible,
  loadStructure,
  openContextMenu,
  selectAll,
  waitForSketcherReady,
  widgetState,
} from './e2e_helpers.js';

const SOURCE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';

test.describe('Playwright bridge top-level Qt surfaces', () => {
  test('reports text and style from the foreground Paste in Text dialog', async ({ page }) => {
    await waitForSketcherReady(page);
    await loadStructure(page, SOURCE);
    await clickMenuButtonRow(page, 'import_btn', 'Paste in Text...');
    // The dialog opens on a later event loop iteration, and until then
    // status_lbl would resolve to a hidden label of the same name.
    await expect.poll(() => isWidgetVisible(page, 'structure_text_edit')).toBe(true);

    const status = await widgetState(page, 'status_lbl');
    expect(status.text).toBe('Specified structure will <b>replace</b> Sketcher content');
    expect(status.styleSheet).toBe('color : #c87c00');
  });

  test('clicks a popup tool embedded in a context menu', async ({ page }) => {
    await waitForSketcherReady(page);
    await loadStructure(page, SOURCE);
    await selectAll(page);
    await openContextMenu(page, 'atom:0', 'Modify Atoms', 'Set Element');

    await expect(page.evaluate(() => Module._sketcher_get_popup_owner('pd_btn'))).resolves.toBe(
      'periodic_table_btn',
    );
    await clickPopupTool(page, 'pd_btn');
    expect(await getExportedSmiles(page)).toContain('[Pd]');
  });
});
