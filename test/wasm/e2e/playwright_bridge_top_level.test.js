import { expect, test } from '@playwright/test';
import {
  activateAction,
  clickPopupTool,
  getExportedSmiles,
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
    await activateAction(page, 'Paste in Text...');

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
