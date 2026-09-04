import { expect, test } from '@playwright/test';

const SOURCE = 'NC(N)=NC(=O)CC1=C(Cl)C=CC=C1Cl';

async function waitForSketcherReady(page) {
  await page.goto('/wasm_shell.html');
  await page.waitForFunction(() =>
    typeof window.Module?._sketcher_get_rect === 'function' &&
    typeof window.Module?.sketcher_import_text === 'function');
}

async function bridgeRect(page, selector) {
  return page.evaluate((value) => JSON.parse(Module._sketcher_get_rect(value)), selector);
}

async function clickBridgeRect(page, selector) {
  const rect = await bridgeRect(page, selector);
  await page.mouse.click(rect.x + rect.width / 2, rect.y + rect.height / 2);
}

test.describe('generic geometry bridge contract', () => {
  test('resolves widget, scene item, menu, state, and popup-owner selectors', async ({ page }) => {
    await waitForSketcherReady(page);

    const view = await bridgeRect(page, 'widget:view');
    expect(view.enabled).toBe(true);
    expect(view.width).toBeGreaterThan(100);

    const selectTool = await bridgeRect(page, 'state:select_tool_btn');
    expect(selectTool.visible).toBe(true);
    expect(typeof selectTool.enabled).toBe('boolean');

    await page.evaluate((text) => Module.sketcher_import_text(text), SOURCE);
    const atom = await bridgeRect(page, 'atom:0');
    const bond = await bridgeRect(page, 'bond:0');
    expect(atom.width).toBeGreaterThan(0);
    expect(bond.width).toBeGreaterThan(0);

    await page.mouse.click(bond.x + bond.width / 2, bond.y + bond.height / 2, {
      button: 'right',
    });
    const deleteAction = await bridgeRect(page, 'menu:Delete');
    expect(deleteAction.width).toBeGreaterThan(0);
    await page.keyboard.press('Escape');

    await expect(page.evaluate(() => Module._sketcher_get_popup_owner('double_btn')))
      .resolves.toBe('bond_order_btn');

    // State selectors deliberately resolve popup tools while their popup is
    // closed; that is how shortcut-driven tool selection is asserted.
    const doubleBond = await bridgeRect(page, 'state:double_btn');
    expect(doubleBond.visible).toBe(false);
    expect(doubleBond.width).toBeGreaterThan(0);
  });

  test('opens and selects a ToolButtonWithPopup child with generic widget geometry', async ({ page }) => {
    await waitForSketcherReady(page);

    // A normal first click selects the displayed parent tool. A second normal
    // click opens its popup, matching the user-facing sequence used by the
    // source-parity wrapper.
    await clickBridgeRect(page, 'widget:bond_order_btn');
    await clickBridgeRect(page, 'widget:bond_order_btn');
    const doubleBond = await bridgeRect(page, 'widget:double_btn');
    expect(doubleBond.width).toBeGreaterThan(0);
    const visibleState = await bridgeRect(page, 'state:double_btn');
    expect(visibleState.visible).toBe(true);
    await page.mouse.click(
      doubleBond.x + doubleBond.width / 2,
      doubleBond.y + doubleBond.height / 2,
    );

    const selected = await bridgeRect(page, 'state:double_btn');
    expect(selected.checked).toBe(true);
  });
});
