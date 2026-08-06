/**
 * Browser-facing Sketcher wrapper used by the Squish-to-Playwright port.
 *
 * The Sketcher UI is Qt rendered into a canvas, so DOM locators cannot reach
 * individual controls.  The WASM test bridge exposes Qt object names through
 * `Module._sketcher_get_widget_rect` and `Module._sketcher_click_button`.
 */

/** @param {import('@playwright/test').Page} page */
export async function openSketcher(page) {
  await page.goto('/wasm_shell.html');
  await page.waitForFunction(() => typeof window.Module !== 'undefined', { timeout: 20000 });
  await page.locator('#screen canvas').waitFor({ state: 'visible', timeout: 10000 });
}

/** @param {import('@playwright/test').Page} page */
export async function focusCanvas(page) {
  await page.locator('#screen canvas').focus();
}

/** @param {import('@playwright/test').Page} page */
export async function drawingAreaCenter(page) {
  const rect = await widgetRect(page, 'view');
  return { x: Math.round(rect.x + rect.width / 2), y: Math.round(rect.y + rect.height / 2) };
}

/**
 * Return a Qt widget's viewport rectangle by its stable objectName.
 * @param {import('@playwright/test').Page} page
 * @param {string} objectName
 */
export async function widgetRect(page, objectName) {
  const rect = await page.evaluate(
    (name) => JSON.parse(Module._sketcher_get_widget_rect(name)),
    objectName,
  );
  if (!rect || rect.width === undefined) {
    throw new Error(`Sketcher widget not found: ${objectName}`);
  }
  return rect;
}

/**
 * Activate a visible Qt button by stable objectName.  Qt renders the controls
 * into its canvas, so the input itself must be sent as a canvas click.
 */
export async function clickWidget(page, objectName) {
  const rect = await widgetRect(page, objectName);
  await page.mouse.click(rect.x + rect.width / 2, rect.y + rect.height / 2);
}

/** Return visible Qt widget state exposed by the Playwright test bridge. */
export async function widgetState(page, objectName) {
  return page.evaluate((name) => JSON.parse(Module._sketcher_get_widget_state(name)), objectName);
}

/** Set text, a combo-box choice, or a spin-box value through the test bridge. */
export async function setWidgetText(page, objectName, text) {
  await page.evaluate(
    ([name, value]) => Module._sketcher_set_widget_text(name, value),
    [objectName, String(text)],
  );
}

/** Trigger a currently visible Qt menu action by objectName or displayed text. */
export async function activateMenuAction(page, objectNameOrText) {
  await page.evaluate((name) => Module._sketcher_activate_menu_action(name), objectNameOrText);
}

export async function clipboardText(page) {
  return page.evaluate(() => Module._sketcher_clipboard_text());
}

export async function setClipboardText(page, text) {
  await page.evaluate((value) => Module._sketcher_set_clipboard_text(value), text);
}

/** @param {import('@playwright/test').Page} page */
export async function clearSketcher(page) {
  await page.evaluate(() => Module.sketcher_clear());
}

/** @param {import('@playwright/test').Page} page @param {string} text */
export async function importText(page, text) {
  await page.evaluate((value) => Module.sketcher_import_text(value), text);
}

/** @param {import('@playwright/test').Page} page @param {string} format */
export async function exportText(page, format = 'SMILES') {
  return page.evaluate((name) => Module.sketcher_export_text(Module.Format[name]), format);
}

/** @param {import('@playwright/test').Page} page */
export async function isEmpty(page) {
  return page.evaluate(() => Module.sketcher_is_empty());
}

/** @param {import('@playwright/test').Page} page @param {string} element */
export async function drawElement(page, element) {
  const center = await drawingAreaCenter(page);
  await focusCanvas(page);
  await page.keyboard.press(element);
  await page.mouse.click(center.x, center.y);
}

/** @param {import('@playwright/test').Page} page */
export async function drawBond(page) {
  const center = await drawingAreaCenter(page);
  await page.mouse.move(center.x, center.y);
  await page.mouse.down();
  await page.mouse.move(center.x + 100, center.y, { steps: 10 });
  await page.mouse.up();
}
