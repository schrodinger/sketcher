/**
 * Navigate to the sketcher page and wait for the WASM module to be fully
 * loaded and the canvas to be visible.
 */
export async function waitForSketcherReady(page) {
  await page.goto('/wasm_shell.html');
  await page.waitForFunction(() => typeof window.Module !== 'undefined', {
    timeout: 20000,
  });
  // Wait for the canvas inside the shadow DOM to be attached and visible
  const canvas = page.locator('#screen canvas');
  await canvas.waitFor({ state: 'visible', timeout: 10000 });
}

/**
 * Focus the shadow DOM canvas so keyboard events reach Qt.
 */
export async function focusCanvas(page) {
  await page.locator('#screen canvas').focus();
}

/**
 * Return the {x, y} center of the canvas bounding box.
 * The center is well past the left toolbar (~90px wide), safely in the
 * drawing area.
 */
export async function getCanvasCenter(page) {
  const box = await page.locator('#screen canvas').boundingBox();
  return { x: Math.round(box.x + box.width / 2), y: Math.round(box.y + box.height / 2) };
}

/**
 * Return the {x, y, width, height} rect of any child widget by its Qt
 * objectName. Uses a single generic C++ function so no new WASM bindings
 * are needed when new widgets are introduced.
 * @param {import('@playwright/test').Page} page
 * @param {string} objectName - Qt objectName of the widget
 */
export async function getWidgetRect(page, objectName) {
  const rect = await page.evaluate(
    (name) => JSON.parse(Module._sketcher_get_widget_rect(name)),
    objectName,
  );
  if (!rect || rect.width === undefined) {
    throw new Error(`Widget "${objectName}" not found`);
  }
  return rect;
}

/**
 * Return the {x, y} center of the drawing area (excluding toolbar and top bar).
 * @param {import('@playwright/test').Page} page
 */
export async function getDrawingAreaCenter(page) {
  const rect = await getWidgetRect(page, 'view');
  return { x: Math.round(rect.x + rect.width / 2), y: Math.round(rect.y + rect.height / 2) };
}

/**
 * Return the current molecule as a SMILES string.
 */
export async function getExportedSmiles(page) {
  return page.evaluate(() => Module.sketcher_export_text(Module.Format.SMILES));
}

/**
 * Return the current molecule as a HELM string.
 */
export async function getExportedHelm(page) {
  return page.evaluate(() => Module.sketcher_export_text(Module.Format.HELM));
}

/**
 * Enable monomeric (peptide/nucleic acid) mode on the sketcher.
 */
export async function enableMonomericMode(page) {
  await page.evaluate(() => Module.sketcher_allow_monomeric(true));
}

/**
 * Select all items on the canvas via Cmd+A / Ctrl+A.
 */
export async function selectAll(page) {
  await focusCanvas(page);
  await page.keyboard.press('ControlOrMeta+a');
}

/**
 * Return whether the sketcher is currently empty.
 */
export async function isSketcherEmpty(page) {
  return page.evaluate(() => Module.sketcher_is_empty());
}

/**
 * Click a toolbar button by its Qt objectName (e.g. "c_btn").
 * Uses the generic getWidgetRect to find the button's position.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} name - Qt objectName of the button
 */
export async function clickWidget(page, name) {
  const rect = await getWidgetRect(page, name);
  await page.mouse.click(rect.x + rect.width / 2, rect.y + rect.height / 2);
}
