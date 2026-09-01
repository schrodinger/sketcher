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
 * Clear all contents from the sketcher canvas.
 */
export async function clearSketcher(page) {
  await page.evaluate(() => Module.sketcher_clear());
}

/**
 * Import a structure string and report the outcome without throwing, so
 * callers can sweep many structures and collect failures. Returns
 * {ok, empty, error} where `ok` is false if the import threw, and `empty`
 * reflects whether anything was actually rendered afterward.
 */
export async function tryImport(page, text) {
  return page.evaluate((t) => {
    try {
      Module.sketcher_import_text(t);
      return { ok: true, empty: Module.sketcher_is_empty(), error: null };
    } catch (e) {
      // Emscripten throws C++ exceptions as a raw pointer (a number).
      // Decode it into the original what() message when possible.
      let error;
      if (typeof e === 'number' && typeof Module.getExceptionMessage === 'function') {
        try {
          error = Module.getExceptionMessage(e).join(': ');
        } catch {
          error = `uncaught C++ exception (ptr ${e})`;
        }
      } else {
        error = String(e && e.message ? e.message : e);
      }
      return { ok: false, empty: true, error };
    }
  }, text);
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

/**
 * Programmatically click a popup button by its Qt objectName.
 * Qt::Popup windows in WASM create a separate canvas whose event listeners
 * can't be targeted by Playwright, so we call into C++ directly.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} name - Qt objectName of the popup button
 * @throws if no button with the given name is found
 */
export async function clickPopupButton(page, name) {
  await page.evaluate((n) => Module._sketcher_click_button(n), name);
}
