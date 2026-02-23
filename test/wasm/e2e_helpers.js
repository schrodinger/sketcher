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
 * Return the {x, y} center of the drawing area (excluding toolbar and top bar).
 * Uses the C++ side to get the actual View widget geometry.
 * @param {import('@playwright/test').Page} page
 */
export async function getDrawingAreaCenter(page) {
  const rect = await page.evaluate(() => JSON.parse(Module._sketcher_get_drawing_area_rect()));
  return { x: Math.round(rect.x + rect.width / 2), y: Math.round(rect.y + rect.height / 2) };
}

/**
 * Return the current molecule as a SMILES string.
 */
export async function getExportedSmiles(page) {
  return page.evaluate(() => Module.sketcher_export_text(Module.Format.SMILES));
}

/**
 * Return whether the sketcher is currently empty.
 */
export async function isSketcherEmpty(page) {
  return page.evaluate(() => Module.sketcher_is_empty());
}

/**
 * Execute an action and wait for Module.sketcher_changed_callback to fire.
 * The callback is invoked by the C++ side via a 100ms setTimeout after any
 * molecule change (see main.cpp sketcher_changed()).
 *
 * @param {import('@playwright/test').Page} page
 * @param {() => Promise<void>} action - The interaction to perform (click, key, etc.)
 * @param {number} timeout - Max time to wait for the callback (ms)
 */
export async function waitForMoleculeChange(page, action, timeout = 5000) {
  // Install a one-shot promise on the page that resolves when the callback fires
  await page.evaluate(() => {
    window.__moleculeChanged = new Promise((resolve) => {
      const prev = Module.sketcher_changed_callback;
      Module.sketcher_changed_callback = () => {
        Module.sketcher_changed_callback = prev;
        resolve();
      };
    });
  });

  await action();

  await page.evaluate(
    (ms) =>
      Promise.race([
        window.__moleculeChanged,
        new Promise((_, reject) =>
          setTimeout(() => reject(new Error('Timed out waiting for molecule change')), ms),
        ),
      ]),
    timeout,
  );
}

/**
 * Return the modifier key for shortcuts (Meta on macOS, Control elsewhere).
 * Qt WASM on macOS expects Cmd (Meta) for undo/redo, not Ctrl.
 */
export function modifierKey() {
  return process.platform === 'darwin' ? 'Meta' : 'Control';
}

/**
 * Query the C++ side for actual toolbar button positions via the WASM API.
 * Returns a map of objectName -> {x, y, text}.
 * @param {import('@playwright/test').Page} page
 */
export async function getToolbarButtonPositions(page) {
  return page.evaluate(() => JSON.parse(Module._sketcher_get_button_positions()));
}

/**
 * Click a toolbar button by its Qt objectName (e.g. "c_btn").
 * Positions are resolved dynamically via sketcher_get_button_positions().
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} name - Qt objectName of the button
 */
export async function clickToolbarButton(page, name) {
  const positions = await getToolbarButtonPositions(page);
  const btn = positions[name];
  if (!btn) {
    throw new Error(
      `Unknown toolbar button "${name}". Available: ${Object.keys(positions).join(', ')}`,
    );
  }
  await page.mouse.click(btn.x, btn.y);
}

/**
 * Wait for a render frame plus a short settling delay.
 */
export async function waitForRender(page, ms = 200) {
  await page.evaluate(
    (delay) => new Promise((resolve) => requestAnimationFrame(() => setTimeout(resolve, delay))),
    ms,
  );
}
