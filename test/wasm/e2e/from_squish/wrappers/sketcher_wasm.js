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

/** Show a test-only cursor marker when PLAYWRIGHT_SHOW_MOUSE=1. */
async function showMouseMarker(page, x, y) {
  if (process.env.PLAYWRIGHT_SHOW_MOUSE !== '1') return;
  await page.evaluate(
    ({ left, top }) => {
      let marker = document.getElementById('playwright-mouse-marker');
      if (!marker) {
        marker = document.createElement('div');
        marker.id = 'playwright-mouse-marker';
        Object.assign(marker.style, {
          background: 'rgba(255, 45, 45, 0.28)',
          border: '2px solid #ff2d2d',
          borderRadius: '50%',
          boxSizing: 'border-box',
          height: '18px',
          left: '0',
          pointerEvents: 'none',
          position: 'fixed',
          top: '0',
          transform: 'translate(-50%, -50%)',
          width: '18px',
          zIndex: '2147483647',
        });
        document.body.append(marker);
      }
      marker.style.display = 'block';
      marker.style.left = `${left}px`;
      marker.style.top = `${top}px`;
    },
    { left: x, top: y },
  );
}

/** Hide the optional marker before a canvas visual checkpoint. */
export async function hideMouseMarker(page) {
  if (process.env.PLAYWRIGHT_SHOW_MOUSE !== '1') return;
  await page.evaluate(() => {
    const marker = document.getElementById('playwright-mouse-marker');
    if (marker) marker.style.display = 'none';
  });
}

/** Move the pointer before pressing and releasing a mouse button. */
export async function mouseClick(page, x, y, { button = 'left', modifiers = [] } = {}) {
  for (const modifier of modifiers) {
    await page.keyboard.down(modifier);
  }
  try {
    await showMouseMarker(page, x, y);
    await page.mouse.move(x, y, { steps: 4 });
    await page.mouse.down({ button });
    await page.waitForTimeout(25);
    await page.mouse.up({ button });
  } finally {
    for (const modifier of [...modifiers].reverse()) {
      await page.keyboard.up(modifier);
    }
  }
}

/** Drag with ordinary pointer events, matching Squish's mouseDrag model. */
export async function mouseDrag(
  page,
  start,
  end,
  { button = 'left', modifiers = [], steps = 12 } = {},
) {
  for (const modifier of modifiers) {
    await page.keyboard.down(modifier);
  }
  try {
    await showMouseMarker(page, start.x, start.y);
    await page.mouse.move(start.x, start.y, { steps: 4 });
    await page.mouse.down({ button });
    for (let step = 1; step <= steps; step += 1) {
      const progress = step / steps;
      const x = start.x + (end.x - start.x) * progress;
      const y = start.y + (end.y - start.y) * progress;
      await showMouseMarker(page, x, y);
      await page.mouse.move(x, y);
    }
    await page.mouse.up({ button });
  } finally {
    for (const modifier of [...modifiers].reverse()) {
      await page.keyboard.up(modifier);
    }
  }
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
  const cached = await page.evaluate(
    (name) => window.__sketcherPlaywrightWidgetRects?.[name],
    objectName,
  );
  if (cached) {
    return cached;
  }
  const rect = await page.evaluate(
    async (name) => JSON.parse(await Module._sketcher_get_widget_rect(name)),
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
  await mouseClick(page, rect.x + rect.width / 2, rect.y + rect.height / 2);
}

/** Return visible Qt widget state exposed by the Playwright test bridge. */
export async function widgetState(page, objectName) {
  return page.evaluate((name) => JSON.parse(Module._sketcher_get_widget_state(name)), objectName);
}

/** Click a visible text control and replace its value through keyboard input. */
export async function setWidgetText(page, objectName, text) {
  await clickWidget(page, objectName);
  await page.keyboard.press('ControlOrMeta+a');
  await page.keyboard.type(String(text), { delay: 10 });
}

/** Return a visible Qt menu action's canvas rectangle by objectName or text. */
export async function menuActionRect(page, objectNameOrText) {
  let rect;
  try {
    rect = await page.waitForFunction(
      async (name) => {
        const value = JSON.parse(await Module._sketcher_get_menu_action_rect(name));
        return value?.width === undefined ? null : value;
      },
      objectNameOrText,
      { timeout: 5000 },
    );
  } catch (error) {
    const available = await page.evaluate(() =>
      Object.keys(window.__sketcherPlaywrightMenuRects || {}),
    );
    throw new Error(
      `Sketcher menu action not found: ${objectNameOrText}; visible actions: ${available.join(', ')}`,
      { cause: error },
    );
  }
  const value = await rect.jsonValue();
  if (!value || value.width === undefined) {
    throw new Error(`Sketcher menu action not found: ${objectNameOrText}`);
  }
  return value;
}

/** Click a visible Qt menu action using the browser's real mouse input. */
export async function clickMenuAction(page, objectNameOrText) {
  const rect = await menuActionRect(page, objectNameOrText);
  await mouseClick(page, rect.x + rect.width / 2, rect.y + rect.height / 2);
}

export async function clipboardText(page) {
  return page.evaluate(() => Module._sketcher_clipboard_text());
}

export async function setClipboardText(page, text) {
  await page.evaluate((value) => Module._sketcher_set_clipboard_text(value), text);
}

/** @param {import('@playwright/test').Page} page */
export async function clearSketcher(page) {
  await clickWidget(page, 'clear_btn');
}

/** @param {import('@playwright/test').Page} page @param {string} text */
export async function importText(page, text) {
  await clickWidget(page, 'import_btn');
  await clickMenuAction(page, 'Paste in Text...');
  await setWidgetText(page, 'structure_text_edit', text);

  // The generated Qt button box has no stable child object name for its OK
  // button. Its right half is the standard OK button in the standalone UI.
  const buttonBox = await widgetRect(page, 'buttonBox');
  await mouseClick(page, buttonBox.x + buttonBox.width * 0.75, buttonBox.y + buttonBox.height / 2);
}

/**
 * Load a structure as fixture setup, bypassing the import dialog.
 *
 * This is intentionally not an assertion of the user-facing Import flow.
 * Tests that cover Import must call importText() above; all other suites can
 * start from deterministic application state without introducing dialog
 * behavior as an unrelated dependency.
 */
export async function loadStructureForTest(page, text) {
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
  await mouseDrag(page, center, { x: center.x + 100, y: center.y });
}
