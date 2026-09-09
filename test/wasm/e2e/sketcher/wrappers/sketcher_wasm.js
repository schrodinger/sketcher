import { readFile } from 'node:fs/promises';
import path from 'node:path';

// The generic selector path has passed the full e2e suite. Keep the proven
// dedicated bridge available as an explicit rollback switch while generic
// geometry is the normal wrapper behavior.
export const USE_GENERIC_GEOMETRY_BRIDGE =
  process.env.SKETCHER_USE_LEGACY_GEOMETRY_BRIDGE !== '1';

/** Return a selector rectangle from the generic C++ geometry bridge. */
export async function genericRect(page, selector, { attempts = 150, interval = 100 } = {}) {
  for (let attempt = 0; attempt < attempts; attempt += 1) {
    const raw = await page.evaluate((value) => Module._sketcher_get_rect(value), selector);
    if (typeof raw === 'string' && raw.trim() !== '') {
      try {
        const rect = JSON.parse(raw);
        if (rect?.width !== undefined) return rect;
      } catch {
        // Qt may expose an incomplete string for one repaint frame.
      }
    }
    await page.waitForTimeout(interval);
  }
  throw new Error(`Sketcher generic geometry not found: ${selector}`);
}

/**
 * Browser-facing Sketcher wrapper used by the Squish-to-Playwright port.
 *
 * The Sketcher UI is Qt rendered into a canvas, so DOM locators cannot reach
 * individual controls.  The WASM test bridge exposes Qt object names through
 * `Module._sketcher_get_widget_rect` and `Module._sketcher_click_button`.
 */

/** @param {import('@playwright/test').Page} page */
export async function openSketcher(page) {
  // QFileDialog::saveFileContent in Qt/WASM ultimately triggers a synthetic
  // browser download.  In headed Linux runs Chromium may show a native Save
  // As chooser, which Playwright cannot inspect or dismiss.  The actual
  // Sketcher Download button is still clicked normally and the test bridge
  // observes its already-created bytes; suppress only this untestable final
  // browser handoff so it does not obscure a visual run.
  await page.addInitScript(() => {
    // Qt 6.7 uses the browser File System Access save picker on platforms
    // that expose it. Supply an in-memory writable handle: the preceding Qt
    // Download click remains real and the test bridge has already recorded
    // the exact bytes at this boundary.
    window.showSaveFilePicker = async () => {
      window.__sketcherPlaywrightSuppressedSavePickers =
        (window.__sketcherPlaywrightSuppressedSavePickers || 0) + 1;
      return {
        createWritable: async () => ({ close: async () => {}, write: async () => {} }),
        kind: 'file',
        name: 'sketcher-playwright-export',
      };
    };
    const nativeClick = HTMLAnchorElement.prototype.click;
    HTMLAnchorElement.prototype.click = function sketcherDownloadShim(...args) {
      if (this.hasAttribute('download')) {
        window.__sketcherPlaywrightSuppressedDownloads =
          (window.__sketcherPlaywrightSuppressedDownloads || 0) + 1;
        return;
      }
      return nativeClick.apply(this, args);
    };
  });
  await page.goto('/wasm_shell.html');
  await page.waitForFunction(() => typeof window.Module !== 'undefined', undefined, {
    timeout: 20000,
  });
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
  if (USE_GENERIC_GEOMETRY_BRIDGE) return genericRect(page, `widget:${objectName}`);
  // Always retrieve the live bridge geometry.  The cached startup geometry is
  // correct for ordinary child widgets, but modular-tool popup controls live
  // in separate top-level Qt/WASM windows.  Treating those as children gives
  // a plausible but wrong canvas point, so a real Playwright click lands on
  // the underlying toolbar instead of the visible popup item.
  // Qt/WASM can return an empty string for one animation frame while a popup
  // closes or the canvas is being repainted.  Wait for a complete bridge
  // response instead of surfacing that transient state as a JSON parse error.
  // Do not use `waitForFunction` with an async bridge call here. Its polling
  // scheduler can begin another Asyncify invocation before Qt/WASM completes
  // the first, which aborts with "multiple async operations in flight" in
  // unoptimized local artifacts. Sequential retries preserve one real bridge
  // operation at a time.
  for (let attempt = 0; attempt < 150; attempt += 1) {
    const raw = await page.evaluate((name) => Module._sketcher_get_widget_rect(name), objectName);
    if (typeof raw === 'string' && raw.trim() !== '') {
      try {
        const rect = JSON.parse(raw);
        if (rect?.width !== undefined) return rect;
      } catch {
        // Qt may expose an incomplete string for one repaint frame.
      }
    }
    await page.waitForTimeout(100);
  }
  throw new Error(`Sketcher widget not found: ${objectName}`);
}

/**
 * Activate a visible Qt button by stable objectName.  Qt renders the controls
 * into its canvas, so the input itself must be sent as a canvas click.
 */
export async function clickWidget(page, objectName) {
  const rect = await widgetRect(page, objectName);
  await mouseClick(page, rect.x + rect.width / 2, rect.y + rect.height / 2);
}

/** Squish-compatible Qt press event for non-browser-rendered QWidgetActions. */
export async function sendWidgetMousePress(page, objectName) {
  await page.evaluate((name) => Module._sketcher_send_mouse_press(name), objectName);
}

/** Test-only cleanup for Qt/WASM QWidgetAction popup bookkeeping. */
export async function closeActiveQtPopups(page) {
  await page.evaluate(() => Module._sketcher_close_active_popups());
}

/**
 * Clear the Playwright-only record made at Qt's QFileDialog save boundary.
 * Qt/WASM does not emit Chromium's download event, so the bridge observes the
 * bytes after the real visible Download click has produced them.
 */
export async function beginBrowserDownloadCapture(page) {
  await page.evaluate(() => {
    window.__sketcherPlaywrightFileExports = [];
  });
}

/** Wait for and return the payload produced by a real Qt/WASM Download click. */
export async function capturedBrowserDownload(page) {
  const result = await page.waitForFunction(
    () => window.__sketcherPlaywrightFileExports?.[0] || null,
    undefined,
    { timeout: 10000 },
  );
  return result.jsonValue();
}

/** Return visible Qt widget state exposed by the Playwright test bridge. */
export async function widgetState(page, objectName) {
  for (let attempt = 0; attempt < 150; attempt += 1) {
    const raw = await page.evaluate((name) => Module._sketcher_get_widget_state(name), objectName);
    if (typeof raw === 'string' && raw.trim() !== '') {
      try {
        return JSON.parse(raw);
      } catch {
        // Qt may expose an incomplete string for one repaint frame.
      }
    }
    await page.waitForTimeout(100);
  }
  throw new Error(`Sketcher widget state not found: ${objectName}`);
}

/** Click a visible text control and replace its value through keyboard input. */
export async function setWidgetText(page, objectName, text) {
  await clickWidget(page, objectName);
  await page.keyboard.press('ControlOrMeta+a');
  await page.keyboard.type(String(text), { delay: 10 });
}

/** Return a visible Qt menu action's canvas rectangle by objectName or text. */
export async function menuActionRect(page, objectNameOrText) {
  if (USE_GENERIC_GEOMETRY_BRIDGE) {
    // Preserve the legacy menu retry budget. A context-menu wrapper needs to
    // abandon a transiently absent nested action quickly, reopen the menu,
    // and repeat the same human hover path; a widget-sized 15-second retry
    // here consumes its entire test timeout before that recovery can occur.
    return genericRect(page, `menu:${objectNameOrText}`, { attempts: 50, interval: 25 });
  }
  for (let attempt = 0; attempt < 50; attempt += 1) {
    const rect = await page.evaluate(async (name) => {
      const raw = await Module._sketcher_get_menu_action_rect(name);
      try {
        const value = JSON.parse(raw);
        return value?.width === undefined ? null : value;
      } catch {
        return null;
      }
    }, objectNameOrText);
    if (rect) return rect;
    await page.waitForTimeout(25);
  }
  throw new Error(`Sketcher menu action not found: ${objectNameOrText}`);
}

/** Click a visible Qt menu action using the browser's real mouse input. */
export async function clickMenuAction(page, objectNameOrText) {
  const rect = await menuActionRect(page, objectNameOrText);
  await mouseClick(page, rect.x + rect.width / 2, rect.y + rect.height / 2);
}

/**
 * Click a row in a Qt/WASM popup through its browser canvas.  Qt gives every
 * popup its own `qt-window-N` canvas; this avoids unsafe QMenu introspection.
 */
export async function clickPopupRow(page, popupIndex, rowCenter) {
  const popup = await popupCanvasGeometry(page, popupIndex);
  await mouseClick(page, popup.x + popup.width / 2, popup.y + rowCenter);
}

/** Return a Qt popup canvas's CSS offset and backing-canvas dimensions. */
export async function popupCanvasGeometry(page, popupIndex) {
  const popups = await page
    .locator('[id^="qt-window-"]')
    .evaluateAll((elements) => {
      return elements
        .slice(1)
        .map((element) => {
          const canvas = element.querySelector('canvas');
          const style = getComputedStyle(element);
          // Qt/WASM keeps a DOM canvas for a closed QComboBox popup.  It is
          // no longer visible to a person, but was still being counted as
          // popup zero and could receive a later menu click.  Restrict the
          // geometry map to canvases that are actually visible in the page.
          if (style.display === 'none' || style.visibility === 'hidden') return null;
          const bounds = element.getBoundingClientRect();
          if (bounds.width === 0 || bounds.height === 0) return null;
          const x = Number.parseFloat(style.left);
          const y = Number.parseFloat(style.top);
          return Number.isFinite(x) && Number.isFinite(y) && canvas?.width && canvas?.height
            ? { height: canvas.height, width: canvas.width, x, y }
            : null;
        })
        .filter(Boolean);
    })
    .catch(() => []);
  // A persistent Qt dialog can retain an older hidden canvas in the DOM
  // after closing.  A newly opened root menu is appended last, so callers
  // that need the currently raised root menu can request `latest`.
  const popup = popupIndex === 'latest' ? popups.at(-1) : popups[popupIndex];
  if (!popup) throw new Error(`Qt popup ${popupIndex} did not expose canvas geometry`);
  return popup;
}

/** Open a cascading Qt submenu by moving the real pointer over its parent row. */
export async function hoverPopupRow(page, popupIndex, rowCenter) {
  const popup = await popupCanvasGeometry(page, popupIndex);
  await showMouseMarker(page, popup.x + popup.width / 2, popup.y + rowCenter);
  await page.mouse.move(popup.x + popup.width / 2, popup.y + rowCenter, { steps: 4 });
  await page.waitForTimeout(250);
}

/** Open a cascading Qt submenu by moving the real pointer over its parent. */
export async function hoverMenuAction(page, objectNameOrText) {
  const rect = await menuActionRect(page, objectNameOrText);
  await showMouseMarker(page, rect.x + rect.width / 2, rect.y + rect.height / 2);
  await page.mouse.move(rect.x + rect.width / 2, rect.y + rect.height / 2, { steps: 4 });
  // Nested Qt/WASM menus are attached on a later browser frame.  Give the
  // ordinary human hover enough time to expose its child popup before the
  // next source-menu path step resolves its geometry.
  await page.waitForTimeout(250);
}

export async function clipboardText(page) {
  return page.evaluate(async () => {
    if (!navigator.clipboard?.readText) {
      throw new Error('The browser Clipboard API is unavailable in this context.');
    }
    return navigator.clipboard.readText();
  });
}

/**
 * Provide the result of the browser-native File System Access picker.
 *
 * Qt/WASM 6.7 invokes `showOpenFilePicker()` rather than an HTML file input,
 * so Playwright cannot use its `filechooser` event.  This replaces only the
 * native-picker result; the Sketcher Import menu is still opened and selected
 * through ordinary pointer input.
 */
export async function setFilePickerResult(page, filePath) {
  const contents = await readFile(filePath);
  await page.evaluate(
    ({ base64, name }) => {
      const bytes = Uint8Array.from(atob(base64), (character) => character.charCodeAt(0));
      const file = new File([bytes], name);
      window.showOpenFilePicker = async () => [{ getFile: async () => file, kind: 'file', name }];
    },
    { base64: contents.toString('base64'), name: path.basename(filePath) },
  );
}

export async function setClipboardText(page, text) {
  await page.evaluate((value) => navigator.clipboard.writeText(value), text);
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
