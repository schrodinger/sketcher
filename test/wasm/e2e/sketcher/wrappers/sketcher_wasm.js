import { readFile } from 'node:fs/promises';
import path from 'node:path';
import {
  clickWidget as clickGenericWidget,
  requireRect as requireBridgeRect,
  widgetState as bridgeWidgetState,
} from '../../e2e_helpers.js';

/** Return a selector rectangle from the shared generic C++ geometry bridge. */
export async function genericRect(page, selector) {
  return requireBridgeRect(page, selector);
}

/**
 * Browser-facing Sketcher wrapper used by the Squish-to-Playwright port.
 *
 * The Sketcher UI is Qt rendered into a canvas, so DOM locators cannot reach
 * individual controls. The shared WASM test bridge resolves stable Qt object
 * names and scene-item selectors to browser coordinates.
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
    window.showSaveFilePicker = async (options = {}) => {
      window.__sketcherPlaywrightSuppressedSavePickers =
        (window.__sketcherPlaywrightSuppressedSavePickers || 0) + 1;
      const filename = options.suggestedName || 'sketcher-playwright-export';
      return {
        createWritable: async () => {
          const chunks = [];
          return {
            write: async (data) => {
              const payload = data?.type === 'write' ? data.data : data;
              chunks.push(payload instanceof Blob ? payload : new Blob([payload]));
            },
            close: async () => {
              const blob = new Blob(chunks);
              const bytes = new Uint8Array(await blob.arrayBuffer());
              let binary = '';
              for (const byte of bytes) binary += String.fromCharCode(byte);
              window.__sketcherPlaywrightFileExportCallbackCount =
                (window.__sketcherPlaywrightFileExportCallbackCount || 0) + 1;
              if (!window.__sketcherPlaywrightFileExports) {
                window.__sketcherPlaywrightFileExports = [];
              }
              window.__sketcherPlaywrightFileExports.push({
                filename,
                contentBase64: btoa(binary),
              });
            },
          };
        },
        kind: 'file',
        name: filename,
      };
    };
    const nativeClick = HTMLAnchorElement.prototype.click;
    HTMLAnchorElement.prototype.click = function sketcherDownloadShim(...args) {
      if (this.hasAttribute('download')) {
        window.__sketcherPlaywrightSuppressedDownloads =
          (window.__sketcherPlaywrightSuppressedDownloads || 0) + 1;
        const filename = this.download || 'sketcher-playwright-export';
        fetch(this.href)
          .then((response) => response.arrayBuffer())
          .then((buffer) => {
            const bytes = new Uint8Array(buffer);
            let binary = '';
            for (const byte of bytes) binary += String.fromCharCode(byte);
            window.__sketcherPlaywrightFileExportCallbackCount =
              (window.__sketcherPlaywrightFileExportCallbackCount || 0) + 1;
            if (!window.__sketcherPlaywrightFileExports) {
              window.__sketcherPlaywrightFileExports = [];
            }
            window.__sketcherPlaywrightFileExports.push({
              filename,
              contentBase64: btoa(binary),
            });
          });
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
  let lastError;
  for (let attempt = 0; attempt < 100; attempt += 1) {
    try {
      return await genericRect(page, `widget:${objectName}`);
    } catch (error) {
      lastError = error;
      await page.waitForTimeout(25);
    }
  }
  throw lastError;
}

/**
 * Activate a visible Qt button by stable objectName.  Qt renders the controls
 * into its canvas, so the input itself must be sent as a canvas click.
 */
export async function clickWidget(page, objectName) {
  await clickGenericWidget(page, objectName);
}

/** Squish-compatible Qt press event for non-browser-rendered QWidgetActions. */
export async function sendWidgetMousePress(page, objectName) {
  const rect = await widgetRect(page, objectName);
  await showMouseMarker(page, rect.x + rect.width / 2, rect.y + rect.height / 2);
  await page.mouse.move(rect.x + rect.width / 2, rect.y + rect.height / 2, { steps: 4 });
  await page.mouse.down();
}

/** Release a press-only gesture and dismiss any remaining Qt popup by clicking outside it. */
export async function closeActiveQtPopups(page) {
  await page.mouse.up();
  // After some Qt/WASM menu actions the popup canvases remain active but do
  // not accept Escape. A normal click outside the menu dismisses the complete
  // cascading stack. Qt consumes that click rather than forwarding it to the
  // canvas underneath, exactly as it does for a person.
  const view = await widgetRect(page, 'view');
  // Use the inert top-bar background rather than the drawing surface. Some
  // Qt/WASM menu stacks let the outside click reach the view, which clears an
  // otherwise valid selection and changes the next context menu's contents.
  await mouseClick(page, view.x + view.width / 2, Math.max(2, view.y - 10));
}

/**
 * Clear the Playwright-only record made at Qt's QFileDialog save boundary.
 * Qt/WASM does not emit Chromium's download event, so the bridge observes the
 * bytes after the real visible Download click has produced them.
 */
export async function beginBrowserDownloadCapture(page) {
  await page.evaluate(() => {
    window.__sketcherPlaywrightFileExports = [];
    window.__sketcherPlaywrightFileExportCallbackCount = 0;
    window.__sketcherPlaywrightFileExportAttemptCount = 0;
  });
}

/** Wait for and return the payload produced by a real Qt/WASM Download click. */
export async function capturedBrowserDownload(page) {
  let result;
  try {
    result = await page.waitForFunction(
      () => window.__sketcherPlaywrightFileExports?.[0] || null,
      undefined,
      { timeout: 10000 },
    );
  } catch (error) {
    const callbackCount = await page.evaluate(
      () => window.__sketcherPlaywrightFileExportCallbackCount || 0,
    );
    const attemptCount = await page.evaluate(
      () => window.__sketcherPlaywrightFileExportAttemptCount || 0,
    );
    throw new Error(
      `Timed out waiting for browser export capture (attempts: ${attemptCount}, bytes recorded: ${callbackCount})`,
      { cause: error },
    );
  }
  return result.jsonValue();
}

/** Return visible Qt widget state exposed by the Playwright test bridge. */
export async function widgetState(page, objectName) {
  return bridgeWidgetState(page, objectName);
}

/** Click a visible text control and replace its value through keyboard input. */
export async function setWidgetText(page, objectName, text) {
  await clickWidget(page, objectName);
  await page.keyboard.press('ControlOrMeta+a');
  await page.keyboard.type(String(text), { delay: 10 });
}

/** Return a visible Qt menu action's canvas rectangle by objectName or text. */
export async function menuActionRect(page, objectNameOrText) {
  let lastError;
  for (let attempt = 0; attempt < 100; attempt += 1) {
    try {
      return await genericRect(page, `menu:${objectNameOrText}`);
    } catch (error) {
      lastError = error;
      await page.waitForTimeout(25);
    }
  }
  throw lastError;
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
  for (let attempt = 0; attempt < 50; attempt += 1) {
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
            // style.left/top are Qt's internal placement values and may be
            // expressed in backing-canvas pixels. Playwright mouse events use
            // CSS viewport pixels, so take the actual rendered DOM bounds.
            return canvas?.width && canvas?.height
              ? {
                  backingHeight: canvas.height,
                  backingWidth: canvas.width,
                  height: bounds.height,
                  width: bounds.width,
                  x: bounds.x,
                  y: bounds.y,
                }
              : null;
          })
          .filter(Boolean);
      })
      .catch(() => []);
    // A persistent Qt dialog can retain an older hidden canvas in the DOM
    // after closing.  A newly opened root menu is appended last, so callers
    // that need the currently raised root menu can request `latest`.
    const popup = popupIndex === 'latest' ? popups.at(-1) : popups[popupIndex];
    if (popup) return popup;
    await page.waitForTimeout(25);
  }
  throw new Error(`Qt popup ${popupIndex} did not expose canvas geometry`);
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
 * Return the first PNG image currently offered by the browser clipboard.
 *
 * This models the certification step of pasting into another application:
 * Playwright reads the same browser clipboard that another web application
 * would receive.  The page never fabricates image bytes; it only observes the
 * result of Sketcher's visible Copy All As -> Image action.
 */
export async function clipboardPngInfo(page) {
  return page.evaluate(async () => {
    if (!navigator.clipboard?.read) {
      throw new Error('The browser Clipboard API cannot read image content in this context.');
    }
    const items = await navigator.clipboard.read();
    for (const item of items) {
      if (!item.types.includes('image/png')) continue;
      const png = await item.getType('image/png');
      const bitmap = await createImageBitmap(png);
      try {
        return {
          height: bitmap.height,
          size: png.size,
          type: png.type,
          width: bitmap.width,
        };
      } finally {
        bitmap.close();
      }
    }
    return null;
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
