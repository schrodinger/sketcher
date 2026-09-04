import { expect } from '@playwright/test';

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
 * Call one of the bridge's exported functions and return its result.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} name - the bound function, e.g. "_sketcher_get_rect"
 * @param {string} arg
 */
async function callBridge(page, name, arg) {
  return page.evaluate(
    async ({ fn, value }) => {
      try {
        // Qt/WASM is built with Asyncify, so a bridge call that unwinds the
        // stack hands back a promise rather than its value; await covers both.
        return await Module[fn](value);
      } catch (e) {
        // A C++ exception reaches JS as an opaque emscripten value, so decode
        // it to keep the reason in the test failure. A WASM trap arrives as an
        // ordinary JS Error instead, and handing one of those to
        // getExceptionMessage() makes it read unrelated memory and take the
        // whole renderer down, so only decode what is actually opaque.
        if (e instanceof Error) {
          throw e;
        }
        let detail;
        try {
          detail = Module.getExceptionMessage(e).join(': ');
        } catch {
          detail = String(e);
        }
        throw new Error(detail);
      }
    },
    { fn: name, value: arg },
  );
}

/**
 * Resolve an object selector to {x, y, width, height, enabled} in page
 * coordinates, or null if nothing visible matches.
 *
 * This is the only way tests can locate anything inside the sketcher: Qt paints
 * the whole application into one canvas, so there are no DOM elements to query.
 * Once you have a rect, drive the application with ordinary Playwright mouse
 * and keyboard events.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} selector - "widget:<objectName>", "atom:<index>", or
 *   "bond:<index>"; indices are 0-based indices into the molecule
 */
async function getRect(page, selector) {
  const rect = JSON.parse(await callBridge(page, '_sketcher_get_rect', selector));
  return rect && rect.width !== undefined ? rect : null;
}

/**
 * Return {x, y, width, height, enabled} for a selector, failing if nothing
 * visible matches.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} selector - see getRect
 */
export async function requireRect(page, selector) {
  const rect = await getRect(page, selector);
  if (rect === null) {
    throw new Error(`Nothing visible matches "${selector}"`);
  }
  return rect;
}

/**
 * Return a widget's state, failing if no widget has that objectName.
 *
 * A button reports "enabled", "checked", "text", and "visible", which is how a
 * test asserts that a shortcut selected the tool it should have. This matches a
 * hidden widget too, since a tool inside a closed popup is still the tool the
 * shortcut is meant to have chosen.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} name - a Qt objectName
 */
export async function widgetState(page, name) {
  return requireRect(page, `state:${name}`);
}

/**
 * Whether a widget is currently on screen, and so can be clicked.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} name - a Qt objectName
 */
export async function isWidgetVisible(page, name) {
  return (await getRect(page, `widget:${name}`)) !== null;
}

/**
 * Return the {x, y} center of the drawing area (excluding toolbar and top bar).
 * @param {import('@playwright/test').Page} page
 */
export async function getDrawingAreaCenter(page) {
  const rect = await requireRect(page, 'widget:view');
  return { x: Math.round(rect.x + rect.width / 2), y: Math.round(rect.y + rect.height / 2) };
}

/**
 * Add a structure to the canvas, bypassing the import UI.
 *
 * This adds rather than replaces, so clear the canvas first if the test wants
 * the structure on its own.
 *
 * This is fixture setup, not an assertion of the Import flow; a test that
 * covers importing should drive the dialog instead.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} text - a structure in any format the sketcher auto-detects
 */
export async function loadStructure(page, text) {
  await page.evaluate((value) => Module.sketcher_import_text(value), text);
}

/**
 * Return the current molecule in the named export format.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} format - a key of Module.Format, e.g. "MDL_MOLV3000"
 */
async function getExportedText(page, format) {
  return page.evaluate((name) => Module.sketcher_export_text(Module.Format[name]), format);
}

/**
 * Return the current molecule as a SMILES string.
 */
export async function getExportedSmiles(page) {
  return getExportedText(page, 'SMILES');
}

/**
 * Return the current molecule as a HELM string.
 */
export async function getExportedHelm(page) {
  return getExportedText(page, 'HELM');
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
 * Return the system clipboard's text/plain content.
 *
 * On WASM the sketcher bypasses QClipboard and talks to navigator.clipboard
 * directly, so the page can read and write what Copy and Paste see. The
 * clipboard permissions this needs are granted in playwright.config.js.
 */
export async function getClipboardText(page) {
  return page.evaluate(() => navigator.clipboard.readText());
}

/**
 * Put text on the system clipboard as text/plain.
 */
export async function setClipboardText(page, text) {
  await page.evaluate((value) => navigator.clipboard.writeText(value), text);
}

/**
 * Show a test-only cursor marker at a page coordinate.
 *
 * Playwright's trace viewer can't show where a click landed inside a canvas,
 * so this draws the pointer position when PLAYWRIGHT_SHOW_MOUSE=1. It is off
 * by default and is only an aid for writing and debugging tests.
 */
async function showMouseMarker(page, x, y) {
  if (process.env.PLAYWRIGHT_SHOW_MOUSE !== '1') {
    return;
  }
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

/**
 * Hide the optional cursor marker, so that it stays out of a screenshot.
 */
export async function hideMouseMarker(page) {
  if (process.env.PLAYWRIGHT_SHOW_MOUSE !== '1') {
    return;
  }
  await page.evaluate(() => {
    const marker = document.getElementById('playwright-mouse-marker');
    if (marker) {
      marker.style.display = 'none';
    }
  });
}

/**
 * Press and release a mouse button at a page coordinate, optionally with
 * keyboard modifiers held down.
 *
 * Playwright's page.mouse.click() takes no modifiers, and Qt wants to see the
 * pointer arrive before the press, so this moves first and holds the modifiers
 * around the whole gesture.
 *
 * @param {import('@playwright/test').Page} page
 * @param {number} x
 * @param {number} y
 * @param {object} [options]
 * @param {'left'|'right'|'middle'} [options.button]
 * @param {string[]} [options.modifiers] - e.g. ['Shift']
 */
export async function mouseClick(page, x, y, { button = 'left', modifiers = [] } = {}) {
  await showMouseMarker(page, x, y);
  for (const modifier of modifiers) {
    await page.keyboard.down(modifier);
  }
  try {
    await page.mouse.move(x, y, { steps: 4 });
    await page.mouse.down({ button });
    await page.mouse.up({ button });
  } finally {
    for (const modifier of [...modifiers].reverse()) {
      await page.keyboard.up(modifier);
    }
  }
}

/**
 * Drag from one page coordinate to another with intermediate moves, so that Qt
 * sees a real drag rather than a teleport.
 *
 * @param {import('@playwright/test').Page} page
 * @param {{x: number, y: number}} start
 * @param {{x: number, y: number}} end
 * @param {object} [options]
 * @param {'left'|'right'|'middle'} [options.button]
 * @param {string[]} [options.modifiers]
 * @param {number} [options.steps]
 */
export async function mouseDrag(
  page,
  start,
  end,
  { button = 'left', modifiers = [], steps = 12 } = {},
) {
  await showMouseMarker(page, start.x, start.y);
  for (const modifier of modifiers) {
    await page.keyboard.down(modifier);
  }
  try {
    await page.mouse.move(start.x, start.y, { steps: 4 });
    await page.mouse.down({ button });
    for (let step = 1; step <= steps; step += 1) {
      const progress = step / steps;
      await page.mouse.move(
        start.x + (end.x - start.x) * progress,
        start.y + (end.y - start.y) * progress,
      );
    }
    await page.mouse.up({ button });
  } finally {
    for (const modifier of [...modifiers].reverse()) {
      await page.keyboard.up(modifier);
    }
  }
  await showMouseMarker(page, end.x, end.y);
}

/**
 * Wait for a selector to resolve to a visible, enabled rect and return it.
 *
 * Qt updates visibility and enabled state in response to model changes that may
 * not have landed yet, and a menu popup is laid out an event-loop turn after
 * the click that opens it. Failing on disabled is deliberate: a real user
 * couldn't have clicked it, and a test that clicks a disabled control then
 * asserts nothing happened would pass either way. To assert unavailability,
 * poll `(await requireRect(page, selector)).enabled` instead.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} selector - see getRect
 */
async function waitForClickable(page, selector) {
  let rect;
  await expect
    .poll(
      async () => {
        rect = await getRect(page, selector);
        if (rect === null) {
          return 'not visible';
        }
        return rect.enabled ? 'clickable' : 'disabled';
      },
      { timeout: 5000, message: `"${selector}"` },
    )
    .toBe('clickable');
  return rect;
}

/**
 * Click whatever a selector resolves to, with a real mouse event so that Qt's
 * own hit-testing runs.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} selector - see getRect
 * @param {object} [options] - see mouseClick
 */
async function click(page, selector, options) {
  const rect = await waitForClickable(page, selector);
  await mouseClick(page, rect.x + rect.width / 2, rect.y + rect.height / 2, options);
}

/**
 * Click a widget by its Qt objectName (e.g. "c_btn").
 *
 * This works for a widget inside a Qt::Popup as well as one in the main window,
 * because the bridge maps rects through global coordinates.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} name - Qt objectName of the widget
 * @param {object} [options] - see mouseClick
 */
export async function clickWidget(page, name, options) {
  await click(page, `widget:${name}`, options);
}

// ToolButtonWithPopup opens its popup 250 ms after the press, and a release
// before then cancels it, so hold well past that.
const POPUP_HOLD_MS = 600;

/**
 * Click a tool that lives in a button's press-and-hold popup.
 *
 * The bridge is asked which button owns the popup so that tests don't have to
 * carry a map of every tool to the button hiding it. The popup is a plain
 * widget rather than a QMenu, so it doesn't run a nested event loop and its
 * contents can be located and clicked normally once it's open.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} name - Qt objectName of the tool inside the popup
 */
export async function clickPopupTool(page, name) {
  const owner = await callBridge(page, '_sketcher_get_popup_owner', name);
  if (!owner) {
    throw new Error(`"${name}" is not visible and is not inside a popup`);
  }
  const rect = await waitForClickable(page, `widget:${owner}`);
  const x = rect.x + rect.width / 2;
  const y = rect.y + rect.height / 2;
  await showMouseMarker(page, x, y);
  await page.mouse.move(x, y, { steps: 4 });
  await page.mouse.down();
  await page.waitForTimeout(POPUP_HOLD_MS);
  await page.mouse.up();
  await click(page, `widget:${name}`);
}

/**
 * Click an atom or bond by its index in the molecule.
 *
 * @param {import('@playwright/test').Page} page
 * @param {'atom'|'bond'} kind - monomers are addressed as 'atom' and monomer
 *   connectors as 'bond'
 * @param {number} index - 0-based index into the molecule's atoms or bonds
 * @param {object} [options] - see mouseClick
 */
export async function clickItem(page, kind, index, options) {
  await click(page, `${kind}:${index}`, options);
}

/**
 * Move the pointer over whatever a selector resolves to.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} selector - see getRect
 */
async function hover(page, selector) {
  const rect = await waitForClickable(page, selector);
  const x = rect.x + rect.width / 2;
  const y = rect.y + rect.height / 2;
  await showMouseMarker(page, x, y);
  await page.mouse.move(x, y, { steps: 4 });
}

/**
 * Open a context menu by right-clicking whatever a selector resolves to, then
 * walk a path of rows: every row but the last is hovered to open its submenu,
 * and the last is clicked.
 *
 * Unlike the More Actions menu, a context menu can be driven for real. The
 * sketcher shows one with QMenu::show() rather than exec(), so it doesn't run
 * the nested event loop that would otherwise suspend the WASM module — see
 * activateAction.
 *
 * The whole gesture is retried, because Qt/WASM can take a frame to put the
 * popup up and a half-open menu resolves no rows. Reopening replays the same
 * user actions rather than falling back to triggering the action directly, so a
 * genuinely unreachable row still fails.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} selector - what to right-click, see getRect
 * @param {...string} path - row labels, outermost first
 */
export async function contextMenuAction(page, selector, ...path) {
  if (path.length === 0) {
    throw new Error('contextMenuAction needs at least one row to click');
  }
  const rect = await waitForClickable(page, selector);
  let lastError;
  for (let attempt = 0; attempt < 3; attempt += 1) {
    await mouseClick(page, rect.x + rect.width / 2, rect.y + rect.height / 2, {
      button: 'right',
    });
    try {
      for (const label of path.slice(0, -1)) {
        await hover(page, `menu:${label}`);
      }
      await click(page, `menu:${path.at(-1)}`);
      return;
    } catch (error) {
      lastError = error;
      await page.keyboard.press('Escape');
    }
  }
  throw lastError;
}

/**
 * Open a context menu and leave it showing, having hovered down to the last row
 * of the path. Use this to screenshot a menu rather than act on it.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} selector - what to right-click, see getRect
 * @param {...string} path - row labels to hover, outermost first
 */
export async function openContextMenu(page, selector, ...path) {
  const rect = await waitForClickable(page, selector);
  let lastError;
  for (let attempt = 0; attempt < 3; attempt += 1) {
    await mouseClick(page, rect.x + rect.width / 2, rect.y + rect.height / 2, {
      button: 'right',
    });
    try {
      for (const label of path) {
        await hover(page, `menu:${label}`);
      }
      return;
    } catch (error) {
      lastError = error;
      await page.keyboard.press('Escape');
    }
  }
  throw lastError;
}

/**
 * Trigger a menu action by objectName or visible text, without opening its
 * menu.
 *
 * This is the one control that can't be clicked for real. Qt runs a nested
 * event loop while a QToolButton's menu is open, and because Qt/WASM is built
 * with Asyncify that leaves the WebAssembly stack suspended — no call into the
 * module returns until the menu closes, so a test can't even ask where a row
 * is. Everything else, popup widgets included, gets a real mouse event.
 *
 * Throws if no action matches or the match is disabled.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} nameOrText - objectName, or the row's visible text
 */
export async function activateAction(page, nameOrText) {
  await callBridge(page, '_sketcher_activate_action', nameOrText);
}

/**
 * Click a visible text control and replace its value through keyboard input.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} name - Qt objectName of the text control
 * @param {string} text
 */
export async function setWidgetText(page, name, text) {
  await clickWidget(page, name);
  await page.keyboard.press('ControlOrMeta+a');
  await page.keyboard.type(String(text), { delay: 10 });
}

/**
 * Place a single atom of the given element in the middle of the drawing area,
 * selecting the element with its keyboard shortcut.
 *
 * @param {import('@playwright/test').Page} page
 * @param {string} element - an element symbol with a shortcut, e.g. "N"
 */
export async function drawElement(page, element) {
  const center = await getDrawingAreaCenter(page);
  await focusCanvas(page);
  await page.keyboard.press(element);
  await mouseClick(page, center.x, center.y);
}

/**
 * Draw one bond with the active bond tool, starting from the middle of the
 * drawing area.
 */
export async function drawBond(page) {
  const center = await getDrawingAreaCenter(page);
  await mouseDrag(page, center, { x: center.x + 100, y: center.y });
}
