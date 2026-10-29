import { expect } from '@playwright/test';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

export const FIXTURES = path.join(path.dirname(fileURLToPath(import.meta.url)), 'fixtures');

/**
 * Navigate to the sketcher page and wait for the WASM module to be fully
 * loaded and the canvas to be visible.
 */
export async function waitForSketcherReady(page) {
  // Window activation must not change the hovered atom/bond for a shortcut.
  await page.addInitScript(() => {
    document.addEventListener('pointermove', ({ clientX, clientY }) => {
      window.__sketcherPlaywrightPointer = { x: clientX, y: clientY };
    });
  });
  await page.goto('/wasm_shell.html');
  await page.waitForFunction(() => typeof window.Module !== 'undefined', undefined, {
    timeout: 20000,
  });
  // Wait for the canvas inside the shadow DOM to be attached and visible
  const canvas = page.locator('#screen canvas');
  await canvas.waitFor({ state: 'visible', timeout: 10000 });
}

/**
 * Activate a Qt window through its inert top edge, preserving the pointer.
 * Browser focus alone does not reactivate a Qt 6.7 window after popup closure.
 */
export async function activateQtWindow(page, qtWindow) {
  const pointer = await page.evaluate(() => window.__sketcherPlaywrightPointer);
  const bounds = await qtWindow.boundingBox();
  await page.mouse.click(bounds.x + bounds.width / 2, bounds.y + 5);
  if (pointer) await page.mouse.move(pointer.x, pointer.y);
}

/** Focus the browser element Qt uses to receive keyboard input. */
export async function focusQtWindow(page, qtWindow) {
  const input = qtWindow.locator('input.qt-window-input-element');
  await ((await input.count()) ? input : qtWindow.locator('canvas')).focus();
}

/**
 * Focus the drawing window so keyboard events reach Qt in both Qt 6.7 and 6.12.
 */
export async function focusCanvas(page) {
  const qtWindow = page.locator('[id^="qt-window-"]').first();
  const alreadyFocused = await qtWindow.evaluate((element) =>
    element.contains(element.getRootNode().activeElement),
  );
  if (!alreadyFocused) await activateQtWindow(page, qtWindow);
  await focusQtWindow(page, qtWindow);
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
 * @param {string} selector - one of:
 *   "widget:<objectName>"  a visible widget, by its Qt objectName
 *   "state:<objectName>"   the same, but matching a hidden widget too, and
 *                          adding "visible" to the result
 *   "atom:<index>"         an atom, by 0-based index into the molecule
 *   "bond:<index>"         a bond, by 0-based index into the molecule
 *   "menu:<name or text>"  a row of an already-open context menu
 *   A widget whose text the user can read also reports "text" — a caption, a
 *   combo box's current item, or a line edit's or spin box's contents — and a
 *   button also reports "checked" and "toolTip".
 */
export async function getRect(page, selector) {
  const encoded = await callBridge(page, '_sketcher_get_rect', selector);
  // Qt/WASM can yield an empty result for one browser turn while it closes a
  // modal or lays out its replacement. Treat that as temporarily unavailable
  // so the normal clickable-state poll can retry instead of parsing it as JSON.
  if (typeof encoded !== 'string' || encoded.trim() === '') return null;
  const rect = JSON.parse(encoded);
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
 * Every widget reports "enabled" and "visible", and one the user can read text
 * from adds "text"; a button adds "checked" and "toolTip", which is how a test
 * asserts that a shortcut selected the tool it should have. This matches a
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
 * Select all items on the canvas and wait for the selection controls to
 * reflect the completed selection.
 */
export async function selectAll(page) {
  await clickWidget(page, 'select_all_btn');
  await expect.poll(async () => (await widgetState(page, 'select_all_btn')).enabled).toBe(false);
  await expect
    .poll(async () => (await widgetState(page, 'clear_selection_btn')).enabled)
    .toBe(true);
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
 *
 * The mouse helpers deliberately leave the marker showing once a gesture ends,
 * so that a trace frame captured after the click still says where the pointer
 * went; hiding it there would erase the only record of it. Call this from a test
 * that is about to take a visual snapshot.
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
export async function waitForClickable(page, selector) {
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
export async function clickSelector(page, selector, options) {
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
  await clickSelector(page, `widget:${name}`, options);
}

// ToolButtonWithPopup::m_popup_delay is 250 ms, and releasing before it elapses
// cancels the popup, so hold well past it. Raise this if that default changes.
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
  // When the owner is embedded in a cascading QMenu, Qt does not create that
  // submenu until its hover delay expires. Waiting before resolving the owner
  // prevents a same-named sidebar button from winning the lookup meanwhile.
  await page.waitForTimeout(500);
  const rect = await waitForClickable(page, `widget:${owner}`);
  const x = rect.x + rect.width / 2;
  const y = rect.y + rect.height / 2;
  await showMouseMarker(page, x, y);
  await page.mouse.move(x, y, { steps: 4 });
  await page.mouse.down();
  try {
    await page.waitForTimeout(POPUP_HOLD_MS);
    const target = await waitForClickable(page, `widget:${name}`);
    const targetX = target.x + target.width / 2;
    const targetY = target.y + target.height / 2;
    await showMouseMarker(page, targetX, targetY);
    // A ToolButtonWithPopup closes its press-and-hold popup if the original
    // press is released before the target is activated. Keep that press held
    // while clicking the desired tool, matching a person's press-drag-click
    // sequence and the interaction Squish performs.
    await page.mouse.move(targetX, targetY, { steps: 4 });
    await page.mouse.down();
    await page.mouse.up();
  } finally {
    await page.mouse.up();
  }
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
  await clickSelector(page, `${kind}:${index}`, options);
}

/**
 * Move the pointer over whatever a selector resolves to.
 *
 * Every selector kind reports "enabled", not just a button, so this waits for an
 * enabled target the way clickSelector() does — its callers hover a menu row to
 * open the submenu underneath it, and a disabled row has none to open.
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
      await clickSelector(page, `menu:${path.at(-1)}`);
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

export async function openSketcher(page) {
  // The wrapper's WASM adapter imports these helpers, so load it on demand.
  const { Sketcher } = await import('./sketcher/wrappers/sketcher.js');
  const sk = new Sketcher(page);
  await sk.open();
  await expect
    .poll(() => page.evaluate(() => typeof Module._sketcher_inspect_scene))
    .toBe('function');
  return sk;
}

export const center = (rect) => ({
  x: rect.centerX ?? rect.x + rect.width / 2,
  y: rect.centerY ?? rect.y + rect.height / 2,
});

// Inspect the foreground scene, including a modal editor, without mutating it.
export async function scene(page) {
  let encoded;
  await expect
    .poll(async () => {
      encoded = await page.evaluate(async () => await Module._sketcher_inspect_scene());
      // Like getRect, wait through the empty result Qt can return while an
      // asynchronous import or dialog transition is still completing.
      return typeof encoded === 'string' && encoded.trim() !== '';
    })
    .toBe(true);
  const state = JSON.parse(encoded);
  state.atoms.sort((a, b) => a.index - b.index);
  return state;
}

export async function waitAtoms(page, count) {
  await expect.poll(async () => (await scene(page)).atoms.length).toBe(count);
  return (await scene(page)).atoms;
}

export async function point(page, index) {
  const atom = (await scene(page)).atoms.find((item) => item.index === index);
  expect(atom, `visible atom ${index}`).toBeTruthy();
  return center(atom);
}

export async function canvasPoint(page, dx = 0, dy = 0) {
  const view = await requireRect(page, 'widget:view');
  return { x: view.x + view.width / 2 + dx, y: view.y + view.height / 2 + dy };
}

export async function place(page, dx = 0, dy = 0) {
  const target = await canvasPoint(page, dx, dy);
  await mouseClick(page, target.x, target.y);
}

export async function dialogTitle(page) {
  return (await getRect(page, 'dialog:title'))?.text ?? null;
}

export async function expectDialog(page, title, timeout = 5000) {
  await expect.poll(() => dialogTitle(page), { timeout }).toBe(title);
}

export async function closeDialog(page, button = 'OK') {
  await clickSelector(page, `dialog:${button}`);
}

export async function pasteText(page, text) {
  const { Sketcher } = await import('./sketcher/wrappers/sketcher.js');
  // Reuse the existing real Import-menu and clipboard workflow. Confirm via
  // the foreground dialog's measured OK button, including in a nested editor.
  await new Sketcher(page).import_menu('paste_in_text', text, false);
  await expect.poll(async () => (await widgetState(page, 'structure_text_edit')).text).toBe(text);
  await closeDialog(page);
}

export async function replaceEditor(page, smiles) {
  if ((await scene(page)).atoms.length > 0) await clear(page);
  await pasteText(page, smiles);
  await clickWidget(page, 'fit_btn');
}

export async function holdPopup(page, button) {
  const target = center(await waitForClickable(page, `widget:${button}`));
  await page.mouse.move(target.x, target.y);
  await page.mouse.down();
  // ToolButtonWithPopup uses a 250 ms timer; this wait is the actual gesture.
  await page.waitForTimeout(POPUP_HOLD_MS);
}

export async function hoverAtom(page, index) {
  const target = await point(page, index);
  // A placement can leave the pointer at this exact coordinate. Leave and
  // re-enter so Qt receives a hover event for the newly created item.
  await page.mouse.move(1100, 40);
  await page.mouse.move(target.x, target.y, { steps: 4 });
}

export async function undo(page) {
  await clickWidget(page, 'undo_btn');
}

export async function redo(page) {
  await clickWidget(page, 'redo_btn');
}

export async function clear(page) {
  await clickWidget(page, 'clear_btn');
  await waitAtoms(page, 0);
}

export async function moveAtom(page, index, dx, dy) {
  if (!(await widgetState(page, 'select_tool_btn')).checked) {
    await clickWidget(page, 'select_tool_btn');
  }
  const target = await point(page, index);
  await mouseClick(page, target.x, target.y);
  await clickWidget(page, 'move_rotate_btn');
  await mouseDrag(page, target, { x: target.x + dx, y: target.y + dy });
  await clickWidget(page, 'clear_selection_btn');
  await clickWidget(page, 'select_tool_btn');
}

export async function popupInfo(page, button) {
  await holdPopup(page, button);
  try {
    let popup;
    await expect
      .poll(async () => {
        popup = (await scene(page)).popup;
        return popup.buttons?.length ?? 0;
      })
      .toBeGreaterThan(0);
    popup.buttons.sort((a, b) => a.y - b.y || a.x - b.x);
    return popup;
  } finally {
    await page.mouse.up();
  }
}

export async function popupButtons(page, button) {
  return (await popupInfo(page, button)).buttons;
}

export async function dismissPopup(page) {
  await page.mouse.up();
  const view = await requireRect(page, 'widget:view');
  await mouseClick(page, view.x + view.width / 2, view.y - 10);
}

export async function choosePopup(page, button, symbol, { scroll = false } = {}) {
  await holdPopup(page, button);
  try {
    let target;
    await expect
      .poll(
        async () => {
          const popup = (await scene(page)).popup;
          target = popup.buttons?.find((item) => item.text === symbol);
          if (!target) return false;
          const pos = center(target);
          const inside =
            pos.y > popup.y &&
            pos.y < popup.y + popup.height &&
            pos.x > popup.x &&
            pos.x < popup.x + popup.width;
          if (!inside && scroll) {
            await page.mouse.move(popup.x + popup.width / 2, popup.y + popup.height / 2);
            await page.mouse.wheel(0, pos.y < popup.y ? -900 : 900);
          }
          return inside;
        },
        { timeout: 15000, message: `${button} contains reachable ${symbol}` },
      )
      .toBe(true);
    const pos = center(target);
    await mouseClick(page, pos.x, pos.y);
  } finally {
    await page.mouse.up();
  }
}

export async function screenshotRegion(page, name, rect) {
  await expect(page).toHaveScreenshot(name, {
    clip: { x: rect.x, y: rect.y, width: rect.width, height: rect.height },
    maxDiffPixelRatio: 0.002,
  });
}

export async function screenshotView(page, name) {
  await screenshotRegion(page, name, await requireRect(page, 'widget:view'));
}

export async function expectCursor(page, name) {
  const pos = await canvasPoint(page);
  await page.mouse.move(pos.x, pos.y);
  let base64;
  await expect
    .poll(async () => {
      base64 = await page.evaluate(() => Module._sketcher_get_cursor_png());
      return base64.length;
    })
    .toBeGreaterThan(0);
  expect(Buffer.from(base64, 'base64')).toMatchSnapshot(name, { maxDiffPixelRatio: 0.002 });
}

export async function relativeLayout(page, order) {
  const atoms = (await scene(page)).atoms;
  const ordered = order ? order.map((index) => atoms.find((atom) => atom.index === index)) : atoms;
  const origin = ordered[0];
  return ordered.map((atom) => ({
    x: atom.sceneX - origin.sceneX,
    y: atom.sceneY - origin.sceneY,
  }));
}

export async function expectLayout(page, expected, order) {
  await expect
    .poll(async () => {
      const actual = await relativeLayout(page, order);
      if (actual.length !== expected.length) return Infinity;
      return Math.max(
        ...actual.map((pos, index) =>
          Math.hypot(pos.x - expected[index].x, pos.y - expected[index].y),
        ),
      );
    })
    .toBeLessThan(0.1);
}
