# Sketcher: Qt Removal & React Migration Plan

- [1. Goals](#1-goals)
- [2. Key decisions](#2-key-decisions)
- [3. Starting point](#3-starting-point)
- [4. Target architecture](#4-target-architecture)
- [5. Phased plan](#5-phased-plan)
- [6. Tech stack](#6-tech-stack)
- [7. Risks](#7-risks)

## 1. Goals

- Remove Qt from the sketcher, apart from a thin shim for Qt hosts.
- Give every target that shows a UI a single frontend that looks and behaves the same everywhere.
- Keep every consumer working:
  - Maestro: its dialogs, and the ligand overlay in the 3D workspace.
  - About 105 Python files that use the PyQt panels through SIP.
  - The drop-into-a-webpage WASM build, which should end up smaller and faster than today's.
  - The standalone desktop app.
- Keep headless image generation, with no browser, Qt or X libraries (or any display server) needed.
- Make it cheap to have many sketchers on one page, as LiveDesign does: the WASM module is downloaded and compiled once, and no sketcher needs its own iframe.

Matching today's Qt look is not required. The rewrite is a chance to refresh the UI if we want one.

## 2. Key decisions

### 2.1 Chemistry stays in C++ (no RDKit.js / MinimalLib)

- **Consumers need C++ RDKit objects, in-process and synchronously.** The public `SketcherWidget` API passes them directly: `getRDKitMolecule()` returns an `RDKit::ROMol`, selection is a set of `RDKit::Atom*`/`Bond*`, and `atomHovered(const RDKit::Atom*)` is a signal. MinimalLib can't provide any of that.
- **MinimalLib's API is too narrow.** It covers whole-molecule I/O, depiction and descriptors. The sketcher edits molecules atom by atom and bond by bond, and relies on substance groups, stereo handling, and HELM, FASTA and the monomer database in `rdkit_extensions`. All of these are built on C++ APIs that MinimalLib doesn't expose.
- **One implementation.** Desktop, Qt hosts and headless rendering need the C++ chemistry anyway. Moving to RDKit.js for the web would mean keeping two chemistry layers in sync.

### 2.2 Web UI in a webview vs a native toolkit (Slint)

A single UI on every target means either a web UI (React, HTML and CSS) in a webview everywhere, or a native cross-platform toolkit like Slint.

- **Why the web UI wins:**
  - **Ecosystem:** mature component libraries for the toolbars, menus and dialogs (palettes, HELM input, monomer browsers), including Crayon, the UI library LiveDesign already uses.
  - **Layout:** browsers are very good at laying out and styling static UI, so we stop hand-building it.
  - **Testability:** a DOM UI is far easier to drive and inspect from Playwright than today's WASM canvas.
  - **Composability:** embedders could supply their own components, such as LiveDesign's monomer picker. This is an opportunity, not a committed extension API.
  - **Speed and hiring:** iteration speed and how easy it is to hire for.
  - **Bundle size:** not a real differentiator, because the RDKit WASM core (about 10–15 MB) dwarfs either UI runtime.
- **What Slint does better:** a smaller desktop footprint, and a cleaner way to embed in Qt hosts.
- **Gate:** whether a webview is acceptable inside Qt hosts is decided by the Qt-host gate in Phase 6.

### 2.3 Interaction lives in a C++ controller

The core owns not just chemistry but the interaction logic. A `controller/` module holds the active tool, handles pointer and key input, does hover and hit-testing, produces drag previews, and owns selection. Scene tools need the C++ model anyway (chain placement, ring fusion, angle snapping, hit priority), so they're ported, not rewritten in TS.

The controller also keeps an **action registry**: each action's id, label, icon, shortcut and enabled/checked state, and which actions belong in the context menu for a given atom, bond, selection or background. It knows nothing about layout or styling.

Frontends are thin:

- They lay out and style the toolbars, menus and dialogs, using the action registry for what to show.
- They paint the display list.
- They forward input events and action invocations to the controller.
- Dialogs are frontend components that send commands. The core doesn't describe dialog layouts.

Because frontends are this thin, a native Qt frontend over the same controller is a realistic fallback for Qt hosts, not a second implementation.

### 2.4 Depiction geometry lives in C++

A new Qt-free `depict/` module computes the geometry: label layout, bond trimming, double-bond offsets, wedges, monomer shapes and hit-test shapes. It emits a **display list** of lines, polygons, paths and positioned text. Text is measured through a `FontMetrics` interface backed by the bundled Arimo TTF, so every target lays out text the same way.

The display list is **retained**, not rebuilt each frame. Primitives are grouped into items keyed by atom, bond, monomer or other scene object, and the core emits diffs (items added, changed or removed). A drag sends the handful of items that moved, not the whole molecule. This is the same pattern as three.js's scene graph, and it keeps the canvas off React's render path.

Thin painters turn the display list into output:

| Painter              | Used for                                |
| -------------------- | --------------------------------------- |
| Browser painter (TS) | The interactive UI on every webview     |
| SVG writer (C++)     | Headless image generation               |
| QPainter (Qt, thin)  | The canvas-only overlay view in Maestro |

The browser painter is either SVG DOM (one `<g>` per item, which fits retained diffs and lets Playwright target individual atoms) or Canvas2D. Both are prototyped in Phase 4 and the choice is made by benchmark.

Keeping geometry in C++ means headless rendering and every UI share a single source, which removes the "it looks different in the sketcher than in the compound image" class of bugs. Putting it in TS would leave no browser-free renderer for the image API.

### 2.5 Alternatives rejected

- **Electron:** its main process is Node, and our core is C++. Electron would need a native module or a child process to reach the core. With `webview/webview`, the C++ binary _is_ the host.
- **CEF:** bundles Chromium, so rendering is identical on every platform, at the cost of about 100 MB. It stays the fallback if the system webviews turn out to render the canvas differently from each other.
- **React Native (desktop):** react-native-windows and react-native-macos have no Linux target, which rules it out for Maestro and the Linux standalone app.

## 3. Starting point

| Layer                             | LOC   | Qt coupling                             | Fate                                      |
| --------------------------------- | ----- | --------------------------------------- | ----------------------------------------- |
| `rdkit_extensions/`               | 16.3K | None                                    | Keep                                      |
| `sketcher/rdkit/`                 | 7.9K  | Minimal                                 | Move to `sketcher_core`                   |
| `sketcher/model/`                 | 7.5K  | `QObject`, `QUndoStack`, a few Qt types | Port to `sketcher_core`                   |
| `sketcher/molviewer/` (rendering) | 11.6K | QPainter/QGraphicsItem throughout       | Geometry → `depict/`; the rest is deleted |
| `sketcher/tool/` (36 scene tools) | 8.0K  | `QObject`                               | Port to `controller/` in the core         |
| `widget/` + `dialog/` + `menu/`   | 11.3K | All QWidget                             | Rewrite in React                          |

- **WASM ships the whole Qt app today:** 44 MB raw, 13 MB gzipped. A lean core should land around 10–15 MB raw, about 3× smaller.
- **Regression net (already in place):**
  - The sketcher C++ tests, including the public-contract tests for the widget and image generation.
  - The Maestro and Python consumer tests in mmshare.
  - The Playwright suite: 41 files and about 147 screenshot checkpoints, written against page-object wrappers.

## 4. Target architecture

Each consumer embeds the same React bundle, and React talks to `sketcher_core`. The one exception is the Maestro overlay, which is select-only (it can select atoms but not edit), has no toolbars, menus or dialogs, and calls the core directly.

```
┌────────────────┐ ┌────────────────┐ ┌─────────────────────────────────────────┐
│ Web page       │ │ Standalone app │ │ Maestro · PyQt panels                   │
│                │ │ (no Qt)        │ │ via the sketcher Qt shim                │
└───────┬────────┘ └───────┬────────┘ └────────┬───────────────────────┬────────┘
┌───────▼────────┐ ┌───────▼────────┐ ┌────────▼─────────┐ ┌───────────▼────────┐
│ browser        │ │ system webview │ │ SketcherWidget   │ │ SketcherView       │
│                │ │ webview/webview│ │ QWebEngineView   │ │ (Maestro overlay)  │
└───────┬────────┘ └───────┬────────┘ └────────┬─────────┘ │ native QWidget,    │
┌───────▼──────────────────▼───────────────────▼─────────┐ │ QPainter painter,  │
│ React UI (TypeScript), one bundle for all webviews     │ │ select-only,       │
│   display-list painter · toolbars, menus, dialogs      │ │ no React           │
│   protocol client (same code on every host)            │ │                    │
└───────┬──────────────────┬───────────────────┬─────────┘ └───────────┬────────┘
┌───────▼────────┐ ┌───────▼────────┐ ┌────────▼─────────┐             │
│ embind         │ │ bind / eval    │ │ QWebChannel      │             │
└───────┬────────┘ └───────┬────────┘ └────────┬─────────┘             │
        └──────────────────┼───────────────────┘                       │
      commands + input in / events + display-list diffs out   direct C++ calls
┌──────────────────────────▼───────────────────────────────────────────▼────────┐
│ sketcher_core (new, no Qt)                                                    │
│   protocol/   command dispatcher (same code on every host)                    │
│   controller/ tools, input, hover/hit-testing, selection, action registry     │
│   model/      MolModel, edit commands, undo stack, events                     │
│   depict/     retained display list, FontMetrics (Arimo), headless SVG        │
│   rdkit/      sketcher chemistry helpers                                      │
└──────────────────────────────────────┬────────────────────────────────────────┘
┌──────────────────────────────────────▼────────────────────────────────────────┐
│ rdkit_extensions (unchanged): HELM, FASTA, SMILES, MOL, monomer DB            │
└──────────────────────────────────────┬────────────────────────────────────────┘
                                  RDKit (C++)
```

**Libraries:** there are three.

- `rdkit_extensions` stays as it is.
- `sketcher_core` has no Qt and is linked by every host. It has no mutable global state, so any number of sketchers can share one process or one WASM instance.
- `sketcher` becomes a thin Qt shim containing `SketcherWidget` and `SketcherView`.

**Protocol:** there is one protocol, with one implementation on each side.

- **C++ side:** a command dispatcher in `sketcher_core`.
- **TS side:** a protocol client in React, written against a small `Transport` interface (send a command, subscribe to events).
- **Per host:** each transport adapter is a few dozen lines on each side.
  - embind (browser) is synchronous and calls into the WASM core.
  - `bind`/`eval` (standalone) is the `webview/webview` API. It exposes the native, in-process C++ core as async JS functions; there's no WASM.
  - QWebChannel (Qt hosts) is asynchronous.
  - The client treats every call as async, so React code is the same on every host.
- **Two channels:**
  - Commands and events are JSON, with TS types generated from the spec.
  - Display-list diffs are a separate stream behind their own interface. They start as JSON and move to a binary layout only if the Phase 4 benchmark needs it. SharedArrayBuffer is out: it needs cross-origin isolation headers on every embedding page. In the browser, JS can read typed-array views straight from WASM memory without copying.
- **Qt hosts' synchronous calls:** calls such as `getRDKitMolecule()` go straight from the shim to the core, not through React.

**Split of responsibilities:**

- **C++ owns** chemistry, edit semantics, undo, selection, the active tool and its interaction logic, hover and hit-testing, drag previews, the action registry and depiction geometry.
- **Frontends own** layout and styling of the toolbars, menus and dialogs, painting the display list, and forwarding input and action invocations.
- **Crossing the boundary:** pointer moves cross into the core. Over embind this is a synchronous call that costs microseconds. Over the async transports, the client coalesces `pointermove` to one per animation frame with at most one in flight.

## 5. Phased plan

Every PR lands on `main` and leaves the existing Qt app working until Phase 8. Phases 1–5 pay off whatever the Qt-host gate in Phase 6 decides.

### Phase 1 — Extract `sketcher_core` (model + rdkit)

Make the chemistry and model layer Qt-free behind a new library, without changing behavior. `SketcherModel` stays in the Qt library for now. Its UI state moves out later: the active tool to the controller in Phase 3, and display settings and preferences to React in Phase 5.

| PR  | Change                                                                                                                                                                                                   |
| --- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | Untangle dependencies: move the chemistry enums, coordinate helpers, monomer utilities and constants out of Qt headers, and move the arrowhead-placement depiction code into molviewer.                  |
| 2   | Remove incidental Qt types from `MolModel`: `QString` becomes `std::string`, `QColor` becomes an RGBA struct, `Q_ASSERT` becomes `assert`. Also remove `QString` from `rdkit/`.                          |
| 3   | Replace `MolModel`'s Qt signals with `boost::signals2`, using `scoped_connection` so listeners disconnect automatically. A thin QObject adapter re-emits the signals so existing consumers don't change. |
| 4   | Replace `QUndoStack`/`QUndoCommand` with a small core undo stack that supports macros, merging and can-undo/redo notifications.                                                                          |
| 5   | Create the `sketcher_core` CMake target with no Qt link, move the files, and add the export macro, install rules and recipe `install_name_tool` lines.                                                   |
| 6   | Move the model and rdkit tests onto `sketcher_core` alone, which proves the core is Qt-free.                                                                                                             |

**Done when:** the app behaves identically, all existing tests pass, the core tests build without Qt, and the core has no mutable globals.

### Phase 2 — Extract depiction into `depict/`

Move the geometry out of the molviewer `paint()` methods into the display-list producer. The existing QGraphicsItems become thin painters of the display list, so today's Qt app is the first consumer. **Port faithfully, with no visual changes.** The existing Playwright screenshots must pass unchanged after every PR.

| PR  | Change                                                                                                                                                                           |
| --- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | Display-list types (items keyed by scene object, with added/changed/removed diffs), the `FontMetrics` interface and the Arimo implementation, with unit tests. No consumers yet. |
| 2   | Atom labels: layout and hit shapes move to `depict/`, and the atom item paints the display list.                                                                                 |
| 3   | Bonds: trimming, double-bond offsets, wedges and hashes.                                                                                                                         |
| 4   | Monomers: shapes, connectors and attachment-point labels.                                                                                                                        |
| 5   | Everything else: arrows, pluses, S-group brackets, and selection and hover highlights.                                                                                           |
| 6   | SVG writer and a Qt-free SVG image API, added next to the existing Qt-typed functions.                                                                                           |

**Done when:** every visual item draws from the display list, headless SVG works without Qt or X libraries, and the screenshots are unchanged.

### Phase 3 — Extract interaction into `controller/`

Port the scene tools into the core. The Qt scene becomes a thin adapter that forwards input to the controller and paints its display list, so today's Qt app is again the first consumer. **No behavior changes.** The Playwright suite must pass unchanged after every PR, and each PR deletes the Qt tool code it replaces.

| PR  | Change                                                                                                                                                                                                                                               |
| --- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | Controller skeleton: input event types (pointer, key, modifiers), the tool interface, active-tool state, hover and hit-testing against the display-list hit shapes, and drag previews as display-list items. The Qt scene forwards its events to it. |
| 2   | Action registry: ids, labels, icons, shortcuts, enabled/checked state, and context-menu contents per context. The Qt toolbars and menus read from it.                                                                                                |
| 3   | Select (rectangle, lasso and fragment modes), move/rotate, and erase.                                                                                                                                                                                |
| 4   | Drawing: atoms, bonds, chains, rings and custom fragments, charge and explicit H.                                                                                                                                                                    |
| 5   | Enumeration and reactions: R-groups, attachment points, atom mapping, reaction arrow and plus.                                                                                                                                                       |
| 6   | Monomers: monomer placement, fragments, connections, and the custom monomer tool.                                                                                                                                                                    |

**Done when:** every scene tool runs in the core, the Qt app behaves identically, and the tool tests build against `sketcher_core` without Qt.

### Phase 4 — Lean WASM core and protocol

Build the Qt-free WASM target and the command/event protocol that every host will use.

| PR  | Change                                                                                                                                                                                                                                                  |
| --- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | Emscripten target for `sketcher_core` + `rdkit_extensions`, built in CI, with the bundle size tracked.                                                                                                                                                  |
| 2   | Protocol spec (JSON commands and events, plus the display-list diff stream), a core-side command dispatcher with C++ tests, and codegen for the TS types.                                                                                               |
| 3   | Embind bindings for the protocol plus headless SVG, and a minimal demo page that edits a molecule through the controller.                                                                                                                               |
| 4   | Painter prototypes: SVG DOM and Canvas2D painters in the demo, with a large-molecule benchmark (for example a 2,000-residue peptide) against a 16 ms frame budget for drags and full redraws. Picks the painter and whether diffs need a binary format. |

**Done when:** the demo page loads and edits a molecule through the protocol, the painter choice is made from the benchmark, and the bundle size is measured against the Qt build.

### Phase 5 — React app (web first)

Build the new frontend in the repo at full feature parity with the Qt app. It replaces the Qt WASM build as "v2 web" only once nothing is missing. Since interaction already lives in the controller, this phase is view work: toolbars, menus, popups and dialogs. The UI refresh design (toolbars, menus, dialogs and depiction styling) is reviewed before the feature PRs start, and Crayon is evaluated as the component library before PR 2.

**Rules for every feature PR:**

- It ports the matching Playwright scenarios.
- It adds Vitest coverage for any new components.
- It's reachable in the dev build only, so the Qt WASM build stays the shipped one.
- State lives in the components that use it. A small context holds anything that's truly app-wide, and the core stays the source of truth. No central store holding all UI state.

PRs 1–2 land in order. PRs 3–9 can land in any order, or in parallel, after that. PRs 10–11 come last.

| PR  | Change                                                                                                                                                                                                                                                                                                                                                              | Ports                                                                                                                                                      |
| --- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | Scaffolding: a Vite + React + TS app in the repo, with linting, Vitest and a CI build that produces a self-contained bundle. Includes the protocol client with its `Transport` interface and the embind transport, plus the canvas view: the chosen painter applying display-list diffs, pan/zoom/fit, and input forwarding to the controller.                      |                                                                                                                                                            |
| 2   | App shell: display settings and persisted preferences; undo/redo wired to the core stack; top and side bars built from the action registry with the refresh design's theme tokens; the atomistic/monomeric mode switch. Also re-points the Playwright page-object wrappers at the React DOM (`data-testid`) and canvas, with screenshot baselines for the new look. | `undo_redo`, `wrapper_contracts`, `toolbar`                                                                                                                |
| 3   | Select UI: select-options widget and selection popup.                                                                                                                                                                                                                                                                                                               | `tst_select_mode_*`, `tst_move_mode`, `tst_move_select_parity`, `tst_erase_mode`                                                                           |
| 4   | Drawing UI: element buttons, periodic table, atom-query popup, set-atom widget, and the bond-order, stereo and query-bond popups.                                                                                                                                                                                                                                   | `drawing`, `tst_tools`, `tst_tools_parity`, `tst_editing_parity`, `tst_query_bond_crash`, `context_set_element_probe`, `context_wildcard_probe`            |
| 5   | Enumeration UI: R-group and attachment-point controls, atom mapping, reaction arrow and plus.                                                                                                                                                                                                                                                                       | `tst_enumeration_details`, `context_rgroup_probe`                                                                                                          |
| 6   | Monomer UI: amino acid and nucleic acid palettes, nucleotide and custom-nucleotide popups, connection options, the custom monomer dialog, and loading a monomer database.                                                                                                                                                                                           | `monomer_mutation`                                                                                                                                         |
| 7   | Top-bar menus and the dialogs they open: import (file, paste in text, replace content), export and save image, More Actions (Modify All, Expand Selection, add custom fragment), view toggles and preferences/rendering settings, and help, getting started and about.                                                                                              | `tst_import_menu`, `import`, `tst_export_menu`, `tst_more_actions_menu`, `tst_configure_view_menu`, `tst_help_menu`                                        |
| 8   | Context menus rendered from the action registry (atom, bond, selection, background, bracket subgroup, attachment point, monomer), clipboard (cut, copy, paste, Copy As formats and image) and keyboard shortcuts, including the hidden ones.                                                                                                                        | `tst_bond_context_menu`, `context_modify_atoms_probe`, `context_modify_bonds_probe`, `context_copy_probe`, `tst_copy_all_as_image`, `tst_hidden_shortcuts` |
| 9   | Structure dialogs: edit atom properties and bracket subgroup.                                                                                                                                                                                                                                                                                                       | `tst_edit_atom_properties`, `edit_atom_properties_probe`, `tst_bracket_subgroup`                                                                           |
| 10  | Embedding JS API (today's public WASM functions, so web embedders don't change) and a parity audit: check every Qt action, tool and dialog against the React app, close any gaps, and port the remaining Playwright files.                                                                                                                                          | `wasm_api`, remaining files                                                                                                                                |
| 11  | Replace the Qt WASM build with v2 web.                                                                                                                                                                                                                                                                                                                              |                                                                                                                                                            |

**Done when:** every tool, dialog and menu action in the Qt app exists in v2 web, and the full ported Playwright suite passes.

### Phase 6 — Qt-host integration

Swap `SketcherWidget` to a `QWebEngineView` shell around the same React bundle. The public C++ and SIP API stays the same, and the core is held in-process, so synchronous calls stay synchronous.

**Gate: QtWebEngine in Qt hosts.** Answer this during Phases 2–4, well before this phase starts. It takes a time-boxed prototype on a branch, and nothing from it merges. The output is a written go/no-go decision on running `QWebEngineView` inside Qt hosts, covering:

- **Packaging:** whether our Qt package can include QtWebEngine. It's currently built with `-skip qtwebengine`, so enabling it means new build deps and a Chromium build.
- **Distribution:** what QtWebEngine adds to Maestro's size and process model.
- **Startup ordering:** QtWebEngine must be initialized before the host creates its `QApplication`, which means changes on the host side.
- **Scale:** startup latency and memory with several sketchers in one PyQt panel.
- **Integration:** focus, shortcuts, clipboard, drag-and-drop and theming across the webview boundary.
- **SIP:** synchronous SIP-exposed calls behave unchanged.

**If the gate fails:** Qt hosts get a thin native Qt frontend over the controller instead. Toolbars and menus come from the action registry and the canvas paints the display list, so only the dialogs are duplicated in Qt. That gives up the single frontend for Qt hosts and needs an explicit decision.

| PR  | Change                                                                                                                                                                      |
| --- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | Enable QtWebEngine in the package-factory Qt recipe.                                                                                                                        |
| 2   | QWebChannel transport adapter, plus a web-backed widget built next to the existing one. The public-contract tests run against both.                                         |
| 3   | `SketcherView`: a native, canvas-only, select-only widget using the controller's select tool and a QPainter painter. A matching mmshare PR moves the Maestro overlay to it. |
| 4   | Swap `SketcherWidget` to the web-backed implementation and rebuild SIP. A preview build goes to the Maestro and Python panel owners before this merges.                     |

**Done when:** the sketcher C++ tests and the mmshare Maestro and Python tests pass, and consumers need no API changes apart from the overlay switch.

### Phase 7 — Qt-free standalone app

| PR  | Change                                                                                                                                                   |
| --- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | Vendor `webview/webview`, either in the repo or as a feedstock, and add the platform packaging deps: WebKitGTK on Linux and the WebView2 SDK on Windows. |
| 2   | New `sketcher_standalone` target: a `main.cpp` that hosts the React bundle in a system webview and links only `sketcher_core` + `rdkit_extensions`.      |
| 3   | Switch the shipped app to `sketcher_standalone`.                                                                                                         |

### Phase 8 — Delete

One PR per area:

- `widget/`, `dialog/`, `menu/` and the Qt scene adapter from Phase 3.
- The molviewer QGraphicsItem code.
- The Qt WASM build and the Playwright Qt bridge.
- Any Qt test infrastructure that's no longer needed.

**End state:** a Qt-free standalone app, a lean Qt-free WASM bundle, and a `sketcher` library whose only Qt is the webview shim and the overlay view.

## 6. Tech stack

| Concern             | Pick                                                | Why                                                                                        |
| ------------------- | --------------------------------------------------- | ------------------------------------------------------------------------------------------ |
| UI framework        | React + TypeScript                                  | Ecosystem and hiring. Its runtime size is negligible next to the WASM core.                |
| Component library   | Crayon (to evaluate)                                | Already used by LiveDesign, so the two look consistent and share design work               |
| Interactive drawing | SVG DOM or Canvas2D painting the C++ display list   | Chosen by the Phase 4 benchmark. WebGL is overkill.                                        |
| Headless rendering  | C++ SVG writer                                      | Works with no browser or display server, from the same display list as the UI              |
| Fonts               | Arimo everywhere; metrics via `stb_truetype`        | C++ layout and browser text rendering agree                                                |
| UI state            | Component-local state; a small context for app-wide | The C++ core is the source of truth, including the active tool, so UI state is small       |
| Core notifications  | `boost::signals2`                                   | Already in the build via `rdkit_extensions`; `scoped_connection` disconnects on scope exit |
| Build / unit tests  | Vite / Vitest                                       | Self-contained bundle the webviews can load from `file://`                                 |
| E2E / visual        | Playwright (extend `test/wasm/`)                    | Reuses the existing scenarios, and covers every webview engine                             |
| WASM bindings       | embind                                              | Already on emscripten                                                                      |
| TS ⇄ C++ types      | Small Python codegen from a JSON spec               | Lighter than protobuf                                                                      |
| Standalone host     | `webview/webview`                                   | Pure C++ on the system webview. CEF is the fallback.                                       |

## 7. Risks

- **QtWebEngine may be unacceptable in Qt hosts.**
  - _Why:_ a new Maestro dependency with startup and memory costs, a Qt recipe change, and `QWebEngineView` doesn't work inside `QGraphicsProxyWidget`.
  - _Mitigation:_ the Qt-host gate in Phase 6, a native `SketcherView` for the overlay, and the fallback to a thin native Qt frontend over the controller.
- **Depiction correctness during extraction** (wedges, stereo annotations, label placement, monomer shapes).
  - _Mitigation:_ Phase 2 is a faithful port checked by the existing screenshots. After the refresh, golden images of the new look take over as the regression net.
- **Interaction regressions while porting tools.**
  - _Mitigation:_ Phase 3 is a faithful port with today's Qt app as the consumer, checked by the existing Playwright suite, before any React work depends on it.
- **Display-list throughput on large molecules** (big peptides, macrocycles).
  - _Mitigation:_ retained, per-item diffs so a drag sends only what moved; the Phase 4 benchmark sets the painter and the wire format before the React app is built on them.
- **Text metrics drift between C++ and browsers.**
  - _Mitigation:_ one bundled font with metrics read from the same TTF, plus Playwright screenshots.
- **Undo across the JS/WASM boundary.**
  - _Mitigation:_ the core owns the undo stack, and every edit is an atomic command. React keeps no undo of its own.
- **Native features from inside a webview** (file dialogs, clipboard, drag-and-drop) are a common source of pain in webview apps.
  - _Mitigation:_ use the web APIs wherever they're enough, and keep the remaining native calls few and behind the `Transport`.
- **UI breadth is the long pole:** 36 tools and 30 dialogs.
  - _Mitigation:_ the tools are ported in C++ rather than rewritten, so the React work is views only. The Qt app keeps shipping until the React app reaches parity, so no partial release is ever needed. The ported Playwright suite tracks how much is left.
- **System webviews render differently on each platform** (affects the standalone app).
  - _Mitigation:_ the bundled font, Playwright runs on all three engines, and CEF as the fallback.
