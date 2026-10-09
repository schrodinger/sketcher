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
- Keep headless image generation, with no browser or Qt needed.

Matching today's Qt look is not required. The rewrite is a chance to refresh the UI if we want one.

## 2. Key decisions

### 2.1 Chemistry stays in C++ (no RDKit.js / MinimalLib)

- **Consumers need C++ RDKit objects, in-process and synchronously.** The public `SketcherWidget` API passes them directly: `getRDKitMolecule()` returns an `RDKit::ROMol`, selection is a set of `RDKit::Atom*`/`Bond*`, and `atomHovered(const RDKit::Atom*)` is a signal. MinimalLib can't provide any of that.
- **MinimalLib's API is too narrow.** It covers whole-molecule I/O, depiction and descriptors. The sketcher edits molecules atom by atom and bond by bond, and relies on substance groups, stereo handling, and HELM, FASTA and the monomer database in `rdkit_extensions`. All of these are built on C++ APIs that MinimalLib doesn't expose.
- **One implementation.** Desktop, Qt hosts and headless rendering need the C++ chemistry anyway. Moving to RDKit.js for the web would mean keeping two chemistry layers in sync.

### 2.2 Webview + React over Slint

A single UI on every target means either a web UI in a webview everywhere, or a native cross-platform toolkit like Slint.

- **Why webview + React wins:**
  - **Ecosystem:** mature component libraries for the chrome (palettes, dialogs, HELM input, monomer browsers).
  - **Speed and hiring:** iteration speed and how easy it is to hire for.
  - **Bundle size:** not a real differentiator, because the RDKit WASM core (about 10–15 MB) dwarfs either UI runtime.
- **What Slint does better:** a smaller desktop footprint, and a cleaner way to embed in Qt hosts.
- **Gate:** whether a webview is acceptable inside Qt hosts is decided by the Qt-host gate in Phase 5.

### 2.3 Depiction geometry lives in C++

A new Qt-free `depict/` module computes the geometry: label layout, bond trimming, double-bond offsets, wedges, monomer shapes and hit-test shapes. It emits a **display list** of lines, polygons, paths and positioned text. Text is measured through a `FontMetrics` interface backed by the bundled Arimo TTF, so every target lays out text the same way.

Thin painters turn the display list into output:

| Painter             | Used for                                |
| ------------------- | --------------------------------------- |
| Canvas2D (TS)       | The interactive UI on every target      |
| SVG writer (C++)    | Headless image generation               |
| QPainter (Qt, thin) | The canvas-only overlay view in Maestro |

Keeping geometry in C++ means headless rendering and every UI share a single source. Putting it in TS would leave no browser-free renderer for the image API.

### 2.4 Alternatives rejected

- **Electron:** its main process is Node, and our core is C++. Electron would need a native module or a child process to reach the core. With `webview/webview`, the C++ binary _is_ the host.
- **CEF:** bundles Chromium, so rendering is identical on every platform, at the cost of about 100 MB. It stays the fallback if the system webviews turn out to render the canvas differently from each other.

## 3. Starting point

| Layer                             | LOC   | Qt coupling                             | Fate                                      |
| --------------------------------- | ----- | --------------------------------------- | ----------------------------------------- |
| `rdkit_extensions/`               | 16.3K | None                                    | Keep                                      |
| `sketcher/rdkit/`                 | 7.9K  | Minimal                                 | Move to `sketcher_core`                   |
| `sketcher/model/`                 | 7.5K  | `QObject`, `QUndoStack`, a few Qt types | Port to `sketcher_core`                   |
| `sketcher/molviewer/` (rendering) | 11.6K | QPainter/QGraphicsItem throughout       | Geometry → `depict/`; the rest is deleted |
| `sketcher/tool/` (36 scene tools) | 8.0K  | `QObject`                               | Rewrite in TS                             |
| `widget/` + `dialog/` + `menu/`   | 11.3K | All QWidget                             | Rewrite in React                          |

- **WASM ships the whole Qt app today:** 44 MB raw, 13 MB gzipped. A lean core should land around 10–15 MB raw, about 3× smaller.
- **Regression net (already in place):**
  - The sketcher C++ tests, including the public-contract tests for the widget and image generation.
  - The Maestro and Python consumer tests in mmshare.
  - The Playwright suite: 41 files and about 147 screenshot checkpoints, written against page-object wrappers.

## 4. Target architecture

Each consumer embeds the same React bundle, and React talks to `sketcher_core`. The one exception is the Maestro overlay, which is select-only, has no chrome, and calls the core directly.

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
│   Canvas2D painter · tools · chrome (bars, dialogs)    │ │ no React           │
│   protocol client (same code on every host)            │ │                    │
└───────┬──────────────────┬───────────────────┬─────────┘ └───────────┬────────┘
┌───────▼────────┐ ┌───────▼────────┐ ┌────────▼─────────┐             │
│ embind         │ │ bind / eval    │ │ QWebChannel      │             │
└───────┬────────┘ └───────┬────────┘ └────────┬─────────┘             │
        └──────────────────┼───────────────────┘                       │
              JSON commands in / events out                  direct C++ calls
┌──────────────────────────▼───────────────────────────────────────────▼────────┐
│ sketcher_core (new, no Qt)                                                    │
│   protocol/ command dispatcher (same code on every host)                      │
│   model/    MolModel, edit commands, undo stack, selection, events            │
│   depict/   display list, FontMetrics (Arimo), headless SVG                   │
│   rdkit/    sketcher chemistry helpers                                        │
└──────────────────────────────────────┬────────────────────────────────────────┘
┌──────────────────────────────────────▼────────────────────────────────────────┐
│ rdkit_extensions (unchanged): HELM, FASTA, SMILES, MOL, monomer DB            │
└──────────────────────────────────────┬────────────────────────────────────────┘
                                  RDKit (C++)
```

**Libraries:** there are three.

- `rdkit_extensions` stays as it is.
- `sketcher_core` has no Qt and is linked by every host.
- `sketcher` becomes a thin Qt shim containing `SketcherWidget` and `SketcherView`.

**Protocol:** there is one protocol, with one implementation on each side.

- **C++ side:** a command dispatcher in `sketcher_core`.
- **TS side:** a protocol client in React, written against a small `Transport` interface (send a command, subscribe to events).
- **Per host:** each transport adapter is a few dozen lines on each side.
  - embind is synchronous; `bind`/`eval` and QWebChannel are asynchronous.
  - The client treats every call as async, so React code is the same on every host.
- **Qt hosts' synchronous calls:** calls such as `getRDKitMolecule()` go straight from the shim to the core, not through React.

**Split of responsibilities:**

- **C++ owns** chemistry, edit semantics, undo, selection and depiction geometry.
- **TS owns** interaction (mouse handling, drag state, tool modes), the chrome, and painting the display list.
- **Crossing the boundary:** only committed edits cross into the core. Drag previews are drawn in TS, so nothing crosses per mouse move.

## 5. Phased plan

Every PR lands on `main` and leaves the existing Qt app working until Phase 7. Phases 1–4 pay off whatever the Qt-host gate in Phase 5 decides.

### Phase 1 — Extract `sketcher_core` (model + rdkit)

Make the chemistry and model layer Qt-free behind a new library, without changing behavior. `SketcherModel` stays in the Qt library: it holds UI state (current tool, display settings) that moves to React in Phase 4.

| PR  | Change                                                                                                                                                                                  |
| --- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | Untangle dependencies: move the chemistry enums, coordinate helpers, monomer utilities and constants out of Qt headers, and move the arrowhead-placement depiction code into molviewer. |
| 2   | Remove incidental Qt types from `MolModel`: `QString` becomes `std::string`, `QColor` becomes an RGBA struct, `Q_ASSERT` becomes `assert`. Also remove `QString` from `rdkit/`.         |
| 3   | Replace `MolModel`'s Qt signals with an observer registry that disconnects automatically (RAII). A thin QObject adapter re-emits the signals so existing consumers don't change.        |
| 4   | Replace `QUndoStack`/`QUndoCommand` with a small core undo stack that supports macros, merging and can-undo/redo notifications.                                                         |
| 5   | Create the `sketcher_core` CMake target with no Qt link, move the files, and add the export macro, install rules and recipe `install_name_tool` lines.                                  |
| 6   | Move the model and rdkit tests onto `sketcher_core` alone, which proves the core is Qt-free.                                                                                            |

**Done when:** the app behaves identically, all existing tests pass, and the core tests build without Qt.

### Phase 2 — Extract depiction into `depict/`

Move the geometry out of the molviewer `paint()` methods into the display-list producer. The existing QGraphicsItems become thin painters of the display list, so today's Qt app is the first consumer. **Port faithfully, with no visual changes.** The existing Playwright screenshots must pass unchanged after every PR.

| PR  | Change                                                                                                           |
| --- | ---------------------------------------------------------------------------------------------------------------- |
| 1   | Display-list types, the `FontMetrics` interface and the Arimo implementation, with unit tests. No consumers yet. |
| 2   | Atom labels: layout and hit shapes move to `depict/`, and the atom item paints the display list.                 |
| 3   | Bonds: trimming, double-bond offsets, wedges and hashes.                                                         |
| 4   | Monomers: shapes, connectors and attachment-point labels.                                                        |
| 5   | Everything else: arrows, pluses, S-group brackets, and selection and hover highlights.                           |
| 6   | SVG writer and a Qt-free SVG image API, added next to the existing Qt-typed functions.                           |

**Done when:** every visual item draws from the display list, headless SVG works without Qt, and the screenshots are unchanged.

### Phase 3 — Lean WASM core and protocol

Build the Qt-free WASM target and the command/event protocol that every host will use.

| PR  | Change                                                                                                                       |
| --- | ---------------------------------------------------------------------------------------------------------------------------- |
| 1   | Emscripten target for `sketcher_core` + `rdkit_extensions`, built in CI, with the bundle size tracked.                       |
| 2   | Protocol spec (JSON commands and events), a core-side command dispatcher with C++ tests, and codegen for the TS types.       |
| 3   | Embind bindings for the protocol plus headless SVG, and a minimal demo page that renders a molecule with a Canvas2D painter. |

**Done when:** the demo page loads and edits a molecule through the protocol, and the bundle size is measured against the Qt build.

### Phase 4 — React app (web first)

Build the new frontend in the repo at full feature parity with the Qt app. It replaces the Qt WASM build as "v2 web" only once nothing is missing. The UI refresh design (chrome and depiction styling) is reviewed before the feature PRs start.

**Rules for every feature PR:**

- It ports the matching Playwright scenarios.
- It adds Vitest coverage for any new components.
- It's reachable in the dev build only, so the Qt WASM build stays the shipped one.

PRs 1–2 land in order. PRs 3–9 can land in any order, or in parallel, after that. PRs 10–11 come last.

| PR  | Change                                                                                                                                                                                                                                                                                                                                                                                                            | Ports                                                                                                                                                      |
| --- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | Scaffolding: a Vite + React + TS app in the repo, with linting, Vitest and a CI build that produces a self-contained bundle. Includes the protocol client with its `Transport` interface and the embind transport, plus the Canvas2D display-list view: pan/zoom/fit, hit-testing against the display-list hit shapes, and hover highlighting.                                                                    |                                                                                                                                                            |
| 2   | App shell: a Zustand store replacing `SketcherModel`'s UI state (active tool, display settings, persisted preferences); undo/redo wired to the core stack; top and side bars with the refresh design's theme tokens; the atomistic/monomeric mode switch. Also re-points the Playwright page-object wrappers at the React DOM (`data-testid`) and canvas coordinates, with screenshot baselines for the new look. | `undo_redo`, `wrapper_contracts`, `toolbar`                                                                                                                |
| 3   | Select (rectangle, lasso and fragment modes, select-options widget, selection popup), move/rotate, and erase.                                                                                                                                                                                                                                                                                                     | `tst_select_mode_*`, `tst_move_mode`, `tst_move_select_parity`, `tst_erase_mode`                                                                           |
| 4   | Drawing: atoms (element buttons, periodic table, atom-query popup, set-atom widget), bonds (bond-order, stereo and query-bond popups, chains), rings and custom fragments, charge and explicit H.                                                                                                                                                                                                                 | `drawing`, `tst_tools`, `tst_tools_parity`, `tst_editing_parity`, `tst_query_bond_crash`, `context_set_element_probe`, `context_wildcard_probe`            |
| 5   | Enumeration: new and existing R-groups, attachment points, atom mapping, reaction arrow and plus.                                                                                                                                                                                                                                                                                                                 | `tst_enumeration_details`, `context_rgroup_probe`                                                                                                          |
| 6   | Monomers: amino acid and nucleic acid palettes, nucleotide and custom-nucleotide popups, monomer fragments, covalent/disulfide and H-bond connections, the custom monomer tool and dialog, and loading a monomer database.                                                                                                                                                                                        | `monomer_mutation`                                                                                                                                         |
| 7   | Top-bar menus and the dialogs they open: import (file, paste in text, replace content), export and save image, More Actions (Modify All, Expand Selection, add custom fragment), view toggles and preferences/rendering settings, and help, getting started and about.                                                                                                                                            | `tst_import_menu`, `import`, `tst_export_menu`, `tst_more_actions_menu`, `tst_configure_view_menu`, `tst_help_menu`                                        |
| 8   | Context menus (atom, bond, selection, background, bracket subgroup, attachment point, monomer), clipboard (cut, copy, paste, Copy As formats and image) and keyboard shortcuts, including the hidden ones.                                                                                                                                                                                                        | `tst_bond_context_menu`, `context_modify_atoms_probe`, `context_modify_bonds_probe`, `context_copy_probe`, `tst_copy_all_as_image`, `tst_hidden_shortcuts` |
| 9   | Structure dialogs: edit atom properties and bracket subgroup.                                                                                                                                                                                                                                                                                                                                                     | `tst_edit_atom_properties`, `edit_atom_properties_probe`, `tst_bracket_subgroup`                                                                           |
| 10  | Embedding JS API (today's public WASM functions, so web embedders don't change) and a parity audit: check every Qt action, tool and dialog against the React app, close any gaps, and port the remaining Playwright files.                                                                                                                                                                                        | `wasm_api`, remaining files                                                                                                                                |
| 11  | Replace the Qt WASM build with v2 web.                                                                                                                                                                                                                                                                                                                                                                            |                                                                                                                                                            |

**Done when:** every tool, dialog and menu action in the Qt app exists in v2 web, and the full ported Playwright suite passes.

### Phase 5 — Qt-host integration

Swap `SketcherWidget` to a `QWebEngineView` shell around the same React bundle. The public C++ and SIP API stays the same, and the core is held in-process, so synchronous calls stay synchronous.

**Gate: QtWebEngine in Qt hosts.** Answer this during Phases 2–3, well before this phase starts. It takes a time-boxed prototype on a branch, and nothing from it merges. The output is a written go/no-go decision on running `QWebEngineView` inside Qt hosts, covering:

- **Packaging:** whether our Qt package can include QtWebEngine. It's currently built with `-skip qtwebengine`, so enabling it means new build deps and a Chromium build.
- **Distribution:** what QtWebEngine adds to Maestro's size and process model.
- **Startup ordering:** QtWebEngine must be initialized before the host creates its `QApplication`, which means changes on the host side.
- **Scale:** startup latency and memory with several sketchers in one PyQt panel.
- **Integration:** focus, shortcuts, clipboard, drag-and-drop and theming across the webview boundary.
- **SIP:** synchronous SIP-exposed calls behave unchanged.

**If the gate fails:** this phase is replaced by a native Qt UI over `sketcher_core` for Qt hosts. That gives up the single frontend for Qt hosts and needs an explicit decision.

| PR  | Change                                                                                                                                                  |
| --- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | Enable QtWebEngine in the package-factory Qt recipe.                                                                                                    |
| 2   | QWebChannel transport adapter, plus a web-backed widget built next to the existing one. The public-contract tests run against both.                     |
| 3   | `SketcherView`: a native, canvas-only, select-only widget with a QPainter painter. A matching mmshare PR moves the Maestro overlay to it.               |
| 4   | Swap `SketcherWidget` to the web-backed implementation and rebuild SIP. A preview build goes to the Maestro and Python panel owners before this merges. |

**Done when:** the sketcher C++ tests and the mmshare Maestro and Python tests pass, and consumers need no API changes apart from the overlay switch.

### Phase 6 — Qt-free standalone app

| PR  | Change                                                                                                                                                   |
| --- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | Vendor `webview/webview`, either in the repo or as a feedstock, and add the platform packaging deps: WebKitGTK on Linux and the WebView2 SDK on Windows. |
| 2   | New `sketcher_standalone` target: a `main.cpp` that hosts the React bundle in a system webview and links only `sketcher_core` + `rdkit_extensions`.      |
| 3   | Switch the shipped app to `sketcher_standalone`.                                                                                                         |

### Phase 7 — Delete

One PR per area:

- `tool/`, `widget/`, `dialog/` and `menu/`.
- The molviewer QGraphicsItem code.
- The Qt WASM build and the Playwright Qt bridge.
- Any Qt test infrastructure that's no longer needed.

**End state:** a Qt-free standalone app, a lean Qt-free WASM bundle, and a `sketcher` library whose only Qt is the webview shim and the overlay view.

## 6. Tech stack

| Concern             | Pick                                         | Why                                                                         |
| ------------------- | -------------------------------------------- | --------------------------------------------------------------------------- |
| UI framework        | React + TypeScript                           | Ecosystem and hiring. Its runtime size is negligible next to the WASM core. |
| Interactive drawing | Canvas2D painting the C++ display list       | Fast enough for hundreds of atoms. WebGL is overkill.                       |
| Headless rendering  | C++ SVG writer                               | Works with no browser, from the same display list as the UI                 |
| Fonts               | Arimo everywhere; metrics via `stb_truetype` | C++ layout and browser text rendering agree                                 |
| UI state            | Zustand                                      | The C++ core is the source of truth, so UI state is small                   |
| Build / unit tests  | Vite / Vitest                                | Self-contained bundle the webviews can load from `file://`                  |
| E2E / visual        | Playwright (extend `test/wasm/`)             | Reuses the existing scenarios, and covers every webview engine              |
| WASM bindings       | embind                                       | Already on emscripten                                                       |
| TS ⇄ C++ types      | Small Python codegen from a JSON spec        | Lighter than protobuf                                                       |
| Standalone host     | `webview/webview`                            | Pure C++ on the system webview. CEF is the fallback.                        |

## 7. Risks

- **QtWebEngine may be unacceptable in Qt hosts.**
  - _Why:_ a new Maestro dependency with startup and memory costs, a Qt recipe change, and `QWebEngineView` doesn't work inside `QGraphicsProxyWidget`.
  - _Mitigation:_ the Qt-host gate in Phase 5, a native `SketcherView` for the overlay, and the fallback to a native Qt UI over the core.
- **Depiction correctness during extraction** (wedges, stereo annotations, label placement, monomer shapes).
  - _Mitigation:_ Phase 2 is a faithful port checked by the existing screenshots. After the refresh, golden images of the new look take over as the regression net.
- **Text metrics drift between C++ and browsers.**
  - _Mitigation:_ one bundled font with metrics read from the same TTF, plus Playwright screenshots.
- **Undo across the JS/WASM boundary.**
  - _Mitigation:_ the core owns the undo stack, and every edit is an atomic command. React keeps no undo of its own.
- **UI breadth is the long pole:** 36 tools and 30 dialogs.
  - _Mitigation:_ the Qt app keeps shipping until the React app reaches parity, so no partial release is ever needed. The ported Playwright suite tracks how much is left.
- **System webviews render differently on each platform** (affects the standalone app).
  - _Mitigation:_ the bundled font, Playwright runs on all three engines, and CEF as the fallback.
