# Sketcher: Qt Removal & React Migration Plan

Drafted 2026-05-16 · Revised 2026-10-08 · Repo: `schrodinger/sketcher` · Owner: Chris Von Bargen

> **Revision 2026-10-08.** Corrected the RDKit constraint (the sketcher _consumes_ RDKit). Added the consumers the first draft missed: PyQt panels via SIP, the headless image API, and Maestro's `QGraphicsScene` overlay. Moved depiction geometry into C++ (a `depict/` module in a new Qt-free `sketcher_core` library) so headless rendering and every UI share one geometry source. Added a go/no-go spike before committing Qt hosts to a webview. Re-estimated effort. The legacy sketcher copy in `mmshare/src/schrodinger/sketcher/` is out of scope; it is being removed in favor of mmshare consuming this repo's build.
>
> **Decisions 2026-10-08.** The public image API becomes a Qt-free SVG API; Qt wrappers will live in mmshare, and retiring the Qt-typed functions is managed separately. The public contract is already covered by existing tests (sketcher C++ tests plus Maestro/Python tests in mmshare), so Phase 0 is done. Matching today's look is not a goal: this is also a UI refresh. One bundled font everywhere: Arimo, which the sketcher already ships.

- [1. Goals & constraints](#1-goals--constraints)
- [2. Key decisions (with rationale)](#2-key-decisions-with-rationale)
- [3. Current state survey](#3-current-state-survey)
- [4. Target architecture](#4-target-architecture)
- [5. Keep / port / replace matrix](#5-keep--port--replace-matrix)
- [6. Phased transition plan](#6-phased-transition-plan)
- [7. Tech stack picks](#7-tech-stack-picks)
- [8. Risks & mitigations](#8-risks--mitigations)
- [9. Rough total effort](#9-rough-total-effort)
- [10. Open questions to resolve before starting](#10-open-questions-to-resolve-before-starting)

## 1. Goals & constraints

Remove the Qt dependency from the sketcher library, subject to:

- **The sketcher consumes RDKit through its standard C++ APIs.** Chemistry, edit semantics and undo stay in C++. No move to RDKit.js / MinimalLib.
- **Qt host consumers keep working, linking this repo's build via mmshare:**
  - Maestro (C++). Embeds `SketcherWidget` in dialogs, and embeds it in the 3D workspace's `QGraphicsScene` as the ligand overlay (`scene()->addWidget()`, select-only mode).
  - ~105 Python files using `schrodinger.ui.sketcher`: PyQt panels that embed or subclass the SIP-wrapped `SketcherWidget`.
- **Headless image generation keeps working without a browser, through a Qt-free SVG API.** The Qt-typed `get_qpicture` / `get_qimage` / `save_image_file` stay for now; their retirement is managed separately. Today's callers (Maestro, Python's `schrodinger/ui/sketcher_image.py`) move to a Qt wrapper over the SVG API that lives in mmshare. The WASM build (`image_generation_from_js`) uses the SVG API directly.
- **The library continues to produce a drop-into-a-webpage WASM build**, and that build should be **smaller and more performant** than today's Qt WASM bundle.
- **A standalone desktop app** (`src/app/main.cpp`) continues to ship.

UX requirement: a **single unified frontend** that looks and behaves identically on every target that shows a UI. This is also a **UI refresh**: matching today's Qt look is not a goal.

## 2. Key decisions (with rationale)

### 2.1 Webview-everywhere over Slint

The unified-FE constraint rules out the "core library + per-platform UI" split. The two viable single-UI strategies are: a web UI hosted in a webview on every target, or a native cross-platform UI framework like Slint (Rust, with C++ bindings and WASM target).

| Dimension                                                              | Webview + React                                                                                                                       | Slint                                                                         |
| ---------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------- |
| Web bundle size                                                        | ~30–100 KB JS on top of the WASM core.                                                                                                | ~1–2 MB Slint runtime on top of the WASM core.                                |
| Desktop footprint                                                      | +100–200 MB Chromium when bundled (CEF/Electron/QtWebEngine); ~50–100 MB with native system webview.                                  | **Wins.** A few MB native; direct C++ calls.                                  |
| Ecosystem for chrome (palettes, dialogs, HELM input, monomer browsers) | **Wins decisively.** Mature component libs for everything.                                                                            | Small widget set; custom work on you.                                         |
| Embedding in the Qt host                                               | QWebEngineView + QWebChannel (async, somewhat clunky). New runtime dependency for Maestro; nothing in mmshare uses QtWebEngine today. | slint-cpp into a QWidget (sync, cleaner). C++ bindings less mature than Rust. |
| Iteration speed & hiring                                               | **Wins.** Hot reload, devtools, common skillset.                                                                                      | Smaller community.                                                            |
| "Looks identical everywhere"                                           | Same Chromium engine (when bundled); some drift with system webviews.                                                                 | Slint owns the renderer end-to-end — strongest guarantee.                     |

**Decision:** webview + React, gated by the Qt-host spike in Phase 4. The deciding factor is the ecosystem advantage for the chrome. Bundle size is _not_ a real differentiator: the RDKit-based WASM core (likely 10–15 MB raw) dwarfs either UI runtime.

### 2.2 React as the frontend framework

Default pick. Largest ecosystem, easiest hiring, runtime size (~45 KB gzipped) is rounding error next to the WASM chem core. The usual re-render complaints don't bite here because the hot interactive loop lives in Canvas2D outside React's reconciler.

> Solid (~7 KB runtime, React-like JSX) or Svelte are reasonable alternatives if bundle size becomes a tier-1 metric. Both give up ecosystem depth for the chrome.

### 2.3 Depiction geometry lives in C++, not TS

Atom-label layout, bond trimming around labels, double-bond offsets, wedge/hash geometry, monomer shapes and hit-test shapes are computed by a new Qt-free C++ module, `depict/`, inside `sketcher_core`. It emits a **display list** of primitives: lines, polygons, paths and positioned text runs. Text measurement goes through an injected `FontMetrics` interface. The default implementation reads metrics directly from the bundled Arimo TTF (e.g. header-only `stb_truetype`, which works in WASM), so every target lays out text the same way.

Several thin **painters** turn the display list into pixels:

- **Canvas2D (TS):** the interactive UI on every target.
- **SVG writer (C++, no Qt):** headless image generation in WASM and on servers.
- **QPainter painter (thin, Qt):** powers only the canvas-only `SketcherView` (see 2.4).

_Why:_ putting geometry in TS (the first draft) left no browser-free renderer for the image API. It would also have forced the layout logic to be duplicated, with drift, for headless output. One geometry source keeps every output consistent.

### 2.4 Hosts — four targets

| Build target                                               | UI                                                                                  | Transport to `sketcher_core`                                                                                                                                         |
| ---------------------------------------------------------- | ----------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Standalone desktop (`main.cpp` in this repo, **no Qt**)    | React in OS-native webview via `webview/webview` (WKWebView / WebView2 / WebKitGTK) | webview `bind` + `eval` → direct C++ calls                                                                                                                           |
| Qt hosts: Maestro dialogs + PyQt panels (`SketcherWidget`) | React in `QWebEngineView` (pending Phase 4 spike)                                   | QWebChannel → C++ adapter holding `sketcher_core` in-process. Synchronous API calls like `getRDKitMolecule()` stay synchronous because they never cross the channel. |
| Maestro ligand overlay (`SketcherView`)                    | Native canvas-only QWidget, display list via QPainter. No chrome.                   | Direct C++ calls                                                                                                                                                     |
| Web / WASM                                                 | React in browser                                                                    | embind + JS postMessage                                                                                                                                              |

Design _one_ protocol (JSON commands in, JSON events out) and ship thin transport adapters. The React bundle is byte-identical across all webview targets.

The overlay gets a native view because it is select-only (it always calls `activateSelectOnlyMode()`) and shows no toolbars, so there is no chrome to keep identical. It also sidesteps `QWebEngineView`'s lack of support inside `QGraphicsProxyWidget`.

### 2.5 Why not Electron

Disk/RAM weight of Electron (~150 MB) is in the same ballpark as the current Qt desktop build, so weight alone isn't the disqualifier. The real reason: **architectural mismatch**. Electron's main process is Node.js. `sketcher_core` is C++. Reaching the core from Electron means either spawning a C++ child process with IPC, or building an N-API native module — more glue, two JS contexts to reason about, no benefit. With `webview/webview` the C++ binary _is_ the host: it owns the webview, holds `sketcher_core` directly, and bridges with a few lines of C++.

Electron is defensible if its tooling/DX is what you want, and a Node wrapper around the C++ core is acceptable. Otherwise, no.

### 2.6 Why not CEF

CEF is the answer if cross-engine drift in system webviews proves unacceptable for canvas rendering — it bundles Chromium for guaranteed parity. Cost: ~100 MB on disk, partially undoing the lean-desktop win. Treat it as the fallback if Playwright visual regression flags real issues. Bundling a single font everywhere (Arimo) removes most of the drift risk up front.

## 3. Current state survey

#### Code composition (wc -l of .cpp/.h, Oct 2026)

| Layer                                               | LOC         | Qt coupling                                                                     |
| --------------------------------------------------- | ----------- | ------------------------------------------------------------------------------- |
| `src/schrodinger/rdkit_extensions/`                 | 16,343      | None                                                                            |
| `src/schrodinger/sketcher/rdkit/`                   | 7,856       | Minimal: 2 files use `QPointF`/`QLineF`/`QGraphicsItem`                         |
| `src/schrodinger/sketcher/molviewer/` (rendering)   | 11,615      | QPainter, QGraphicsItem throughout                                              |
| `src/schrodinger/sketcher/model/`                   | 7,521       | Light: `QObject`, `QUndoStack`/`QUndoCommand`, `QPointF` in a couple of headers |
| `src/schrodinger/sketcher/tool/` (36 scene tools)   | 8,045       | AbstractSceneTool : QObject                                                     |
| `src/schrodinger/sketcher/widget/` (30 .ui files)   | 5,642       | All QWidget                                                                     |
| `src/schrodinger/sketcher/dialog/`                  | 3,676       | All QDialog                                                                     |
| `src/schrodinger/sketcher/menu/`                    | 2,013       | All QMenu                                                                       |
| **Total sketcher library** (excl. rdkit_extensions) | **~49,300** |                                                                                 |

#### Key findings

- **RDKit integration is clean.** `RWMol` is wrapped, not subclassed. Conversions live in `rdkit_extensions`; sketcher-specific chemistry helpers in `sketcher/rdkit/` are nearly Qt-free.
- **Public API is two surfaces, both Qt-typed.**
  - `sketcher_widget.h`: `SketcherWidget : public QWidget`, ~60 methods. Signals include `atomHovered(const RDKit::Atom*)` and `bondHovered(const RDKit::Bond*)`.
  - `image_generation.h`: returns `QPicture`, `QImage` and `QByteArray`.
  - Plus `public_constants.h`, `image_constants.h` and `definitions.h`.
- **Consumers are broader than one host.**
  - Maestro: dialogs, `nonstandard_mutate_dlg`, and the ligand overlay, which subclasses `SketcherWidget`, installs it as a scene event filter and uses `get_qpicture`.
  - ~105 Python files via SIP.
  - The WASM app.
- **The model port is the easy part.** Qt in `model/` is mostly `QObject` + `QUndoStack`; the logic is already separable.
- **Rendering is entirely QPainter / QGraphicsScene.** No OpenGL, no custom rasterization. ~72 visual item `paint()` methods. The geometry inside them is the valuable part to port.
- **WASM today ships the entire Qt app.** Local build: `Sketcher.wasm` 44 MB raw / 13 MB gzipped. For comparison, `smiles_to_helm.wasm` (RDKit subset + Qt Core) is 12 MB / 3.75 MB gzipped. Expect a lean core to land around 10–15 MB raw: a real gain, roughly 3×, not 10×.
- **Test infra to reuse.**
  - **Playwright:** 41 test files (~4K lines) in `test/wasm/`, with ~147 screenshot checkpoints. Tests are written against page-object wrappers (`sk.click_tool('erase')`, `sk.import_menu(...)`), and the wrappers reach Qt through `src/app/playwright_test_bridge.cpp` (`objectName` lookups + canvas coordinates).
  - **Public contract:** `test_sketcher_widget.cpp` (53 cases) and `test_image_generation.cpp` here, plus Maestro and Python consumer tests in mmshare.
  - Most other C++ tests exercise chemistry/model logic without Qt widgets.

## 4. Target architecture

```
┌──────────────────────────────────────────────────────────────┐
│  React + TypeScript frontend                                 │
│  - Canvas2D painter (draws the display list)                 │
│  - Tool layer (36 interaction handlers, in TS)               │
│  - Chrome: top bar, palettes, dialogs, menus, popups         │
└──────────────────┬───────────────────────────────────────────┘
                   │  message-passing protocol (commands + events)
   ┌───────────────┼──────────────────────────────┐
┌──▼─────┐ ┌───────▼──────────────────────────┐ ┌─▼──────────────┐
│ WASM   │ │ libsketcher  (Qt shim)           │ │ standalone app │
│ embind │ │ SketcherWidget    SketcherView   │ │ webview/webview│
│        │ │ (QWebEngineView   (QPainter,     │ │ bind (no Qt)   │
│        │ │  + QWebChannel)    overlay only) │ │                │
└──┬─────┘ └───────┬──────────────────────────┘ └─┬──────────────┘
   └───────────────┼──────────────────────────────┘
┌──────────────────▼───────────────────────────────────────────┐
│  libsketcher_core  (NEW, no Qt)                              │
│  depict/   label layout, bond trimming, double-bond offsets, │
│            wedges, monomer shapes, hit-test shapes           │
│            display list out; FontMetrics in (Arimo TTF)      │
│            headless SVG painter → public image API           │
│  model/    MolModel, edit commands, undo stack, selection,   │
│            coordgen layout, callback-registry events         │
│  rdkit/    sketcher chemistry helpers (moved from sketcher/) │
└──────────────────┬───────────────────────────────────────────┘
┌──────────────────▼───────────────────────────────────────────┐
│  librdkit_extensions  (KEEP)                                 │
│  HELM, FASTA, SMILES, MOL conversions, monomer DB            │
└──────────────────┬───────────────────────────────────────────┘
                   ▼
             RDKit (C++ API)
```

### Libraries

Three link targets. The boxes inside `sketcher_core` are modules (directories), not separate libraries.

| Library               | Contents                                                                                    | Qt? | Linked by                                        |
| --------------------- | ------------------------------------------------------------------------------------------- | --- | ------------------------------------------------ |
| `rdkit_extensions`    | Unchanged; already a standalone public library                                              | No  | Everything, including Python/LiveDesign directly |
| `sketcher_core` (NEW) | `rdkit/` (moved from `sketcher/rdkit/`), `model/`, `depict/`                                | No  | WASM, standalone app, Qt shim                    |
| `sketcher` (Qt shim)  | `SketcherWidget` (QWebEngineView + QWebChannel adapter), `SketcherView` (QPainter, overlay) | Yes | Maestro, PyQt panels via SIP                     |

_Why not separate core and depict libraries:_ every consumer needs both, so a split adds a target without adding value. The one boundary that must be a separate library is Qt vs. no-Qt. The model-must-not-include-depict rule is enforced by directory convention (or a CI include check if it ever matters).

### Core split

C++ owns chemistry, edit semantics, undo, selection and depiction geometry. TS owns interaction (mouse events, drag state, tool modes), chrome, and painting the display list. Only commit-style operations (add bond, delete atom) cross into the core; the display list is recomputed once per commit. Drag previews are drawn by TS on top of the last display list, so there is no IPC per mouse move.

## 5. Keep / port / replace matrix

| Layer                                     | LOC            | Action                                                                                                                                                                         |
| ----------------------------------------- | -------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `rdkit_extensions/`                       | 16K            | **KEEP** Verbatim                                                                                                                                                              |
| `sketcher/rdkit/`                         | 8K             | **KEEP** Move into `sketcher_core` (`rdkit/`); replace the 2 Qt geometry usages with RDGeom types                                                                              |
| `model/`                                  | 7.5K           | **PORT** to `sketcher_core` (`model/`): strip QObject, callback registry instead of signals, custom undo stack instead of `QUndoStack`                                         |
| `molviewer/` geometry                     | ~6–8K of 11.6K | **PORT** to `sketcher_core` (`depict/`)                                                                                                                                        |
| `molviewer/` drawing + QGraphicsItem glue | remainder      | **REPLACE** with painters (Canvas2D, SVG, thin QPainter)                                                                                                                       |
| `image_generation`                        | —              | **PORT** Reimplement as a Qt-free SVG API on the display list. mmshare wraps SVG for Qt callers; Qt-typed functions retired later (managed separately).                        |
| `tool/`                                   | 8K             | **REPLACE** Reimplement in TS                                                                                                                                                  |
| `widget/` + `dialog/` + `menu/`           | 11K            | **REPLACE** Reimplement in React                                                                                                                                               |
| Chemistry/model tests                     | —              | **PORT** Re-target at `sketcher_core`                                                                                                                                          |
| Rendering tests                           | —              | **PORT** Display-list assertions + golden images of the new look. During Phase 2, the existing Playwright screenshot checkpoints confirm the port changed nothing visually.    |
| Playwright e2e tests (41 files)           | ~4K            | **PORT** Keep scenarios; reimplement the page-object wrappers against the React DOM (e.g. `data-testid`), drop the Qt bridge, regenerate screenshot baselines for the new look |
| Qt widget unit tests                      | —              | **REPLACE** Vitest for React components                                                                                                                                        |

## 6. Phased transition plan

Each phase is independently shippable, and the existing Qt app keeps working until Phase 8. Phases 0–3 pay off whatever is decided about Qt hosts.

### Phase 0 — Public contract coverage · already in place

- Covered by existing tests: `test_sketcher_widget.cpp` and `test_image_generation.cpp` in this repo, plus Maestro and Python consumer tests in mmshare.
- These run in every later phase; together with the Playwright suite they are the regression net for the transition.

### Phase 1 — Extract `sketcher_core/` · ~3–4 weeks

- New `sketcher_core` CMake library, no Qt link. Starts with `sketcher/rdkit/` (moved as `rdkit/`) and `model/`.
- Move MolModel, edit operations, undo and selection. Use plain classes instead of `: public QObject`, a `std::function` observer registry instead of `Q_SIGNALS`, and a small custom undo stack instead of `QUndoStack`.
- Port model tests to the new target.
- **Packaging:** add `sketcher_core` to the hand-maintained `install(FILES ...)` list in `CMakeLists.txt` (library + linker file); headers under `include/schrodinger/` are installed by directory. In `package-factory_recipes/feedstocks/sketcher/recipe.yaml`, add macOS `install_name_tool -change` lines for `libschrodinger_sketcher_core.dylib` (zlib, plus sqlite if linked directly — confirm with `otool -L`). `build.py` and requirements don't change.
- _Optional:_ split the recipe into two packages — Qt-free `sketcher-core` (`rdkit_extensions` + `sketcher_core`) and `sketcher` (Qt shim) — so Qt-free consumers don't pull in Qt.
- **Validation:** the existing Qt sketcher links `sketcher_core` + a Qt adapter and behaves identically; the sketcher C++ tests and mmshare consumer tests pass.

### Phase 2 — Extract depiction into `sketcher_core/depict/` · ~6–10 weeks

- The riskiest new piece, so it comes early.
- Move geometry out of molviewer `paint()` methods into a Qt-free display-list producer with a `FontMetrics` interface. The molviewer items become thin QPainter painters of the display list, so the existing Qt app is the first real consumer.
- Implement the default `FontMetrics` from the bundled Arimo TTF.
- Add the Qt-free SVG image API on the display list. The Qt-typed functions stay in place; retiring them is managed separately.
- **Port faithfully, with no visual changes in this phase.** Matching the old look isn't a long-term goal, but keeping it during the port means the old output can be pixel-compared to catch regressions. The refresh happens afterward, in display-list styling and the React chrome.
- **Validation:** the existing Playwright screenshot checkpoints (~147) pass unchanged, along with the image-generation tests.

### Phase 3 — Lean WASM bindings · ~2–3 weeks

- New emscripten target: `sketcher_core` + `rdkit_extensions` with embind.
- Expose ~30–50 functions: load/save (HELM, SMILES, MOL), atom/bond CRUD, undo/redo, selection, get-display-list, headless SVG.
- **Validation:** a minimal HTML page with a Canvas2D painter renders a molecule. Measure WASM size against the Qt build (expect ~10–15 MB raw vs 44 MB).

### Phase 4 — Qt-host spike (parallel with Phases 2–3) · ~2 weeks

Go/no-go gate for `QWebEngineView` in Qt hosts. Answer:

- Whether the Qt package can include QtWebEngine at all. `package-factory_recipes/feedstocks/qt/build.py` builds with `-skip qtwebengine` ("requires nss, gperf, html5lib"). Enabling it means new build deps, a Chromium build inside the Qt build and a much larger Qt package. `qtpdf` is already skipped over Chromium compile problems.
- What QtWebEngine adds to Maestro's distribution (size, process model, platform support).
- Startup ordering: QtWebEngine must be initialized before the host creates its `QApplication` (`QtWebEngineQuick::initialize()` / `AA_ShareOpenGLContexts` in C++; import QtWebEngineWidgets first in Python). Maestro and the Python apps create their own `QApplication`, so this is a host-side change.
- Startup latency and memory with N `QWebEngineView` sketchers in one PyQt panel.
- Keyboard focus and shortcuts, clipboard, drag-and-drop and theming (e.g. dark mode) across the webview boundary.
- That SIP-exposed synchronous calls behave unchanged with the core held in-process by the shim.

**If the gate fails:** Qt hosts keep a native Qt UI over `sketcher_core`. This gives up "identical UI" for Qt hosts and needs an explicit decision.

### Phase 5 — React app, web-only first · ~3–5 months

- Vite + React + TS. Canvas2D painter reading the display list from WASM.
- **Port the Playwright suite:** reimplement the page-object wrappers against the React DOM so the existing scenarios become the behavior spec for the rewrite. Regenerate screenshot baselines once the refreshed look settles.
- **UI refresh:** design the new chrome and depiction styling up front, before building out the tool and dialog breadth. Use the Monomer Sketcher PRD as input for monomer UX.
- Tier-1 scope first: ~8–10 tools (draw, erase, select, lasso, rotate, charge, atom-type, bond-type) and ~5 dialogs (atom props, bond props, save image, import text, settings). Then the remainder, in priority order: monomer palette, HELM input, everything else.
- Tool layer: each `AbstractSceneTool` subclass maps to a TS class implementing `onMouseDown/Move/Up` and producing edit commands for the core.
- **Validation:** deploy the React+WASM build alongside the Qt WASM build as "v2 web."

### Phase 6 — Qt-host integration · ~4–6 weeks

- `SketcherWidget` becomes a `QWebEngineView` shell around an in-process `sketcher_core` + QWebChannel adapter, with the same public API and signals. Hover signals map display-list hit IDs back to `RDKit::Atom*`/`Bond*`.
- Rebuild the SIP bindings against the new widget; run the Python panel test suites.
- Add `SketcherView`, the native canvas-only select-only widget, and move Maestro's ligand overlay to it (small Maestro-side change).
- Consumers will see the refreshed UI. Give Maestro and Python panel owners a preview build early.
- **Validation:** sketcher C++ tests and mmshare Maestro/Python tests pass; Maestro and Python panels work with no API changes beyond the overlay switch and the image-API migration.

### Phase 7 — Standalone Qt-free `main.cpp` · ~1–2 weeks

- Add `webview/webview` (https://github.com/webview/webview), vendored in the repo or as its own feedstock — package-factory builds shouldn't download sources via FetchContent.
- **Packaging:** Linux needs WebKitGTK (gtk3 + webkit2gtk) added to the recipe requirements; Windows needs the WebView2 SDK at build time and the WebView2 runtime on user machines; macOS needs nothing extra (WKWebView is in the OS). The app no longer links Qt.
- Rewrite `src/app/main.cpp` to instantiate `SketcherApp` + a webview window — no Qt in this binary.
- New CMake target `sketcher_standalone` linking `sketcher_core` + `rdkit_extensions` + webview.

```cpp
// main.cpp sketch
#include "webview.h"
#include "sketcher_core/sketcher_app.h"

int main() {
    webview::webview w(true, nullptr);
    w.set_title("Sketcher");
    w.set_size(1200, 800, WEBVIEW_HINT_NONE);

    SketcherApp app;  // model + depiction from sketcher_core

    w.bind("sketcher_call", [&](const std::string& req) {
        return app.handleCommand(req);     // JSON in, JSON out
    });
    app.onEvent([&](const std::string& ev) {
        w.eval("window.__sketcher_event(" + ev + ")");
    });

    w.navigate("file:///path/to/react-bundle/index.html");
    w.run();
}
```

### Phase 8 — Delete · ~2 weeks

- Remove `tool/`, `widget/`, `dialog/`, `menu/` and the molviewer QGraphicsItem code.
- Remove Qt MOC from test infrastructure for the dropped tests.

#### End state

- `sketcher_standalone` — Qt-free desktop binary (webview/webview hosts React + core)
- `sketcher.wasm` — lean web bundle (Qt-free), including headless SVG
- `sketcher` library (Qt shim) — the only Qt left: the QWebEngineView shim and the native `SketcherView` (QPainter painter) for the overlay
- Public image API is Qt-free SVG; Qt conversions live in mmshare

## 7. Tech stack picks

| Concern                 | Pick                                                               | Why                                                                                               |
| ----------------------- | ------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------- |
| UI framework            | React + TypeScript                                                 | Ecosystem + hiring; runtime size irrelevant next to WASM core                                     |
| Interactive rendering   | Canvas2D painting a C++ display list                               | Right perf envelope for 10s–100s of atoms. WebGL is overkill and harder to debug.                 |
| Headless rendering      | C++ SVG writer (Qt-free)                                           | Image API must work with no browser; same display list as the UI. Qt callers wrap SVG in mmshare. |
| Fonts                   | Arimo (already bundled) on all targets; metrics via `stb_truetype` | C++ `FontMetrics` and browser text rendering agree; removes most cross-engine drift               |
| State (UI)              | Zustand                                                            | C++ core is the source of truth; UI state is small                                                |
| Tool state machines     | Plain TS classes, XState only if a tool gets gnarly                | Most won't need it                                                                                |
| Build                   | Vite                                                               | Output a self-contained bundle the desktop webview can `file://` load                             |
| Unit tests              | Vitest                                                             | Fast, Vite-native                                                                                 |
| E2E + visual regression | Playwright (extend existing `test/wasm/`)                          | Critical for "identical look across platforms" — screenshot tests across webview engines          |
| WASM bindings           | embind                                                             | Already using emscripten flow; embind is the next step                                            |
| TS ⇄ C++ types          | Small Python codegen from a shared JSON spec                       | Don't reach for protobuf; a hand-rolled spec is plenty                                            |
| Standalone desktop host | `webview/webview`                                                  | Pure C++, system webview, ~100 KB glue. CEF if cross-engine drift becomes a real problem.         |

## 8. Risks & mitigations

> **1. QtWebEngine in Qt hosts may be unacceptable.** It's a new dependency for Maestro, with real startup cost and memory per instance, and many PyQt panels create sketchers. `QWebEngineView` does not work inside `QGraphicsProxyWidget`. Today's Qt package is built with `-skip qtwebengine`, so enabling it is a Qt recipe change with new build deps and a Chromium build. Hosts must also initialize QtWebEngine before creating their `QApplication`. _Mitigation:_ the Phase 4 spike gates Phase 6; the overlay uses native `SketcherView`; the fallback is a native Qt UI over the shared core.

> **2. Depiction correctness during extraction.** Matching today's look isn't a goal, but chemistry conventions must stay correct: wedge bonds, stereo annotations, label placement, monomer shapes. _Mitigation:_ Phase 2 is a faithful port inside the existing Qt app, checked by the existing Playwright screenshot checkpoints. After the refresh, golden images of the new look become the regression net.

> **3. Text metrics differ between C++ and browsers.** If `FontMetrics` disagrees with Canvas2D text rendering, labels overlap bonds. _Mitigation:_ one bundled font everywhere (Arimo), with C++ metrics read from the same TTF. Playwright screenshots catch any remaining mismatch.

> **4. Undo across the JS/WASM boundary is subtle.** Each edit must be a command object the core can apply/revert atomically; the _core_ owns the stack. _Mitigation:_ no React-side undo — it must be a C++ concern.

> **5. Tools crossing C++/TS every mouse move = jank.** _Mitigation:_ tools live in TS; only commit-style operations cross into the core. Display list recomputed once per commit, not per mouse move.

> **6. Consumer behavior change (Maestro + ~105 Python files).** _Mitigation:_ existing sketcher C++ tests and mmshare Maestro/Python consumer tests, run in every phase.

> **7. 36 scene tools + 30 dialogs is the actual long pole.** Phase 5 will slip without disciplined prioritization. _Mitigation:_ ship tier-1 scope first; everything else follows behind "v2 web."

> **8. Cross-engine webview drift (standalone desktop).** `webview/webview` uses WKWebView / WebView2 / WebKitGTK. _Mitigation:_ bundled font; Playwright visual regression across all three engines in CI. If drift is unacceptable, switch standalone build to CEF.

## 9. Rough total effort

**~9–14 months for one focused engineer (~5–8 months with two).** This is a judgment call, not a measurement. The largest items are Phase 5 (UI breadth, now including design work for the refresh) and Phase 2 (depiction extraction). The first draft's 3–6 months omitted the depiction layer, the headless renderer and Python/SIP compatibility.

- ~12–15K LOC of new TypeScript
- ~8–10K LOC of new/ported C++ (core API, depiction layer, painters, transport adapters)
- ~25–30K LOC of Qt code deleted

Milestones: the depiction layer running inside today's Qt app (end of Phase 2), then the lean WASM target (Phase 3, ~3–4 months in) as the first visible win.

## 10. Open questions to resolve before starting

1. **QtWebEngine in Maestro.** Is it acceptable to ship (size, licensing, platform support)? It also means un-skipping `qtwebengine` in the package-factory Qt recipe (new build deps nss/gperf/html5lib, Chromium build). Who owns that recipe, and will they take it on? Decides the Phase 4 gate.
2. **Tool scope cut.** Confirm the tier-1 tools and dialogs for v2 web.
3. **Branch strategy.** Recommended: strangler fig — incremental merges to `main` behind the core/depict adapters. Phases 0–3 fit this naturally.

#### Resolved 2026-10-08

- **Qt-typed image API:** replaced by a Qt-free SVG API; Qt wrappers live in mmshare. Retiring the Qt-typed functions is managed separately.
- **Public contract coverage:** already in place (sketcher C++ tests + Maestro/Python tests in mmshare). Phase 0 is done.
- **Render parity bar:** none. This is a UI refresh; the old look is used only as a regression check during the Phase 2 port.
- **Font strategy:** one bundled font everywhere, Arimo (already shipped by the sketcher), unless the refresh design picks a different one.

---

Sources: codebase survey of the sketcher repo (May 2026, re-checked Oct 2026), plus consumer survey of `mmshare/maestro`, `mmshare/maestrolibs` and `mmshare/python`. The legacy copy in `mmshare/src/schrodinger/sketcher/` is intentionally excluded. Decisions reflect the 2026-10-08 review.
