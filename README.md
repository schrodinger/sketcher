[![sketcher](https://github.com/schrodinger/sketcher/blob/main/.github/schrodinger-sketcher-logo.png)](https://www.schrodinger.com/2dsketcher)<br />

![build-status](https://github.com/schrodinger/sketcher/actions/workflows/nightly-build.yml/badge.svg?branch=main)
[![latest release](https://img.shields.io/github/v/release/schrodinger/sketcher)](https://github.com/schrodinger/sketcher/releases)

The Schrödinger Sketcher is an open-source application for drawing and editing chemical structures. Sketcher is entirely built on the RDKit open source toolkit, which serves as its underlying chemical model. It provides a user-friendly interface for creating molecules, reactions, and other chemical diagrams, which can be used independently or integrated into other cheminformatics workflows.

This project is released by Schrödinger, Inc. and is available under an open-source license to foster collaboration and development within the scientific community.

## Features

- Intuitive drawing tools for atoms, bonds, rings, and functional groups.
- Support for many chemical file formats.
- Cleanup and layout algorithms for generating clear 2D representations.
- Support for stereochemistry representation.
- Integration capabilities for use in other applications, including as a web component.

## Access and Training

**[Open Access Online Version](https://www.schrodinger.com/2dsketcher)** -- Explore Sketcher directly in your web browser.

**[Training and Video Walkthrough](https://www.schrodinger.com/sites/default/files/s3/public/2D-Sketcher/2023-2/Content/Resources/Videos/2D_Sketcher.mp4)** -- Learn how to use Sketcher with this guided video.

## Embedding the HTML page

Embed the WebAssembly distribution's `wasm_shell.html` in an iframe. The parent
can get and set molecules across origins using `postMessage`:

```html
<iframe
  id="sketcher"
  src="https://sketcher.example/wasm_shell.html"
  width="800"
  height="600"
></iframe>
<script>
  const iframe = document.getElementById('sketcher');
  const sketcherOrigin = new URL(iframe.src).origin;
  const pending = new Map();
  let nextId = 0;

  window.addEventListener('message', (event) => {
    if (event.source !== iframe.contentWindow || event.origin !== sketcherOrigin) return;
    const response = event.data;
    if (!response || response.type !== 'schrodinger.sketcher.response') return;
    const call = pending.get(response.id);
    if (!call) return;
    pending.delete(response.id);
    clearTimeout(call.timer);
    if (response.error) call.reject(new Error(response.error.message));
    else call.resolve(response.result);
  });

  function sketcherRequest(method, params = {}) {
    return new Promise((resolve, reject) => {
      const id = String(++nextId);
      const timer = setTimeout(() => {
        pending.delete(id);
        reject(new Error('Sketcher request timed out'));
      }, 60000);
      pending.set(id, { resolve, reject, timer });
      iframe.contentWindow.postMessage(
        {
          ...params,
          type: 'schrodinger.sketcher.request',
          id,
          method,
        },
        sketcherOrigin,
      );
    });
  }

  iframe.addEventListener('load', async () => {
    try {
      await sketcherRequest('setMolecule', { text: 'CCO' });
      const smiles = await sketcherRequest('getMolecule');
      const molfile = await sketcherRequest('getMolecule', { format: 'MDL_MOLV3000' });
      const png = await sketcherRequest('getImage');
      const image = document.createElement('img');
      image.src = `data:image/png;base64,${png}`;
      document.body.appendChild(image);
      console.log(smiles, molfile);
    } catch (error) {
      console.error(error);
    }
  });
</script>
```

Replace the example URL with your sketcher deployment. Send requests after the
iframe's `load` event; requests then wait for WASM initialization internally.
Each request has type `schrodinger.sketcher.request`, a unique string `id`, and
one of these methods:

| Method        | Parameters                                                  | Result                     |
| ------------- | ----------------------------------------------------------- | -------------------------- |
| `getImage`    | Optional `format`: `PNG` (default) or `SVG`                 | Base64-encoded image bytes |
| `getMolecule` | Optional `format`, default `SMILES`                         | Molecule text              |
| `setMolecule` | Required string `text`, input format detected automatically | `null`                     |

Replies have type `schrodinger.sketcher.response`, the request's `id`, and either
`result` or `error: { message: "..." }`. Loading failures, invalid molecule input,
and unknown or unsupported methods/formats return errors. Format names come
from `Module.Format` for molecules, for example `SMILES`, `MDL_MOLV3000`, or `MAESTRO`.

`setMolecule` replaces the canvas contents; an empty or whitespace-only string
clears it. The canvas is cleared before importing, so invalid molecule text
leaves it empty.

`getImage` exports the current canvas through the WASM `sketcher_export_image`
binding to `SketcherWidget::getImageBytes`. Its result is base64 without a data
URL prefix. Use `data:image/png;base64,` for PNG or
`data:image/svg+xml;base64,` for SVG when displaying it in an image element.

Only the immediate embedding parent may send requests. Replies target the
sender's exact origin. Unrelated messages and requests without a string ID are
ignored. Serve both applications over HTTP(S); opaque parent origins such as
sandboxed pages without `allow-same-origin` are not supported. Hosting headers
must permit framing by your application.

Same-origin hosts can also call the promise-based `window.getMolecule(format)`, `window.getImage(format)`, and `window.setMolecule(text)` methods directly. The existing `window.Module`
API remains available for advanced use.

Older WASM distributions without the exception cleanup exports remain compatible,
but cannot release native exception storage through this interface. Use a build
with the updated exports for long-lived embeds that may receive invalid input.

## Build Prerequisites

- A C++ compiler supporting C++20 or later
- CMake (version 3.24 or later)
- All [required dependencies](https://github.com/schrodinger/sketcher/blob/main/external/versions.json) installed and accessible to CMake

## Support and Community

For questions, support, or discussions related to the Schrödinger Sketcher, please use the [GitHub issue tracker](https://github.com/schrodinger/sketcher/issues).

To contribute code, please follow these steps:

- **Fork the repository:** Fork the `schrodinger/sketcher` repository to your own GitHub account.
- **Create a branch:** Create a new branch in your forked repository for your changes.
- **Implement your changes:** Make your code changes or additions.
- **Test your changes:** Ensure your changes build correctly and pass any existing tests. Add new tests for new features if applicable.
- **Submit a Pull Request:** Open a pull request from your branch to the `main` branch of the `schrodinger/sketcher` repository. Please provide a clear description of your changes.

All contributions are subject to the terms of the BSD license.

## License

The code is released under the [BSD 3-Clause License](https://github.com/schrodinger/sketcher/blob/master/LICENSE).
