# Embedding the WebAssembly sketcher

Embed the WebAssembly distribution's `wasm_shell.html` in an iframe. This lets
another web application display the sketcher and exchange molecules or rendered
images with it, even when the two applications are hosted on different origins.

```html
<iframe id="sketcher" src="https://sketcher.example/wasm_shell.html"></iframe>
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
        { ...params, type: 'schrodinger.sketcher.request', id, method },
        sketcherOrigin,
      );
    });
  }

  iframe.addEventListener('load', async () => {
    await sketcherRequest('setMolecule', { text: 'CCO' });
    const smiles = await sketcherRequest('getMolecule');
    const png = await sketcherRequest('getImage');
    console.log(smiles, png);
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

Replies have type `schrodinger.sketcher.response`, the request's `id`, and
either `result` or `error: { message: "..." }`. Loading failures, invalid
molecule input, and unknown or unsupported methods or formats return errors.
Molecule format names come from `Module.Format`, for example `SMILES`,
`MDL_MOLV3000`, or `MAESTRO`.

`setMolecule` replaces the canvas contents; an empty or whitespace-only string
clears it. The canvas is cleared before importing, so invalid molecule text
leaves it empty.

`getImage` returns base64 without a data URL prefix. Use
`data:image/png;base64,` for PNG or `data:image/svg+xml;base64,` for SVG when
displaying it in an image element.

Only the immediate embedding parent may send requests. This prevents unrelated
windows and nested frames from controlling the sketcher. Replies target the
sender's exact origin. Unrelated messages and requests without a string ID are
ignored. Serve both applications over HTTP(S); opaque parent origins such as
sandboxed pages without `allow-same-origin` are not supported. Hosting headers
must permit framing by your application.

Same-origin hosts can call `window.getMolecule(format)`,
`window.getImage(format)`, and `window.setMolecule(text)` directly. These methods
return promises. The existing `window.Module` API remains available for advanced
use.
