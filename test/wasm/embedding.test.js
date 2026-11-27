import { expect, test } from '@playwright/test';
import { createServer } from 'node:http';

let parentServer;
let parentOrigin;

test.beforeAll(async () => {
  parentServer = createServer((request, response) => {
    response.setHeader('Content-Type', 'text/html');
    response.end('<!doctype html><title>Embedding host</title>');
  });
  await new Promise((resolve) => parentServer.listen(0, '127.0.0.1', resolve));
  parentOrigin = `http://127.0.0.1:${parentServer.address().port}`;
});

test.afterAll(async () => {
  await new Promise((resolve, reject) =>
    parentServer.close((error) => (error ? reject(error) : resolve())),
  );
});

async function openEmbedding(page, baseURL) {
  await page.goto(parentOrigin);
  await page.evaluate(async (sketcherURL) => {
    const iframe = document.createElement('iframe');
    const origin = new URL(sketcherURL).origin;
    let nextId = 0;
    const pending = new Map();
    window.addEventListener('message', (event) => {
      if (event.source !== iframe.contentWindow || event.origin !== origin) return;
      const response = event.data;
      if (response?.type !== 'schrodinger.sketcher.response') return;
      const call = pending.get(response.id);
      if (!call) return;
      clearTimeout(call.timer);
      pending.delete(response.id);
      if (response.error) call.reject(new Error(response.error.message));
      else call.resolve(response.result);
    });
    function request(method, params = {}) {
      return new Promise((resolve, reject) => {
        const id = String(++nextId);
        const timer = setTimeout(() => {
          pending.delete(id);
          reject(new Error('Sketcher request timed out'));
        }, 20000);
        pending.set(id, { resolve, reject, timer });
        iframe.contentWindow.postMessage(
          {
            type: 'schrodinger.sketcher.request',
            id,
            method,
            ...params,
          },
          origin,
        );
      });
    }
    const loaded = new Promise((resolve) =>
      iframe.addEventListener('load', resolve, { once: true }),
    );
    iframe.src = sketcherURL;
    document.body.appendChild(iframe);
    await loaded;
    window.sketcherRequest = request;
  }, `${baseURL}/wasm_shell.html`);
}

test('cross-origin parent gets molecules and images through messages', async ({
  page,
  baseURL,
}) => {
  await openEmbedding(page, baseURL);
  const result = await page.evaluate(async () => {
    const request = window.sketcherRequest;
    const iframe = document.querySelector('iframe');
    let directAccessBlocked = false;
    try {
      void iframe.contentWindow.document;
    } catch {
      directAccessBlocked = true;
    }
    await request('setMolecule', { text: 'CCO' });
    const [smiles, molfile] = await Promise.all([
      request('getMolecule'),
      request('getMolecule', { format: 'MDL_MOLV3000' }),
    ]);
    await request('setMolecule', { text: molfile });
    const roundTrip = await request('getMolecule');
    const png = await request('getImage');
    const svg = await request('getImage', { format: 'SVG' });
    const imageSizes = [];
    for (const [mime, data] of [
      ['image/png', png],
      ['image/svg+xml', svg],
    ]) {
      const image = new Image();
      image.src = `data:${mime};base64,${data}`;
      await image.decode();
      imageSizes.push([image.naturalWidth, image.naturalHeight]);
    }
    let imageError;
    try {
      await request('getImage', { format: 'JPEG' });
    } catch (err) {
      imageError = err.message;
    }
    let moleculeFormatError;
    try {
      await request('getMolecule', { format: 'NOT_A_FORMAT' });
    } catch (err) {
      moleculeFormatError = err.message;
    }
    let moleculeTypeError;
    try {
      await request('setMolecule', { text: null });
    } catch (err) {
      moleculeTypeError = err.message;
    }
    let methodError;
    try {
      await request('unknownMethod');
    } catch (err) {
      methodError = err.message;
    }
    let error;
    try {
      await request('setMolecule', { text: 'foobar' });
    } catch (err) {
      error = err.message;
    }
    await request('setMolecule', { text: '' });
    const empty = await request('getMolecule');
    return {
      directAccessBlocked,
      smiles,
      molfile,
      roundTrip,
      error,
      empty,
      png,
      svg,
      imageSizes,
      imageError,
      moleculeFormatError,
      moleculeTypeError,
      methodError,
    };
  });
  expect(result.directAccessBlocked).toBe(true);
  expect(result.smiles).toBe('CCO');
  expect(result.molfile).toContain('V3000');
  expect(result.roundTrip).toBe('CCO');
  expect(result.error).toBe('Unable to determine format');
  expect(result.empty).toBe('');
  expect([...Buffer.from(result.png, 'base64').subarray(0, 8)]).toEqual([
    0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a,
  ]);
  expect(Buffer.from(result.svg, 'base64').toString('utf8')).toContain('<svg');
  for (const [width, height] of result.imageSizes) {
    expect(width).toBeGreaterThan(0);
    expect(height).toBeGreaterThan(0);
  }
  expect(result.imageError).toBe('Unknown image format: JPEG');
  expect(result.moleculeFormatError).toBe('Unknown molecule format: NOT_A_FORMAT');
  expect(result.moleculeTypeError).toBe('Molecule must be a string');
  expect(result.methodError).toBe('Unknown sketcher method: unknownMethod');
});

test('ignores requests from a sibling frame', async ({ page, baseURL }) => {
  await openEmbedding(page, baseURL);
  const sketcherFrame = page.frames().find((frame) => frame.url().endsWith('/wasm_shell.html'));
  await page.evaluate(async () => {
    await window.sketcherRequest('setMolecule', { text: 'CCO' });
    const sibling = document.createElement('iframe');
    const sketcherOrigin = new URL(document.querySelector('iframe').src).origin;
    sibling.srcdoc = `<!doctype html><script>
      parent.postMessage({ type: 'sibling-ready' }, '*');
      setTimeout(() => parent.frames[0].postMessage({
          type: 'schrodinger.sketcher.request',
          id: 'sibling-request',
          method: 'setMolecule',
          text: 'CCN',
        }, '${sketcherOrigin}'), 100);
    <\/script>`;
    document.body.appendChild(sibling);
    await new Promise((resolve) => {
      function onMessage(event) {
        if (event.data?.type === 'sibling-ready') {
          removeEventListener('message', onMessage);
          resolve();
        }
      }
      addEventListener('message', onMessage);
    });
  });

  const result = await sketcherFrame.evaluate(async () => {
    const original = await window.getMolecule();
    const response = new Promise((resolve) => {
      const timer = setTimeout(() => resolve(null), 500);
      window.addEventListener(
        'message',
        (event) => {
          if (
            event.data?.type === 'schrodinger.sketcher.response' &&
            event.data.id === 'sibling-request'
          ) {
            clearTimeout(timer);
            resolve(event.data);
          }
        },
        { once: true },
      );
    });
    const reply = await response;
    return { reply, molecule: await window.getMolecule(), original };
  });

  expect(result.reply).toBeNull();
  expect(result.original).toBe('CCO');
  expect(result.molecule).toBe('CCO');
});

test('ignores requests from an opaque-origin parent', async ({ page, baseURL }) => {
  await page.goto(baseURL);
  await page.evaluate(async (sketcherURL) => {
    const opaqueParent = document.createElement('iframe');
    opaqueParent.sandbox = 'allow-scripts';
    opaqueParent.srcdoc = `<!doctype html><iframe src="${sketcherURL}"></iframe><script>
      const sketcher = document.querySelector('iframe');
      sketcher.addEventListener('load', () => {
        sketcher.contentWindow.postMessage({
          type: 'schrodinger.sketcher.request',
          id: 'opaque-request',
          method: 'setMolecule',
          text: 'CCN',
        }, '*');
        parent.postMessage({ type: 'opaque-request-sent' }, '*');
      });
      addEventListener('message', (event) => {
        if (event.source === sketcher.contentWindow && event.data?.id === 'opaque-request') {
          parent.postMessage({ type: 'opaque-response', request: event.data }, '*');
        }
      });
    <\/script>`;
    document.body.appendChild(opaqueParent);
    await new Promise((resolve) => {
      function onMessage(event) {
        if (event.data?.type === 'opaque-request-sent') {
          removeEventListener('message', onMessage);
          resolve();
        }
      }
      addEventListener('message', onMessage);
    });
  }, `${baseURL}/wasm_shell.html`);

  const result = await page.evaluate(async () => {
    const response = new Promise((resolve) => {
      const timer = setTimeout(() => resolve(null), 250);
      window.addEventListener(
        'message',
        (event) => {
          if (event.data?.type === 'opaque-response') {
            clearTimeout(timer);
            resolve(event.data);
          }
        },
        { once: true },
      );
    });
    const reply = await response;
    return { reply };
  });

  expect(result.reply).toBeNull();
});

const malformedSmiles = [
  { name: 'unclosed ring', text: 'C1CC' },
  { name: 'unclosed branch', text: 'CC(O' },
  { name: 'unclosed atom bracket', text: 'C[NH3+' },
  { name: 'dangling bond', text: 'CC=' },
];

for (const { name, text } of malformedSmiles) {
  test(`cross-origin setMolecule rejects SMILES with ${name} and recovers`, async ({
    page,
    baseURL,
  }) => {
    await openEmbedding(page, baseURL);
    const result = await page.evaluate(async (text) => {
      const request = window.sketcherRequest;
      await request('setMolecule', { text: 'CCO' });
      let error;
      try {
        await request('setMolecule', { text });
      } catch (err) {
        error = { name: err.name, message: err.message };
      }
      const afterFailure = await request('getMolecule');
      await request('setMolecule', { text: 'CCN' });
      const recovered = await request('getMolecule');
      const image = await request('getImage');
      return { error, afterFailure, recovered, image };
    }, text);
    expect(result.error).toEqual({ name: 'Error', message: 'Unable to determine format' });
    expect(result.afterFailure).toBe('');
    expect(result.recovered).toBe('CCN');
    expect([...Buffer.from(result.image, 'base64').subarray(0, 8)]).toEqual([
      0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a,
    ]);
  });
}
