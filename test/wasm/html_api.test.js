import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import vm from 'node:vm';

import { test } from '@playwright/test';

const html = readFileSync(new URL('../../wasm/public/wasm_shell.html', import.meta.url), 'utf8');
const script = html.match(/<script type="text\/javascript">([\s\S]*?)<\/script>/)[1];

function createPage(qtLoad) {
  const loadModule = async (...args) => {
    const instance = await qtLoad(...args);
    instance.stackSave ??= () => 0;
    instance.stackRestore ??= () => {};
    instance.incrementExceptionRefcount ??= () => {};
    instance.decrementExceptionRefcount ??= () => {};
    return instance;
  };
  const context = vm.createContext({
    window: {
      devicePixelRatio: 1,
      parent: { postMessage() {} },
      addEventListener(type, listener) {
        this[type] = listener;
      },
    },
    document: { querySelector: () => ({ style: {}, addEventListener() {} }) },
    console: { error() {} },
    schrodinger_sketcher_entry() {},
    qtLoad: loadModule,
  });
  vm.runInContext(script, context);
  return context;
}

test('HTML API waits for loading and replaces, exports, and clears molecules', async () => {
  let finishLoading;
  let molecule = 'old';
  const instance = {
    Format: { SMILES: {}, MDL_MOLV3000: {} },
    sketcher_clear() {
      molecule = '';
    },
    sketcher_import_text(text) {
      molecule += text;
    },
    sketcher_export_text(format) {
      return format === this.Format.SMILES ? molecule : `molfile:${molecule}`;
    },
  };
  const page = createPage(
    () =>
      new Promise((resolve) => {
        finishLoading = resolve;
      }),
  );
  const loading = page.init();
  const setting = page.window.setMolecule('CCO');
  assert.equal(molecule, 'old');
  finishLoading(instance);
  await loading;
  await setting;
  assert.equal(page.window.Module, instance);
  assert.equal(await page.window.getMolecule(), 'CCO');
  assert.equal(await page.window.getMolecule('MDL_MOLV3000'), 'molfile:CCO');
  await assert.rejects(page.window.getMolecule('toString'), /Unknown molecule format/);
  await assert.rejects(page.window.setMolecule(null), /Molecule must be a string/);
  assert.equal(molecule, 'CCO');
  await page.window.setMolecule('  ');
  assert.equal(await page.window.getMolecule(), '');
});

test('HTML API translates native exceptions', async () => {
  const page = createPage(async () => ({
    sketcher_clear() {},
    sketcher_import_text() {
      throw 123;
    },
    getExceptionMessage(pointer) {
      assert.equal(pointer, 123);
      return ['exception', 'Unable to determine format'];
    },
  }));
  await page.init();
  await assert.rejects(page.window.setMolecule('invalid'), /Unable to determine format/);
});

test('HTML API rejects pending calls when loading fails', async () => {
  const page = createPage(async () => {
    throw new Error('Load failed');
  });
  const result = page.window.getMolecule();
  await page.init();
  await assert.rejects(result, /Load failed/);
  await assert.rejects(page.window.setMolecule('CC'), /Load failed/);
});

function request(page, data, overrides = {}) {
  const replies = [];
  page.window.parent.postMessage = (message, origin) => {
    replies.push({ message, origin });
  };
  return page.window
    .message({
      source: page.window.parent,
      origin: 'https://embedding.example',
      data: { type: 'schrodinger.sketcher.request', id: 'request-1', ...data },
      ...overrides,
    })
    .then(() => replies);
}

test('messages wait for loading and correlate get/set responses to the parent origin', async () => {
  let finishLoading;
  let molecule = 'old';
  const page = createPage(
    () =>
      new Promise((resolve) => {
        finishLoading = resolve;
      }),
  );
  const loading = page.init();
  const pending = request(page, { method: 'setMolecule', text: 'CCO' });
  assert.equal(molecule, 'old');
  finishLoading({
    Format: { SMILES: 'smiles', MDL_MOLV3000: 'molfile' },
    sketcher_clear() {
      molecule = '';
    },
    sketcher_import_text(text) {
      molecule += text;
    },
    sketcher_export_text(format) {
      return `${format}:${molecule}`;
    },
  });
  await loading;
  const [setting] = await pending;
  assert.equal(setting.message.type, 'schrodinger.sketcher.response');
  assert.equal(setting.message.id, 'request-1');
  assert.equal(setting.message.result, null);
  assert.equal(setting.origin, 'https://embedding.example');
  const [getting] = await request(page, { id: 'request-2', method: 'getMolecule' });
  assert.equal(getting.message.id, 'request-2');
  assert.equal(getting.message.result, 'smiles:CCO');
  const [molfile] = await request(page, { method: 'getMolecule', format: 'MDL_MOLV3000' });
  assert.equal(molfile.message.result, 'molfile:CCO');
  await request(page, { method: 'setMolecule', text: '' });
  assert.equal(molecule, '');
});

test('messages return errors for invalid requests and native or loading failures', async () => {
  const page = createPage(async () => ({
    Format: { SMILES: {} },
    sketcher_clear() {},
    sketcher_import_text() {
      throw 123;
    },
    getExceptionMessage() {
      return ['exception', 'Invalid molecule'];
    },
  }));
  await page.init();
  for (const [data, message] of [
    [{ method: 'setMolecule', text: 'bad' }, 'Invalid molecule'],
    [{ method: 'setMolecule', text: null }, 'Molecule must be a string'],
    [{ method: 'getMolecule', format: 'bogus' }, 'Unknown molecule format: bogus'],
    [{ method: 'sketcher_clear' }, 'Unknown sketcher method: sketcher_clear'],
  ]) {
    const [reply] = await request(page, data);
    assert.equal(reply.message.error.message, message);
    assert.equal(reply.message.id, 'request-1');
    assert.equal('result' in reply.message, false);
  }
  const failed = createPage(async () => {
    throw new Error('Load failed');
  });
  const pending = request(failed, { method: 'getMolecule' });
  await failed.init();
  // The mock Error comes from outside the VM and fails its instanceof Error check.
  assert.equal((await pending)[0].message.error.message, 'Error: Load failed');
});

test('messages ignore unrelated senders, opaque origins, and malformed envelopes', async () => {
  const page = createPage(() => {
    throw new Error('Must not load');
  });
  for (const overrides of [
    { source: {} },
    { source: page.window },
    { source: null },
    { origin: 'null' },
    { data: null },
    { data: { type: 'other', id: '1' } },
    { data: { type: 'schrodinger.sketcher.request' } },
  ]) {
    assert.deepEqual(await request(page, { method: 'getMolecule' }, overrides), []);
  }
});

test('image messages wait for loading, select formats, and report export errors', async () => {
  let finishLoading;
  const page = createPage(
    () =>
      new Promise((resolve) => {
        finishLoading = resolve;
      }),
  );
  const loading = page.init();
  const pending = request(page, { method: 'getImage' });
  let fail = false;
  finishLoading({
    ImageFormat: { PNG: 1, SVG: 2 },
    sketcher_export_image(format) {
      if (fail) throw 123;
      return format === 1 ? 'png-base64' : 'svg-base64';
    },
    getExceptionMessage() {
      return ['exception', 'Image export failed'];
    },
  });
  await loading;
  assert.equal((await pending)[0].message.result, 'png-base64');
  const [svg] = await request(page, { method: 'getImage', format: 'SVG' });
  assert.equal(svg.message.result, 'svg-base64');
  const [invalid] = await request(page, { method: 'getImage', format: 'toString' });
  assert.equal(invalid.message.error.message, 'Unknown image format: toString');
  fail = true;
  const [failed] = await request(page, { method: 'getImage' });
  assert.equal(failed.message.error.message, 'Image export failed');
});

test('native errors release exception storage and restore the WASM stack', async () => {
  const calls = [];
  const page = createPage(async () => ({
    stackSave() {
      calls.push('save');
      return 42;
    },
    stackRestore(stack) {
      assert.equal(stack, 42);
      calls.push('restore');
    },
    sketcher_clear() {},
    sketcher_import_text() {
      throw 123;
    },
    incrementExceptionRefcount(pointer) {
      assert.equal(pointer, 123);
      calls.push('retain');
    },
    decrementExceptionRefcount(pointer) {
      assert.equal(pointer, 123);
      calls.push('release');
    },
    getExceptionMessage() {
      calls.push('message');
      return ['exception', 'Invalid molecule'];
    },
  }));
  await page.init();
  await assert.rejects(page.window.setMolecule('bad'), /Invalid molecule/);
  assert.deepEqual(calls, ['save', 'retain', 'message', 'release', 'restore']);
  calls.length = 0;
  page.window.Module.getExceptionMessage = () => {
    throw new Error('Message lookup failed');
  };
  await assert.rejects(page.window.setMolecule('bad'), /Message lookup failed/);
  assert.deepEqual(calls, ['save', 'retain', 'release', 'restore']);
  calls.length = 0;
  await page.window.setMolecule('');
  assert.deepEqual(calls, ['save', 'restore']);
});
