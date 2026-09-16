import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { dirname, resolve } from 'node:path';
import test from 'node:test';
import { fileURLToPath } from 'node:url';
import vm from 'node:vm';
import { gzipSync } from 'node:zlib';

const HERE = dirname(fileURLToPath(import.meta.url));
const LOADER = resolve(HERE, '../../../wasm/compressed_wasm_loader.js');
const SOURCE = new TextEncoder().encode('compiled Sketcher WebAssembly');
const LOADER_SOURCE = await readFile(LOADER, 'utf8');

function installLoader(moduleConfig, fetchImpl, overrides = {}) {
  const wasm = {
    compileStreaming: async (response) => new Uint8Array(await response.arrayBuffer()),
    compile: async (bytes) => new Uint8Array(bytes),
    instantiate: async (compiledModule) => ({ exports: { compiledModule } }),
    ...overrides.WebAssembly,
  };
  const context = vm.createContext({
    console,
    DecompressionStream,
    Error,
    fetch: fetchImpl,
    Module: moduleConfig,
    Response,
    WebAssembly: wasm,
    ...overrides.globals,
  });
  vm.runInContext(LOADER_SOURCE, context, { filename: LOADER });
  return context.Module;
}

async function instantiate(moduleConfig) {
  return new Promise((resolveResult, reject) => {
    const timeout = setTimeout(() => reject(new Error('instantiateWasm timed out')), 1000);
    const returnValue = moduleConfig.instantiateWasm({}, (instance, compiledModule) => {
      clearTimeout(timeout);
      resolveResult({ compiledModule, instance, returnValue });
    });
  });
}

test('loads and expands a static gzip sibling while preserving the cache query', async () => {
  const urls = [];
  // The production post-processing adds the query before this loader runs.
  const source = LOADER_SOURCE.replace(
    "const wasmPath = 'Sketcher.wasm';",
    "const wasmPath = 'Sketcher.wasm?cache_bust=abc123';",
  );
  const context = vm.createContext({
    console,
    DecompressionStream,
    Error,
    fetch: async (url) => {
      urls.push(url);
      return new Response(gzipSync(SOURCE));
    },
    Module: {
      locateFile: (filename, prefix) => `${prefix}sketcher-assets/${filename}`,
    },
    Response,
    scriptDirectory: '/app/',
    WebAssembly: {
      compileStreaming: async (response) => new Uint8Array(await response.arrayBuffer()),
      compile: async (bytes) => new Uint8Array(bytes),
      instantiate: async (compiledModule) => ({ exports: { compiledModule } }),
    },
  });
  vm.runInContext(source, context, { filename: LOADER });

  const { compiledModule, returnValue } = await instantiate(context.Module);
  assert.deepEqual(Array.from(compiledModule), Array.from(SOURCE));
  assert.deepEqual(urls, ['/app/sketcher-assets/Sketcher.wasm.gz?cache_bust=abc123']);
  assert.equal(Object.keys(returnValue).length, 0);
});

test('accepts a gzip response already decoded by Content-Encoding', async () => {
  class UnexpectedDecompressionStream {
    constructor() {
      throw new Error('manual decompression should not run');
    }
  }

  const moduleConfig = installLoader(
    {},
    async () => new Response(SOURCE, { headers: { 'Content-Encoding': 'gzip' } }),
    { globals: { DecompressionStream: UnexpectedDecompressionStream } },
  );
  const { compiledModule } = await instantiate(moduleConfig);
  assert.deepEqual(Array.from(compiledModule), Array.from(SOURCE));
});

test('falls back to raw WASM when the gzip sibling is absent', async () => {
  const urls = [];
  const moduleConfig = installLoader({}, async (url) => {
    urls.push(url);
    return url.includes('.gz')
      ? new Response(null, { status: 404 })
      : new Response(SOURCE, { headers: { 'Content-Type': 'application/octet-stream' } });
  });

  const { compiledModule } = await instantiate(moduleConfig);
  assert.deepEqual(Array.from(compiledModule), Array.from(SOURCE));
  assert.deepEqual(urls, ['Sketcher.wasm.gz', 'Sketcher.wasm']);
});

test('preserves a consumer-provided instantiateWasm hook', () => {
  const existing = () => ({ custom: true });
  const moduleConfig = installLoader({ instantiateWasm: existing }, async () => {
    throw new Error('fetch should not run');
  });
  assert.equal(moduleConfig.instantiateWasm, existing);
});

test('aborts module startup when neither compressed nor raw WASM can load', async () => {
  let reportAbort;
  const aborted = new Promise((resolveAbort) => {
    reportAbort = resolveAbort;
  });
  const moduleConfig = installLoader(
    {},
    async () => new Response(null, { status: 404, statusText: 'Not Found' }),
    {
      globals: {
        abort: (message) => {
          reportAbort(message);
          throw new Error(message);
        },
      },
    },
  );

  moduleConfig.instantiateWasm({}, () => {
    throw new Error('success callback should not run');
  });
  assert.match(
    await aborted,
    /Failed to load Sketcher WebAssembly: fetch Sketcher\.wasm failed: 404 Not Found/,
  );
});
