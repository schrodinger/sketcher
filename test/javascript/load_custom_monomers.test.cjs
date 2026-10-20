const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const test = require('node:test');
const vm = require('node:vm');

function setup({ failFirstPostMessage = false } = {}) {
  const workers = [];
  const installed = [];
  class Worker {
    constructor(url) {
      assert.equal(url, '/assets/custom_monomers.worker.js');
      this.requests = [];
      workers.push(this);
    }
    postMessage(json) {
      if (failFirstPostMessage && this === workers[0]) {
        throw new Error('Could not send input to worker');
      }
      this.requests.push(json);
    }
    terminate() {
      this.terminated = true;
    }
  }
  const Module = {
    sketcher_install_custom_monomers(database) {
      installed.push(database);
      return '';
    },
  };
  vm.runInNewContext(
    fs.readFileSync(path.join(__dirname, '../../wasm/load_custom_monomers.js'), 'utf8'),
    { Module, Worker, locateFile: (file) => `/assets/${file}`, Promise, Error },
  );
  return { Module, workers, installed };
}

// Drain promise jobs before inspecting requests submitted through the queue.
const flush = () => new Promise((resolve) => setImmediate(resolve));

test('imports run in submission order and install only completed databases', async () => {
  const { Module, workers, installed } = setup();
  const first = Module.sketcher_load_custom_monomers('first');
  const second = Module.sketcher_load_custom_monomers('second');
  assert.ok(first instanceof Promise);
  assert.deepEqual(installed, []);
  await flush();
  assert.equal(workers[0].requests[0].json, 'first');
  assert.equal(workers[0].requests[0].compiledModule, undefined);
  const compiledModule = new WebAssembly.Module(new Uint8Array([0, 97, 115, 109, 1, 0, 0, 0]));
  workers[0].onmessage({ data: {
    database: 'db1', succeeded: ['a'], failed: [], compiledModule,
  } });
  const result = await first;
  assert.deepEqual(Array.from(result.succeeded), ['a']);
  assert.equal(workers[0].terminated, true);
  await flush();
  assert.equal(workers.length, 2);
  assert.equal(workers[0].requests.length, 1);
  assert.equal(workers[1].requests[0].json, 'second');
  assert.equal(workers[1].requests[0].compiledModule, compiledModule);
  workers[1].onmessage({ data: { database: 'db2', succeeded: [], failed: ['b'] } });
  await second;
  assert.equal(workers[1].terminated, true);
  assert.deepEqual(installed, ['db1', 'db2']);
});

test('parse failures reject with Error and leave the queue and database intact', async () => {
  const { Module, workers, installed } = setup();
  const invalid = Module.sketcher_load_custom_monomers('invalid');
  const rejected = assert.rejects(invalid, /Invalid JSON/);
  const next = Module.sketcher_load_custom_monomers('valid');
  await flush();
  const compiledModule = new WebAssembly.Module(new Uint8Array([0, 97, 115, 109, 1, 0, 0, 0]));
  workers[0].onmessage({ data: { error: 'Invalid JSON', compiledModule } });
  await rejected;
  assert.equal(workers[0].terminated, true);
  assert.deepEqual(installed, []);
  await flush();
  assert.equal(workers.length, 2);
  assert.equal(workers[1].requests[0].compiledModule, compiledModule);
  workers[1].onmessage({ data: { database: 'valid', succeeded: [], failed: [] } });
  await next;
  assert.equal(workers[1].terminated, true);
  assert.deepEqual(installed, ['valid']);
});

for (const event of ['onerror', 'onmessageerror', 'initialization failure']) {
  test(`${event} terminates the worker and allows the next import to retry`, async () => {
    const { Module, workers } = setup();
    const failed = Module.sketcher_load_custom_monomers('first');
    const rejected = assert.rejects(failed, Error);
    const next = Module.sketcher_load_custom_monomers('second');
    await flush();
    if (event === 'initialization failure') {
      workers[0].onmessage({ data: { error: 'Worker failed' } });
    } else {
      workers[0][event]({ message: 'Worker failed' });
    }
    await rejected;
    assert.equal(workers[0].terminated, true);
    await flush();
    assert.equal(workers.length, 2);
    assert.equal(workers[1].requests[0].compiledModule, undefined);
    workers[1].onmessage({ data: { database: 'second', succeeded: [], failed: [] } });
    await next;
    assert.equal(workers[1].terminated, true);
  });
}

test('installation errors reject with a JavaScript Error', async () => {
  const { Module, workers } = setup();
  Module.sketcher_install_custom_monomers = () => 'Invalid database image';
  const pending = Module.sketcher_load_custom_monomers('first');
  const rejected = assert.rejects(pending, /Invalid database image/);
  await flush();
  workers[0].onmessage({ data: { database: 'invalid' } });
  await rejected;
  assert.equal(workers[0].terminated, true);
});

test('sending input failures terminate the worker and allow subsequent imports', async () => {
  const { Module, workers } = setup({ failFirstPostMessage: true });
  const first = Module.sketcher_load_custom_monomers('first');
  const rejected = assert.rejects(first, /Could not send input to worker/);
  const next = Module.sketcher_load_custom_monomers('second');
  await rejected;
  assert.equal(workers[0].terminated, true);
  await flush();
  workers[1].onmessage({ data: { database: 'second', succeeded: [], failed: [] } });
  await next;
  assert.equal(workers[1].terminated, true);
});

function setupWorker({ mimeType = 'application/wasm', failFetch = false } = {}) {
  const replies = [];
  const calls = { fetch: 0, streaming: 0, fallback: 0, instances: 0 };
  const self = { postMessage: (data, transfers) => replies.push({ data, transfers }) };
  vm.runInNewContext(
    fs.readFileSync(path.join(__dirname, '../../wasm/public/custom_monomers.worker.js'), 'utf8'),
    {
      self, Error,
      importScripts: (url) => assert.equal(url, 'MonomerLoader.js'),
      fetch: async (url, options) => {
        ++calls.fetch;
        assert.equal(url, 'MonomerLoader.wasm');
        assert.equal(options.credentials, 'same-origin');
        if (failFetch) throw new Error('Download failed');
        return new Response(new Uint8Array([0, 97, 115, 109, 1, 0, 0, 0]), {
          headers: { 'Content-Type': mimeType },
        });
      },
      WebAssembly: {
        Instance: WebAssembly.Instance,
        compileStreaming: async (response) => {
          ++calls.streaming;
          return WebAssembly.compileStreaming(response);
        },
        compile: async (bytes) => {
          ++calls.fallback;
          return WebAssembly.compile(bytes);
        },
      },
      createMonomerLoader: async (options) => {
        let receivedInstance;
        const exports = options.instantiateWasm({}, (instance, module) => {
          assert.ok(instance instanceof WebAssembly.Instance);
          assert.ok(module instanceof WebAssembly.Module);
          receivedInstance = instance;
          ++calls.instances;
        });
        assert.equal(exports, receivedInstance.exports);
        return {
          prepare_custom_monomers(json) {
            if (json === 'invalid') throw new Error('Invalid JSON');
            return { database: new Uint8Array([1, 2, 3]), succeeded: [], failed: [] };
          },
        };
      },
    },
  );
  return { self, replies, calls };
}

test('a fresh worker instantiates the shared module without fetching or compiling it', async () => {
  const first = setupWorker();
  await first.self.onmessage({ data: { json: 'first' } });
  assert.deepEqual(first.calls, { fetch: 1, streaming: 1, fallback: 0, instances: 1 });
  const { compiledModule, database } = first.replies[0].data;
  assert.ok(compiledModule instanceof WebAssembly.Module);
  assert.equal(first.replies[0].transfers[0], database.buffer);

  const second = setupWorker();
  await second.self.onmessage({ data: { json: 'second', compiledModule } });
  assert.deepEqual(second.calls, { fetch: 0, streaming: 0, fallback: 0, instances: 1 });
  assert.equal(second.replies[0].data.compiledModule, compiledModule);
});

test('workers return compiled code even when the import fails', async () => {
  const worker = setupWorker();
  await worker.self.onmessage({ data: { json: 'invalid' } });
  assert.equal(worker.replies[0].data.error, 'Invalid JSON');
  assert.ok(worker.replies[0].data.compiledModule instanceof WebAssembly.Module);
});

test('a download failure returns no compiled module to cache', async () => {
  const worker = setupWorker({ failFetch: true });
  await worker.self.onmessage({ data: { json: 'first' } });
  assert.equal(worker.replies[0].data.error, 'Download failed');
  assert.equal(worker.replies[0].data.compiledModule, undefined);
  assert.equal(worker.calls.instances, 0);
});

test('incorrect WASM MIME types still support non-streaming compilation', async () => {
  const worker = setupWorker({ mimeType: 'application/octet-stream' });
  await worker.self.onmessage({ data: { json: 'first' } });
  assert.deepEqual(worker.calls, { fetch: 1, streaming: 1, fallback: 1, instances: 1 });
  assert.ok(worker.replies[0].data.compiledModule instanceof WebAssembly.Module);
});
