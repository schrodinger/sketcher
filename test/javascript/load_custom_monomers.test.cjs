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
  assert.deepEqual(workers[0].requests, ['first']);
  workers[0].onmessage({ data: { database: 'db1', succeeded: ['a'], failed: [] } });
  const result = await first;
  assert.deepEqual(Array.from(result.succeeded), ['a']);
  assert.equal(workers[0].terminated, true);
  await flush();
  assert.equal(workers.length, 2);
  assert.deepEqual(workers[0].requests, ['first']);
  assert.deepEqual(workers[1].requests, ['second']);
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
  workers[0].onmessage({ data: { error: 'Invalid JSON' } });
  await rejected;
  assert.equal(workers[0].terminated, true);
  assert.deepEqual(installed, []);
  await flush();
  assert.equal(workers.length, 2);
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
