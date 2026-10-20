// Included in the main Emscripten module with --post-js.
let customMonomerQueue = Promise.resolve();
// Retain executable code only; each worker owns its instance and database.
let customMonomerCompiledModule;

Module['sketcher_load_custom_monomers'] = function (json) {
  const request = customMonomerQueue.then(async () => {
    const worker = new Worker(locateFile('custom_monomers.worker.js'));
    try {
      const data = await new Promise((resolve, reject) => {
        worker.onmessage = ({ data }) => {
          if (data.compiledModule) {
            customMonomerCompiledModule = data.compiledModule;
          }
          if (data.error !== undefined) {
            reject(new Error(data.error));
          } else {
            resolve(data);
          }
        };
        worker.onerror = (event) => {
          reject(new Error(event.message || 'Custom monomer worker failed to load'));
        };
        worker.onmessageerror = () => {
          reject(new Error('Could not receive the custom monomer database'));
        };
        worker.postMessage({ json, compiledModule: customMonomerCompiledModule });
      });
      const error = Module['sketcher_install_custom_monomers'](data.database);
      if (error) {
        throw new Error(error);
      }
      return { succeeded: data.succeeded, failed: data.failed };
    } finally {
      worker.terminate();
    }
  });
  // Failed requests must not prevent subsequent imports from running.
  customMonomerQueue = request.catch(() => {});
  return request;
};
