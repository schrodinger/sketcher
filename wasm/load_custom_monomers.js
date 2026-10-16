// Included in the main Emscripten module with --post-js.
let customMonomerWorker;
let customMonomerQueue = Promise.resolve();

Module['sketcher_load_custom_monomers'] = function (json) {
  const request = customMonomerQueue.then(
    () =>
      new Promise((resolve, reject) => {
        try {
          if (!customMonomerWorker) {
            customMonomerWorker = new Worker(locateFile('custom_monomers.worker.js'));
          }
          const failWorker = (message) => {
            customMonomerWorker.terminate();
            customMonomerWorker = undefined;
            reject(new Error(message));
          };
          customMonomerWorker.onmessage = ({ data }) => {
            if (data.error !== undefined) {
              if (data.fatal) {
                failWorker(data.error);
              } else {
                reject(new Error(data.error));
              }
              return;
            }
            try {
              const error = Module['sketcher_install_custom_monomers'](data.database);
              if (error) {
                throw new Error(error);
              }
              resolve({ succeeded: data.succeeded, failed: data.failed });
            } catch (error) {
              reject(error);
            }
          };
          customMonomerWorker.onerror = (event) => {
            failWorker(event.message || 'Custom monomer worker failed to load');
          };
          customMonomerWorker.onmessageerror = () => {
            failWorker('Could not receive the custom monomer database');
          };
          customMonomerWorker.postMessage(json);
        } catch (error) {
          reject(error);
        }
      }),
  );
  // Failed requests must not prevent subsequent imports from running.
  customMonomerQueue = request.catch(() => {});
  return request;
};
