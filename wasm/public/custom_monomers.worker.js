importScripts('MonomerLoader.js');

const loader = createMonomerLoader();

self.onmessage = async ({ data: json }) => {
  try {
    const module = await loader;
    const result = module.prepare_custom_monomers(json);
    self.postMessage(result, result.database ? [result.database.buffer] : []);
  } catch (error) {
    self.postMessage({
      error: error instanceof Error ? error.message : String(error),
    });
  }
};
