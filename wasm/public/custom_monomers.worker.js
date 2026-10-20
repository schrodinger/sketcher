importScripts('MonomerLoader.js');

async function compileMonomerLoader() {
  const response = await fetch('MonomerLoader.wasm', { credentials: 'same-origin' });
  if (!response.ok) {
    throw new Error(`Could not load monomer loader WASM: HTTP ${response.status}`);
  }
  if (WebAssembly.compileStreaming) {
    try {
      return await WebAssembly.compileStreaming(response.clone());
    } catch {
      // Preserve the generated loader's fallback for servers with a wrong MIME type.
    }
  }
  return WebAssembly.compile(await response.arrayBuffer());
}

self.onmessage = async ({ data: { json, compiledModule } }) => {
  try {
    compiledModule ??= await compileMonomerLoader();
    const module = await createMonomerLoader({
      instantiateWasm(imports, receiveInstance) {
        const instance = new WebAssembly.Instance(compiledModule, imports);
        receiveInstance(instance, compiledModule);
        return instance.exports;
      },
    });
    const result = module.prepare_custom_monomers(json);
    self.postMessage(
      { ...result, compiledModule },
      result.database ? [result.database.buffer] : [],
    );
  } catch (error) {
    self.postMessage({
      error: error instanceof Error ? error.message : String(error),
      compiledModule,
    });
  }
};
