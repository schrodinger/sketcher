/*
 * Emscripten --pre-js hook for the main Sketcher module.
 *
 * Published npm packages contain Sketcher.wasm.gz instead of the much larger
 * raw WebAssembly file. Development builds keep the raw file, so try the gzip
 * sibling first and fall back without requiring special static-host headers.
 */
(() => {
  if (Module['instantiateWasm']) {
    return;
  }

  // update_cache_bust.py replaces this filename with a hashed URL in the
  // generated Sketcher.js. gzipSiblingUrl keeps that query on the .gz request.
  const wasmPath = 'Sketcher.wasm';
  const fetchOptions = { credentials: 'same-origin' };

  function gzipSiblingUrl(url) {
    const fragmentIndex = url.indexOf('#');
    const fragment = fragmentIndex === -1 ? '' : url.slice(fragmentIndex);
    const withoutFragment = fragmentIndex === -1 ? url : url.slice(0, fragmentIndex);
    const queryIndex = withoutFragment.indexOf('?');
    const query = queryIndex === -1 ? '' : withoutFragment.slice(queryIndex);
    const pathname = queryIndex === -1 ? withoutFragment : withoutFragment.slice(0, queryIndex);
    return `${pathname}.gz${query}${fragment}`;
  }

  function hasGzipContentEncoding(response) {
    return response.headers
      .get('content-encoding')
      ?.split(',')
      .some((encoding) => encoding.trim().toLowerCase() === 'gzip');
  }

  async function compileResponse(response, decompress) {
    if (!response.body) {
      throw new Error('WebAssembly response has no body');
    }

    let body = response.body;
    if (decompress && !hasGzipContentEncoding(response)) {
      if (typeof globalThis.DecompressionStream !== 'function') {
        throw new Error('DecompressionStream is unavailable');
      }
      body = body.pipeThrough(new globalThis.DecompressionStream('gzip'));
    }

    // Static hosts often serve .gz files as application/gzip (or raw .wasm
    // files as octet-stream). Give compileStreaming the MIME type it requires.
    const wasmResponse = new Response(body, {
      headers: { 'Content-Type': 'application/wasm' },
    });

    if (typeof WebAssembly.compileStreaming === 'function') {
      return WebAssembly.compileStreaming(wasmResponse);
    }
    return WebAssembly.compile(await wasmResponse.arrayBuffer());
  }

  async function compileSketcherWasm() {
    const scriptPrefix = typeof scriptDirectory === 'string' ? scriptDirectory : '';
    const wasmUrl = Module['locateFile']
      ? Module['locateFile'](wasmPath, scriptPrefix)
      : `${scriptPrefix}${wasmPath}`;

    try {
      const gzipResponse = await fetch(gzipSiblingUrl(wasmUrl), fetchOptions);
      if (gzipResponse.ok) {
        return await compileResponse(gzipResponse, true);
      }
    } catch (_) {
      // Existing deployments only have the raw asset. Preserve that path when
      // the gzip sibling is missing, cannot be decoded, or cannot be compiled.
    }

    const rawResponse = await fetch(wasmUrl, fetchOptions);
    if (!rawResponse.ok) {
      throw new Error(
        `fetch ${wasmUrl} failed: ${rawResponse.status} ${rawResponse.statusText}`,
      );
    }
    return compileResponse(rawResponse, false);
  }

  Module['instantiateWasm'] = (imports, successCallback) => {
    compileSketcherWasm()
      .then(async (compiledModule) => {
        const instance = await WebAssembly.instantiate(compiledModule, imports);
        successCallback(instance, compiledModule);
      })
      .catch((error) => {
        const message = error instanceof Error ? error.message : String(error);
        const failure = `Failed to load Sketcher WebAssembly: ${message}`;
        if (typeof abort === 'function') {
          // Emscripten's abort rejects the module-ready promise as well as
          // notifying Qt. Swallow its deliberate throw inside this callback.
          try {
            abort(failure);
          } catch (_) {
            // abort already reported the failure and rejected module startup.
          }
        } else {
          Module['onAbort']?.(failure);
        }
      });
    return {};
  };
})();
