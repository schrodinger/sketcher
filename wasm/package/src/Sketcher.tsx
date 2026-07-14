import type { ComponentProps } from 'react';
import { forwardRef, useCallback, useEffect, useImperativeHandle, useRef, useState } from 'react';

import type { Merge, Tagged } from 'type-fest';

declare global {
  interface Window {
    Module?: SketcherWASM;
  }
}

export const SketcherImageFormat = {
  PNG: 'PNG',
  SVG: 'SVG',
} as const;
export type SketcherImageFormat = (typeof SketcherImageFormat)[keyof typeof SketcherImageFormat];

/**
 * @enum {string}
 *
 * NOTE: This mirrors the C++ `Format` enum which is the source of truth -
 * @see {@link ../../../include/schrodinger/rdkit_extensions/file_format.h}
 * It is exposed to the WASM module via the `EMSCRIPTEN_BINDINGS` block in
 * @see {@link ../../../src/app/main.cpp}
 */
export const RepresentationFormat = {
  AUTO_DETECT: 'AUTO_DETECT',
  RDMOL_BINARY_BASE64: 'RDMOL_BINARY_BASE64',
  SMILES: 'SMILES',
  EXTENDED_SMILES: 'EXTENDED_SMILES',
  SMARTS: 'SMARTS',
  EXTENDED_SMARTS: 'EXTENDED_SMARTS',
  MDL_MOLV2000: 'MDL_MOLV2000',
  MDL_MOLV3000: 'MDL_MOLV3000',
  MAESTRO: 'MAESTRO',
  INCHI: 'INCHI',
  INCHI_KEY: 'INCHI_KEY',
  PDB: 'PDB',
  MOL2: 'MOL2',
  XYZ: 'XYZ',
  MRV: 'MRV',
  CDXML: 'CDXML',
  HELM: 'HELM',
  FASTA_PEPTIDE: 'FASTA_PEPTIDE',
  FASTA_DNA: 'FASTA_DNA',
  FASTA_RNA: 'FASTA_RNA',
  FASTA: 'FASTA',
} as const;
export type RepresentationFormat = (typeof RepresentationFormat)[keyof typeof RepresentationFormat];

// This is what the actual emscripten bound enum looks like, we use a "Tagged"
// type to prevent misuse of object literals and ensure only a Format is passed
// to the sketcher methods
type SketcherWasmFormat = Tagged<{ value: number }, 'SketcherWasmFormat'>;
type SketcherWasmImageFormat = Tagged<{ value: number }, 'SketcherWasmImageFormat'>;

export type SketcherWASM = {
  Format: { [K in keyof typeof RepresentationFormat]: SketcherWasmFormat };
  ImageFormat: { [K in keyof typeof SketcherImageFormat]: SketcherWasmImageFormat };
  sketcher_import_text: (text: string) => void;
  sketcher_export_text: (format: SketcherWasmFormat) => string;
  sketcher_export_image: (format: SketcherWasmImageFormat) => string;
  sketcher_clear: () => void;
  sketcher_is_empty: () => boolean;
  sketcher_has_monomers: () => boolean;
  sketcher_allow_monomeric: (allowMonomeric: boolean) => void;
  sketcher_changed_callback?: () => void;
  getExceptionMessage: (pointer: number) => string[];
};

export type SketcherRef = {
  /**
   * Provides a reference to the instance of the WASM sketcher application
   */
  getInstance: () => Promise<SketcherWASM>;
  sketcherExportStructure: (format: RepresentationFormat) => string | undefined;
  sketcherImportText: (representation: string, clearSketcher: boolean) => void;
};

export type SketcherProps = Merge<
  ComponentProps<'iframe'>,
  {
    representation?: string;
    onChange?: (sketcher: SketcherWASM) => void;
    onError?: (message: string) => void;
  }
>;

const Sketcher = forwardRef<SketcherRef, SketcherProps>(function Sketcher(props, ref) {
  const {
    representation,
    onChange,
    onError,
    height = '400',
    width = '100%',
    ...iframeProps
  } = props;
  const [sketcherInstance, setSketcherInstance] = useState<SketcherWASM | null>(null);

  const initializeSketcherRef = useCallback(async (iframe: HTMLIFrameElement | null) => {
    if (!iframe) {
      return;
    }

    let instance: SketcherWASM;
    while (!(instance = iframe.contentWindow?.Module)) {
      await new Promise((resolve) => setTimeout(resolve, 100));
    }
    setSketcherInstance(instance);
  }, []);

  const sketcherInstanceRef = useRef(sketcherInstance);
  sketcherInstanceRef.current = sketcherInstance;
  useImperativeHandle(
    ref,
    (): SketcherRef => ({
      getInstance: async () => {
        while (!sketcherInstanceRef.current) {
          await new Promise((resolve) => setTimeout(resolve, 100));
        }
        return sketcherInstanceRef.current;
      },
      sketcherExportStructure(format) {
        if (!sketcherInstanceRef.current) {
          return;
        }

        const { sketcher_export_text, Format } = sketcherInstanceRef.current;
        return sketcher_export_text(Format[format]);
      },
      sketcherImportText(representation, clearStructure = true) {
        if (!sketcherInstanceRef.current) {
          return;
        }
        const { sketcher_clear, sketcher_import_text } = sketcherInstanceRef.current;
        if (clearStructure) {
          sketcher_clear();
        }
        sketcher_import_text(representation);
      },
    }),
    [],
  );

  /**
   * Sync onChange callback from props to the sketcher's change callback
   */
  useEffect(() => {
    if (!sketcherInstance) {
      return;
    }

    sketcherInstance.sketcher_changed_callback = () => onChange?.(sketcherInstance);
    return () => {
      delete sketcherInstance.sketcher_changed_callback;
    };
  }, [onChange, sketcherInstance]);

  /**
   * Sync the representation from props with the contents of the sketcher canvas
   */
  useEffect(
    () => {
      if (!sketcherInstance) {
        return;
      }
      const { sketcher_clear, sketcher_import_text, getExceptionMessage } = sketcherInstance;

      sketcher_clear();
      if (!representation?.trim()) {
        return;
      }

      try {
        sketcher_import_text(representation);
      } catch (e) {
        if (typeof e === 'number') {
          onError?.(`Error importing to sketcher [${getExceptionMessage(e).join(': ')}]`);
        } else {
          onError?.(`Error importing to sketcher [${e}]`);
        }
      }
    },
    // omitting onError here because the error handler changing shouldn't re-trigger a sketcher
    // update, we just want to use the onError defined at the time of the representation change
    // render
    [representation, sketcherInstance],
  );

  return (
    <iframe
      ref={initializeSketcherRef}
      title="Maestro Sketcher"
      height={height}
      width={width}
      src={`static/wasm_shell.html?cache_bust=__WASM_HASH__`}
      {...iframeProps}
    />
  );
});

export default Sketcher;
