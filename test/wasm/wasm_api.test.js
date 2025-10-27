import { expect, test } from '@playwright/test';

test.describe('WASM Sketcher API', () => {
  test.beforeEach(async ({ page }) => {
    // Navigate to the page that loads the WASM module and
    // wait for the WASM module to be fully loaded and available
    await page.goto('/wasm_shell.html');
    await page.waitForFunction(() => typeof window.Module !== 'undefined');
  });

  // Test import and export for all formats
  const FORMATS = [
    { format: 'AUTO_DETECT', skip: [true, "Doesn't make sense to test here"] },
    { format: 'RDMOL_BINARY_BASE64' },
    { format: 'SMILES' },
    { format: 'EXTENDED_SMILES' },
    { format: 'SMARTS' },
    { format: 'EXTENDED_SMARTS' },
    { format: 'MDL_MOLV2000', exportUnsupported: true },
    { format: 'MDL_MOLV3000' },
    { format: 'MAESTRO' },
    { format: 'INCHI' },
    { format: 'INCHI_KEY', importUnsupported: true },
    { format: 'PDB' },
    { format: 'MOL2', exportUnsupported: true },
    { format: 'XYZ' },
    { format: 'MRV' },
    { format: 'CDXML', skip: [true, "Format doesn't import correctly in WASM builds"] },
    { format: 'HELM' },
    { format: 'FASTA_PEPTIDE', exportUnsupported: true },
    { format: 'FASTA_DNA', exportUnsupported: true },
    { format: 'FASTA_RNA', exportUnsupported: true },
    { format: 'FASTA', importUnsupported: true },
    { format: 'FMP', skip: [true, 'Format not supported'] },
    { format: 'CUSTOM_ENTITY', skip: [true, 'Format not supported'] },
  ];

  FORMATS.forEach(({ format, skip, exportUnsupported }) => {
    test(`importing SMILES and exporting ${format}`, async ({ page }) => {
      if (skip) {
        test.skip(...skip);
      } else {
        test.skip(!!exportUnsupported, `${format} is import only`);
      }
      const exportedText = await page.evaluate((format) => {
        Module.sketcher_clear();
        Module.sketcher_import_text('C[C@H](N)C=O');
        const exported = Module.sketcher_export_text(Module.Format[format]);
        return exported;
      }, format);
      expect(exportedText).toMatchSnapshot(`export_text_${format}.txt`);
    });
  });

  test('clearing the sketcher', async ({ page }) => {
    const isEmptyOnLoad = await page.evaluate(() => Module.sketcher_is_empty());
    expect(isEmptyOnLoad).toBe(true);

    const isEmptyAfterImport = await page.evaluate(() => {
      Module.sketcher_import_text('C');
      return Module.sketcher_is_empty();
    });
    expect(isEmptyAfterImport).toBe(false);

    const isEmptyAfterClear = await page.evaluate(() => {
      Module.sketcher_clear();
      return Module.sketcher_is_empty();
    });
    expect(isEmptyAfterClear).toBe(true);
  });

  test('checking if molecule has monomers', async ({ page }) => {
    const hasMonomersOnLoad = await page.evaluate(() => Module.sketcher_has_monomers());
    expect(hasMonomersOnLoad).toBe(false);

    const hasMonomersAfterSmilesImport = await page.evaluate(() => {
      Module.sketcher_import_text('c1ccccc1');
      return Module.sketcher_has_monomers();
    });
    expect(hasMonomersAfterSmilesImport).toBe(false);

    const hasMonomersAfterHelmImport = await page.evaluate(() => {
      Module.sketcher_clear();
      Module.sketcher_import_text('PEPTIDE1{A.S.D.F.G.H.W}$$$$V2.0');
      return Module.sketcher_has_monomers();
    });
    expect(hasMonomersAfterHelmImport).toBe(true);

    const hasMonomersAfterClear = await page.evaluate(() => {
      Module.sketcher_clear();
      return Module.sketcher_has_monomers();
    });
    expect(hasMonomersAfterClear).toBe(false);
  });

  // Test image export for all formats
  ['SVG', 'PNG'].forEach((imageFormat) => {
    test(`exporting a ${imageFormat} image`, async ({ page }) => {
      const base64Content = await page.evaluate((imageFormat) => {
        Module.sketcher_import_text('C=O');
        return Module.sketcher_export_image(Module.ImageFormat[imageFormat]);
      }, imageFormat);
      const buffer = Buffer.from(base64Content, 'base64');

      expect(buffer).toMatchSnapshot(`export_image.${imageFormat.toLowerCase()}`);
    });
  });
});
