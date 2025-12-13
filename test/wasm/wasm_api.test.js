import { expect, test } from '@playwright/test';
import fs from 'fs';
import path from 'path';

test.beforeEach(async ({ page }) => {
  // Navigate to the page that loads the WASM module and
  // wait for the WASM module to be fully loaded and available
  await page.goto('/wasm_shell.html');
  await page.waitForFunction(() => typeof window.Module !== 'undefined', {
    timeout: 20000,
  });
  // Enable monomer support in tests
  await page.evaluate(() => {
    Module.sketcher_allow_monomeric();
  });
});

test.describe('WASM Sketcher UI', () => {
  test('sketcher UI loads correctly', async ({ page }) => {
    expect(await page.screenshot()).toMatchSnapshot('sketcher-load.png');
  });
});

test.describe('WASM Sketcher API', () => {
  // Test import and export for all formats
  const FORMATS = [
    { format: 'AUTO_DETECT', importUnsupported: true, exportUnsupported: true },
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
    { format: 'CDXML', exportUnsupported: true },
    { format: 'HELM' },
    { format: 'FASTA_PEPTIDE', exportUnsupported: true },
    { format: 'FASTA_DNA', exportUnsupported: true },
    { format: 'FASTA_RNA', exportUnsupported: true },
    { format: 'FASTA', importUnsupported: true, exportUnsupported: true },
  ];

  // Test exporting for formats that support export
  FORMATS.filter(({ exportUnsupported }) => !exportUnsupported).forEach(({ format }) => {
    test(`exporting ${format}`, async ({ page }) => {
      const exportedText = await page.evaluate((format) => {
        Module.sketcher_clear();
        Module.sketcher_import_text('C[C@H](N)C=O');
        const exported = Module.sketcher_export_text(Module.Format[format]);
        return exported;
      }, format);
      expect(exportedText).toMatchSnapshot(`text-${format}.txt`);
    });
  });

  // Test importing for formats that support import
  FORMATS.filter(({ importUnsupported }) => !importUnsupported).forEach(({ format }) => {
    test(`importing ${format}`, async ({ page }, testInfo) => {
      const fixturePath = path.join(
        __dirname,
        '..',
        'testfiles',
        'autodetect',
        `text-${format.replace(/_/g, '-')}.txt`,
      );
      const inputText = fs.readFileSync(fixturePath, 'utf8');

      const importSuccessful = await page.evaluate(
        ({ inputText, format }) => {
          // Import with AUTO_DETECT
          Module.sketcher_clear();
          Module.sketcher_import_text(inputText);
          const autoDetectWorks = !Module.sketcher_is_empty();
          // Import with explicit format parameter
          Module.sketcher_clear();
          Module.sketcher_import_text(inputText, Module.Format[format]);
          const explicitFormatWorks = !Module.sketcher_is_empty();
          return { autoDetectWorks, explicitFormatWorks };
        },
        { inputText, format },
      );
      expect(importSuccessful.autoDetectWorks).toBe(true);
      expect(importSuccessful.explicitFormatWorks).toBe(true);
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

      let actualImage = buffer;
      // For SVGs, we render them in the browser and take a screenshot so we can do actual visual
      // comparison instead of comparing the markup (which might differ slightly by platform)
      if (imageFormat === 'SVG') {
        const svgDataUri = `data:image/svg+xml;charset=utf-8;base64,${buffer.toString('base64')}`;
        await page.goto(svgDataUri);
        actualImage = await page
          .locator('svg')
          .screenshot({ type: 'png', omitBackground: true, scale: 'css' });
      }

      expect(actualImage).toMatchSnapshot(`export-image-${imageFormat}.png`);
    });
  });
});
