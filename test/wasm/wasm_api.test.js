import { expect, test } from '@playwright/test';

test.beforeEach(async ({ page }) => {
  // Navigate to the page that loads the WASM module and
  // wait for the WASM module to be fully loaded and available
  await page.goto('/wasm_shell.html');
  await page.waitForFunction(() => typeof window.Module !== 'undefined', {
    timeout: 20000,
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
    {
      format: 'CDXML',
      skip: [true, "Format doesn't import correctly in WASM builds"],
    },
    { format: 'HELM' },
    { format: 'FASTA_PEPTIDE', exportUnsupported: true },
    { format: 'FASTA_DNA', exportUnsupported: true },
    { format: 'FASTA_RNA', exportUnsupported: true },
    { format: 'FASTA', importUnsupported: true },
  ];

  FORMATS.forEach(({ format, skip, exportUnsupported, importUnsupported }) => {
    test(`exporting/importing ${format}`, async ({ page }) => {
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

      // Test importing the exported value (round-trip test)
      if (!importUnsupported) {
        const importSuccessful = await page.evaluate((exportedText) => {
          Module.sketcher_clear();
          Module.sketcher_allow_monomeric();
          Module.sketcher_import_text(exportedText);
          return !Module.sketcher_is_empty();
        }, exportedText);
        expect(importSuccessful).toBe(true);
      }
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
      Module.sketcher_allow_monomeric();
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

  test('Import CDXML from ketcher file', async ({ page }) => {
    const cdxmlInput = `<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd">
<CDXML BondLength="30.000000" LabelFont="3" CaptionFont="4">
    <fonttable id="1">
        <font id="2" charset="utf-8" name="Arial"/>
        <font id="3" charset="utf-8" name="Times New Roman"/>
    </fonttable>
    <colortable>
        <color r="1" g="1" b="1"/>
        <color r="0" g="0" b="0"/>
        <color r="1" g="0" b="0"/>
        <color r="1" g="1" b="0"/>
        <color r="0" g="1" b="0"/>
        <color r="0" g="1" b="1"/>
        <color r="0" g="0" b="1"/>
        <color r="1" g="0" b="1"/>
        <color r="0.5" g="0.5" b="0.5"/>
    </colortable>
    <page HeightPages="1" WidthPages="1" BoundingBox="240.295486 161.250992 292.204529 101.249008">
        <fragment id="4">
            <n id="5" p="0.000000 15.003225"/>
            <n id="6" p="51.909027 14.988670"/>
            <n id="7" p="26.003637 -0.000000"/>
            <n id="8" p="51.909027 45.016960"/>
            <n id="9" p="0.000000 45.151600"/>
            <n id="10" p="26.069126 60.001980"/>
            <b id="11" B="7" E="5" Order="2"/>
            <b id="12" B="5" E="9"/>
            <b id="13" B="9" E="10" Order="2"/>
            <b id="14" B="10" E="8"/>
            <b id="15" B="8" E="6" Order="2"/>
            <b id="16" B="6" E="7"/>
        </fragment>
    </page>
</CDXML>`;

    // Import the CDXML and extract SMILES
    const exportedSmiles = await page.evaluate((cdxml) => {
      Module.sketcher_import_text(cdxml);
      return Module.sketcher_export_text(Module.Format.SMILES);
    }, cdxmlInput);

    expect(exportedSmiles).toBe('C1=CC=CC=C1');
  });
});
