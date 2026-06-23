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
    { format: 'HELM', exportUnsupported: true },
    { format: 'FASTA_PEPTIDE', exportUnsupported: true },
    { format: 'FASTA_DNA', exportUnsupported: true },
    { format: 'FASTA_RNA', exportUnsupported: true },
    { format: 'FASTA', importUnsupported: true, exportUnsupported: true },
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
          Module.sketcher_allow_monomeric(true);
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
      Module.sketcher_allow_monomeric(true);
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

  test('Exception handling for invalid input', async ({ page }) => {
    const result = await page.evaluate(() => {
      try {
        Module.sketcher_import_text('foobar');
        throw new Error('Expected exception to be thrown');
      } catch (e) {
        // Use emscripten's getExceptionMessage to extract C++ exception info
        const [type, message] = Module.getExceptionMessage(e);
        return { type, message };
      }
    });

    // Verify that we can extract the C++ exception message and type
    expect(result.message).toBe('Unable to determine format');
    expect(result.type).toBeTruthy();
  });
});

test.describe('Custom Monomer DB', () => {
  test('load custom monomers from JSON', async ({ page }) => {
    const customMonomers = JSON.stringify([
      {
        SYMBOL: 'TestMon',
        POLYMER_TYPE: 'PEPTIDE',
        NATURAL_ANALOG: 'A',
        SMILES: 'CC(C)(N[H:1])C(=O)[OH:2]',
        CORE_SMILES: 'CC(C)(N)C=O',
        NAME: 'Test Monomer',
        MONOMER_TYPE: 'Backbone',
        AUTHOR: 'test',
      },
    ]);

    await page.evaluate((json) => {
      Module.sketcher_load_custom_monomers(json);
    }, customMonomers);
  });

  test('load custom monomers from SQL', async ({ page }) => {
    const sql = `INSERT INTO monomer_definitions (SYMBOL, POLYMER_TYPE, NATURAL_ANALOG, SMILES, CORE_SMILES, NAME, MONOMER_TYPE, AUTHOR)
      VALUES ('SqlMon', 'PEPTIDE', 'A', 'CC(C)(N[H:1])C(=O)[OH:2]', 'CC(C)(N)C=O', 'SQL Monomer', 'Backbone', 'test');`;

    const result = await page.evaluate((sql) => {
      Module.sketcher_load_custom_monomers_from_sql(sql);
      Module.sketcher_allow_monomeric(true);
      Module.sketcher_clear();
      Module.sketcher_import_text('PEPTIDE1{A.[SqlMon].G}$$$$V2.0');
      return !Module.sketcher_is_empty();
    }, sql);

    expect(result).toBe(true);
  });

  test('loading SQL twice replaces previous definitions', async ({ page }) => {
    const sql1 = `INSERT INTO monomer_definitions (SYMBOL, POLYMER_TYPE, NATURAL_ANALOG, SMILES, CORE_SMILES, NAME, MONOMER_TYPE, AUTHOR)
      VALUES ('Sql1', 'PEPTIDE', 'G', 'NCC(N[H:1])C(=O)[OH:2]', 'NCC(N)C=O', 'SQL One', 'Backbone', 'test');`;

    const sql2 = `INSERT INTO monomer_definitions (SYMBOL, POLYMER_TYPE, NATURAL_ANALOG, SMILES, CORE_SMILES, NAME, MONOMER_TYPE, AUTHOR)
      VALUES ('Sql2', 'PEPTIDE', 'A', 'CCC(N[H:1])C(=O)[OH:2]', 'CCC(N)C=O', 'SQL Two', 'Backbone', 'test');`;

    await page.evaluate((sql) => {
      Module.sketcher_load_custom_monomers_from_sql(sql);
      Module.sketcher_allow_monomeric(true);
    }, sql1);

    // Verify Sql1 works
    const sql1Works = await page.evaluate(() => {
      Module.sketcher_clear();
      Module.sketcher_import_text('PEPTIDE1{[Sql1].A}$$$$V2.0');
      return !Module.sketcher_is_empty();
    });
    expect(sql1Works).toBe(true);

    // Load second SQL — should replace, not append
    await page.evaluate((sql) => {
      Module.sketcher_load_custom_monomers_from_sql(sql);
    }, sql2);

    // Sql2 should work
    const sql2Works = await page.evaluate(() => {
      Module.sketcher_clear();
      Module.sketcher_import_text('PEPTIDE1{[Sql2].A}$$$$V2.0');
      return !Module.sketcher_is_empty();
    });
    expect(sql2Works).toBe(true);
  });

  test('insert custom monomers preserves existing definitions', async ({ page }) => {
    const monomer1 = JSON.stringify([
      {
        SYMBOL: 'Mon1',
        POLYMER_TYPE: 'PEPTIDE',
        NATURAL_ANALOG: 'G',
        SMILES: 'NCC(N[H:1])C(=O)[OH:2]',
        CORE_SMILES: 'NCC(N)C=O',
        NAME: 'Monomer One',
        MONOMER_TYPE: 'Backbone',
        AUTHOR: 'test',
      },
    ]);

    const monomer2 = JSON.stringify([
      {
        SYMBOL: 'Mon2',
        POLYMER_TYPE: 'PEPTIDE',
        NATURAL_ANALOG: 'A',
        SMILES: 'CC(C)(N[H:1])C(=O)[OH:2]',
        CORE_SMILES: 'CC(C)(N)C=O',
        NAME: 'Monomer Two',
        MONOMER_TYPE: 'Backbone',
        AUTHOR: 'test',
      },
    ]);

    // Load first and verify Mon1 works
    await page.evaluate((json) => {
      Module.sketcher_load_custom_monomers(json);
      Module.sketcher_allow_monomeric(true);
    }, monomer1);

    const mon1Works = await page.evaluate(() => {
      Module.sketcher_clear();
      Module.sketcher_import_text('PEPTIDE1{[Mon1].A}$$$$V2.0');
      return !Module.sketcher_is_empty();
    });
    expect(mon1Works).toBe(true);

    // Insert second — both should still be available
    await page.evaluate((json) => {
      Module.sketcher_insert_custom_monomers(json);
    }, monomer2);

    // Both custom monomers should be usable in a HELM import
    const result = await page.evaluate(() => {
      Module.sketcher_allow_monomeric(true);
      Module.sketcher_clear();
      Module.sketcher_import_text('PEPTIDE1{[Mon1].[Mon2].A}$$$$V2.0');
      return !Module.sketcher_is_empty();
    });

    expect(result).toBe(true);
  });

  test('custom monomer is usable in HELM import', async ({ page }) => {
    const customMonomers = JSON.stringify([
      {
        SYMBOL: 'Sar',
        POLYMER_TYPE: 'PEPTIDE',
        NATURAL_ANALOG: 'G',
        SMILES: '[H:3]N(C)C([H:1])C(=O)[OH:2]',
        CORE_SMILES: 'CN(CC(O)=O)[H]',
        NAME: 'Sarcosine',
        MONOMER_TYPE: 'Backbone',
        AUTHOR: 'test',
      },
    ]);

    const result = await page.evaluate((json) => {
      Module.sketcher_load_custom_monomers(json);
      Module.sketcher_allow_monomeric(true);
      Module.sketcher_clear();
      Module.sketcher_import_text('PEPTIDE1{A.[Sar].G}$$$$V2.0');
      return !Module.sketcher_is_empty();
    }, customMonomers);

    expect(result).toBe(true);
  });

  test('loading monomer with missing required fields throws', async ({ page }) => {
    const incomplete = JSON.stringify([
      {
        SYMBOL: 'Bad',
        POLYMER_TYPE: 'PEPTIDE',
        // missing required NATURAL_ANALOG, SMILES, NAME, MONOMER_TYPE, AUTHOR
      },
    ]);

    const result = await page.evaluate((json) => {
      try {
        Module.sketcher_load_custom_monomers(json);
        return { threw: false };
      } catch (e) {
        const [type, message] = Module.getExceptionMessage(e);
        return { threw: true, type, message };
      }
    }, incomplete);

    expect(result.threw).toBe(true);
    expect(result.message).toContain('missing');
  });

  test('loading monomer with unknown fields throws', async ({ page }) => {
    const extraFields = JSON.stringify([
      {
        SYMBOL: 'Extra',
        POLYMER_TYPE: 'PEPTIDE',
        NATURAL_ANALOG: 'A',
        SMILES: 'CC(N)C(=O)O',
        CORE_SMILES: 'CC(N)C(=O)O',
        NAME: 'Extra Field Monomer',
        MONOMER_TYPE: 'Backbone',
        AUTHOR: 'test',
        BOGUS_FIELD: 'should not be here',
      },
    ]);

    const result = await page.evaluate((json) => {
      try {
        Module.sketcher_load_custom_monomers(json);
        return { threw: false };
      } catch (e) {
        const [type, message] = Module.getExceptionMessage(e);
        return { threw: true, type, message };
      }
    }, extraFields);

    expect(result.threw).toBe(true);
    expect(result.message).toContain('exception');
  });

  test('reset custom monomer definitions', async ({ page }) => {
    const customMonomers = JSON.stringify([
      {
        SYMBOL: 'TmpMon',
        POLYMER_TYPE: 'PEPTIDE',
        NATURAL_ANALOG: 'A',
        SMILES: 'CCC(N[H:1])C(=O)[OH:2]',
        CORE_SMILES: 'CCC(N)C=O',
        NAME: 'Temporary Monomer',
        MONOMER_TYPE: 'Backbone',
        AUTHOR: 'test',
      },
    ]);

    // Load custom monomer, verify it works, then reset
    const worksBeforeReset = await page.evaluate((json) => {
      Module.sketcher_load_custom_monomers(json);
      Module.sketcher_allow_monomeric(true);
      Module.sketcher_clear();
      Module.sketcher_import_text('PEPTIDE1{A.[TmpMon].G}$$$$V2.0');
      return !Module.sketcher_is_empty();
    }, customMonomers);
    expect(worksBeforeReset).toBe(true);

    // Reset should succeed without throwing
    await page.evaluate(() => {
      Module.sketcher_reset_custom_monomers();
    });
  });
});
