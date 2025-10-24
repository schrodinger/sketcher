const { test, expect } = require("@playwright/test");

test.describe("WASM Sketcher API", () => {
  test.beforeEach(async ({ page }) => {
    // Navigate to the page that loads the WASM module and
    // wait for the WASM module to be fully loaded and available
    await page.goto("/wasm_shell.html");
    // wait for the sketcher WASM module to be loaded
    await page.waitForFunction(() => typeof window.Module !== undefined);
  });

  test("Roundtrip SMILES through the sketcher", async ({ page }) => {
    const smiles_input = "c1ccccc1";
    const exportedText = await page.evaluate(async (text) => {
      Module.sketcher_import_text(text);
      return Module.sketcher_export_text(Module.Format.SMILES);
    }, smiles_input);
    // Verify that the exported text matches the input
    expect(exportedText).toBe(smiles_input);
  });

  test("Clear the sketcher", async ({ page }) => {
    // Confirm the sketcher starts empty
    const is_empty = await page.evaluate(() => Module.sketcher_is_empty());
    expect(is_empty).toBe(true);
    // Import a molecule and confirm it's no longer empty
    await page.evaluate(() => {
      Module.sketcher_import_text("C");
    });
    const has_content = await page.evaluate(() => !Module.sketcher_is_empty());
    expect(has_content).toBe(true);
    // Clear the sketcher and confirm it's empty again
    await page.evaluate(() => {
      Module.sketcher_clear();
    });
    const is_cleared = await page.evaluate(() => Module.sketcher_is_empty());
    expect(is_cleared).toBe(true);
  });

  test("Export as image from the sketcher", async ({ page }) => {
    const smiles_input = "C=O";
    const base64Content = await page.evaluate((smiles) => {
      Module.sketcher_import_text(smiles);
      return Module.sketcher_export_image(Module.ImageFormat.SVG);
    }, smiles_input);

    // SVG validation using regex pattern similar to Python version
    const SVG_REGEX =
      /(?:<\?xml\b[^>]*>[^<]*)?(?:<!--.*?-->[^<]*)*(?:<svg|<!DOCTYPE svg)\b/s;
    function isValidSvg(svgData) {
      return SVG_REGEX.test(svgData);
    }

    // Decode base64 to get the actual SVG content
    const svgContent = Buffer.from(base64Content, "base64").toString("utf8");
    expect(isValidSvg(svgContent)).toBe(true);
  });

  test("Import CDXML from ketcher file", async ({ page }) => {
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

    // Import the CDXML, verify it's not empty, and extract SMILES
    const exportedSmiles = await page.evaluate((cdxml) => {
      Module.sketcher_import_text(cdxml);
      return Module.sketcher_export_text(Module.Format.SMILES);
    }, cdxmlInput);

    const is_empty = await page.evaluate(() => Module.sketcher_is_empty());
    expect(is_empty).toBe(false);
    expect(exportedSmiles).toBe("C1=CC=CC=C1");
  });
});
