const { test, expect } = require("@playwright/test");

test.describe("WASM Sketcher API", () => {
  test.beforeEach(async ({ page }) => {
    // Navigate to the page that loads the WASM module and
    // wait for the WASM module to be fully loaded and available
    await page.goto("/wasm_shell.html");
    // wait for the sketcher WASM module to be loaded
    await page.waitForFunction(() => typeof window.Module !== undefined);
  });

  // Test import and export for all formats
  const allFormats = [
    // AUTO_DETECT doesn't make sense to test here
    { name: "RDMOL_BINARY_BASE64" },
    { name: "SMILES" },
    { name: "EXTENDED_SMILES" },
    { name: "SMARTS" },
    { name: "EXTENDED_SMARTS" },
    // MDL_MOLV2000 is import only
    { name: "MDL_MOLV3000" },
    { name: "MAESTRO" },
    { name: "INCHI" },
    // INCHI_KEY is export only
    { name: "PDB" },
    // MOL2 is import only
    { name: "XYZ" },
    { name: "MRV" },
    // CDXML does not import correctly in WASM builds
    { name: "HELM" },
    // FASTA_PEPTIDE, FASTA_DNA, FASTA_RNA are import-only formats
    // FASTA is an export-only format
    // FMP and CUSTOM_ENTITY are not supported
  ];

  allFormats.forEach(({ name, expectedError }) => {
    test(`Import/Export Format ${name}`, async ({ page }) => {
      // Export a molecule in the specified format
      const exported = await page.evaluate(
        (data) => {
          const { formatName } = data;
          Module.sketcher_clear();
          // Use alanine - works for both atomistic and biologics formats
          Module.sketcher_import_text("C[C@H](N)C=O");
          const exported = Module.sketcher_export_text(
            Module.Format[formatName],
          );
          return exported;
        },
        { formatName: name },
      );

      expect(exported).toBeTruthy();
      expect(typeof exported).toBe("string");
      expect(exported.length).toBeGreaterThan(0);

      // Roundtrip it back in
      await page.evaluate(
        (data) => {
          const { exportedText } = data;
          Module.sketcher_clear();
          Module.sketcher_import_text(exportedText);
        },
        { exportedText: exported },
      );
    });
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

  test("Check if molecule has monomers", async ({ page }) => {
    // Empty sketcher should return false
    const empty_has_monomers = await page.evaluate(() =>
      Module.sketcher_has_monomers(),
    );
    expect(empty_has_monomers).toBe(false);

    // SMILES should return false
    await page.evaluate(() => {
      Module.sketcher_import_text("c1ccccc1");
    });
    const smiles_has_monomers = await page.evaluate(() =>
      Module.sketcher_has_monomers(),
    );
    expect(smiles_has_monomers).toBe(false);

    // HELM should return true
    await page.evaluate(() => {
      Module.sketcher_clear();
      Module.sketcher_import_text("PEPTIDE1{A.S.D.F.G.H.W}$$$$V2.0");
    });
    const helm_has_monomers = await page.evaluate(() =>
      Module.sketcher_has_monomers(),
    );
    expect(helm_has_monomers).toBe(true);

    // After clearing, should return false again
    await page.evaluate(() => {
      Module.sketcher_clear();
    });
    const cleared_has_monomers = await page.evaluate(() =>
      Module.sketcher_has_monomers(),
    );
    expect(cleared_has_monomers).toBe(false);
  });

  // Parameterized test for image export formats
  const imageFormats = ["SVG", "PNG"];

  imageFormats.forEach((imageFormat) => {
    test(`Export as ${imageFormat} image from the sketcher`, async ({
      page,
    }) => {
      const smiles_input = "C=O";
      const base64Content = await page.evaluate(
        (data) => {
          const { smiles, format } = data;
          Module.sketcher_import_text(smiles);
          return Module.sketcher_export_image(Module.ImageFormat[format]);
        },
        { smiles: smiles_input, format: imageFormat },
      );

      // Verify we got non-empty base64 content
      expect(base64Content).toBeTruthy();
      expect(typeof base64Content).toBe("string");
      expect(base64Content.length).toBeGreaterThan(0);

      if (imageFormat === "SVG") {
        const SVG_REGEX =
          /(?:<\?xml\b[^>]*>[^<]*)?(?:<!--.*?-->[^<]*)*(?:<svg|<!DOCTYPE svg)\b/s;
        function isValidSvg(svgData) {
          return SVG_REGEX.test(svgData);
        }
        // Decode base64 to get the actual SVG content
        const svgContent = Buffer.from(base64Content, "base64").toString(
          "utf8",
        );
        expect(isValidSvg(svgContent)).toBe(true);
      } else if (imageFormat === "PNG") {
        const pngBuffer = Buffer.from(base64Content, "base64");
        expect(pngBuffer.length).toBeGreaterThan(8);
        // PNG files start with specific magic bytes (89 50 4E 47 0D 0A 1A 0A)
        const pngSignature = [0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a];
        for (let i = 0; i < pngSignature.length; i++) {
          expect(pngBuffer[i]).toBe(pngSignature[i]);
        }
      }
    });
  });
});
