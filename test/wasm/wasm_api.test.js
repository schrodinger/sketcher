const { test, expect } = require('@playwright/test');

test.describe('WASM Sketcher API', () => {

    test.beforeEach(async ({ page }) => {
        // Navigate to the page that loads the WASM module and
        // wait for the WASM module to be fully loaded and available
        await page.goto('/wasm_shell.html');
        // wait for the sketcher WASM module to be loaded
        await page.waitForFunction(() => typeof window.Module !== undefined);
    });

    test('Roundtrip SMILES through the sketcher', async ({ page }) => {
        const smiles_input = "c1ccccc1";
        const exportedText = await page.evaluate(async (text) => {
            Module.sketcher_import_text(text);
            return Module.sketcher_export_text(Module.Format.SMILES);
        }, smiles_input);
        // Verify that the exported text matches the input
        expect(exportedText).toBe(smiles_input);
    });

    test('Clear the sketcher', async ({ page }) => {
        // Confirm the sketcher starts empty
        const is_empty = await page.evaluate(() => Module.sketcher_is_empty());
        expect(is_empty).toBe(true);
        // Import a molecule and confirm it's no longer empty
        await page.evaluate(() => {
            Module.sketcher_import_text("C");
        });
        const is_empty = await page.evaluate(() => Module.sketcher_is_empty());
        expect(is_empty).toBe(false);
        // Clear the sketcher and confirm it's empty again
        await page.evaluate(() => {
            Module.sketcher_clear();
        });
        const is_cleared = await page.evaluate(() => Module.sketcher_is_empty());
        expect(is_cleared).toBe(true);
    });

    test('Export as image from the sketcher', async ({ page }) => {
        const smiles_input = "C=O";
        const base64Content = await page.evaluate((smiles) => {
            Module.sketcher_import_text(smiles);
            return Module.sketcher_export_image(Module.ImageFormat.SVG);
        }, smiles_input);

        // SVG validation using regex pattern similar to Python version
        const SVG_REGEX = /(?:<\?xml\b[^>]*>[^<]*)?(?:<!--.*?-->[^<]*)*(?:<svg|<!DOCTYPE svg)\b/s;
        function isValidSvg(svgData) {
            return SVG_REGEX.test(svgData);
        }

        // Decode base64 to get the actual SVG content
        const svgContent = Buffer.from(base64Content, 'base64').toString('utf8');
        expect(isValidSvg(svgContent)).toBe(true);
    });

});