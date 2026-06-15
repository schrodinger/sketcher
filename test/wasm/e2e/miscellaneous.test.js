import { readFileSync, readdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { expect, test } from '@playwright/test';
import { waitForSketcherReady, clearSketcher, tryImport } from './e2e_helpers.js';

// Ported from `suite_2D_sketcher_new/tst_miscellaneous`, which imported the
// SHARED-6974 regression structures (one.mae, two.mae) and saved a reference
// image of each. The pixel comparison isn't portable to the WASM renderer, so
// we keep the regression intent: each structure must import and render without
// error.

const FIXTURES_DIR = join(
  dirname(fileURLToPath(import.meta.url)),
  'fixtures',
  'miscellaneous_structures',
);

const MAE_FILES = readdirSync(FIXTURES_DIR).filter((f) => f.endsWith('.mae'));

test.beforeEach(async ({ page }) => {
  await waitForSketcherReady(page);
});

test.describe('Miscellaneous regression structures (SHARED-6974)', () => {
  for (const fileName of MAE_FILES) {
    test(`import ${fileName}`, async ({ page }) => {
      const text = readFileSync(join(FIXTURES_DIR, fileName), 'utf8');
      await clearSketcher(page);
      const result = await tryImport(page, text);
      expect(result.ok, `${fileName} import threw: ${result.error}`).toBe(true);
      expect(result.empty, `${fileName} produced an empty canvas`).toBe(false);
    });
  }
});
