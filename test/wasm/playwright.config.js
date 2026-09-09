import { defineConfig } from '@playwright/test';

const buildDir = process.env.SKETCHER_WASM_BUILD_DIR || 'build';
const slowMo = Number(process.env.PLAYWRIGHT_SLOW_MO || 0);

export default defineConfig({
  // Expand the small, versioned Squish-derived structure/reference bundle
  // before test discovery.  The extracted files are ignored so PRs show one
  // intentional fixture update rather than a long list of data-file changes.
  globalSetup: './e2e/sketcher/scripts/prepare_runtime_fixtures.mjs',
  webServer: {
    command: `python3 -m http.server 8000 --directory ../../${buildDir}/sketcher_app`,
    url: 'http://localhost:8000',
    reuseExistingServer: !process.env.CI,
  },
  use: {
    baseURL: 'http://localhost:8000',
    headless: process.env.PLAYWRIGHT_HEADED !== '1',
    launchOptions: { slowMo },
    // Qt/WASM forwards Copy/Cut/Paste through the browser Clipboard API.  Give
    // every Playwright context both permissions before the page loads so
    // Chromium does not show its "wants to see text ..." permission prompt.
    permissions: ['clipboard-read', 'clipboard-write'],
    screenshot: 'only-on-failure',
    viewport: { width: 1280, height: 720 },
    actionTimeout: 10000,
  },
  retries: process.env.CI ? 1 : 0,
  snapshotPathTemplate: '{testDir}/{testFileDir}/__snapshots__/{testFileName}/{arg}{ext}',
  expect: {
    toMatchSnapshot: { maxDiffPixelRatio: 0.1 },
  },
  // Output to the same place as executable tests for GitHub actions
  reporter: [['junit', { outputFile: `../../${buildDir}/junit-report.xml` }]],
});
