import { defineConfig } from '@playwright/test';

const buildDir = process.env.SKETCHER_WASM_BUILD_DIR || 'build';

export default defineConfig({
  webServer: {
    command: `python3 -m http.server 8000 --directory ../../${buildDir}/sketcher_app`,
    url: 'http://localhost:8000',
    reuseExistingServer: !process.env.CI,
  },
  use: {
    baseURL: 'http://localhost:8000',
    screenshot: 'only-on-failure',
  },
  snapshotPathTemplate: '{testDir}/{testFileDir}/__snapshots__/{testFileName}/{arg}{ext}',
  expect: {
    toMatchSnapshot: { maxDiffPixelRatio: 0.1 },
  },
  // Output to the same place as executable tests for GitHub actions
  reporter: [['junit', { outputFile: `../../${buildDir}/junit-report.xml` }]],
});
