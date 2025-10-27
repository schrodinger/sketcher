import { defineConfig } from '@playwright/test';

export default defineConfig({
  webServer: {
    command: 'python3 -m http.server 8000 --directory ../../sketcher-wasm',
    url: 'http://localhost:8000',
    reuseExistingServer: !process.env.CI,
  },
  use: {
    baseURL: 'http://localhost:8000',
    screenshot: 'only-on-failure',
  },
  snapshotPathTemplate: '{testDir}/{testFileDir}/__snapshots__/{testFileName}/{arg}{ext}',
  // Output to the same place as executable tests for GitHub actions
  reporter: [['junit', { outputFile: '../../build/junit-report.xml' }]],
});
