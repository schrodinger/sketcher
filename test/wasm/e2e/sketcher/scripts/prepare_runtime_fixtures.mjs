import { createHash } from 'node:crypto';
import { existsSync, mkdirSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import extractSnapshotArchives from './extract_snapshot_archives.mjs';

const fromSquishDir = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const archivePath = path.join(
  fromSquishDir,
  'fixtures',
  'parity_runtime_fixtures.tar.gz',
);
const outputPath = path.join(fromSquishDir, '.runtime_fixtures');
const stampPath = path.join(outputPath, '.archive-sha256');
const archiveHash = createHash('sha256').update(readFileSync(archivePath)).digest('hex');

// This is intentionally a small runtime subset: browser import inputs and the
// text/SDF goldens used by the export parity test.  Historical Squish image
// references are not needed once Playwright's own visual baselines exist.
export default function prepareRuntimeFixtures() {
  extractSnapshotArchives();
  if (!existsSync(stampPath) || readFileSync(stampPath, 'utf8') !== archiveHash) {
    rmSync(outputPath, { force: true, recursive: true });
    mkdirSync(outputPath, { recursive: true });
    execFileSync('tar', ['-xzf', archivePath, '-C', outputPath]);
    writeFileSync(stampPath, archiveHash);
  }
}
