import { execFileSync } from 'node:child_process';
import { existsSync, mkdirSync, readdirSync, renameSync, rmSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const sketcherDir = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const e2eDir = path.resolve(sketcherDir, '..');
const snapshotRoot = path.join(e2eDir, '__snapshots__');
const archiveRoot = path.join(e2eDir, 'reference_archives');

if (!existsSync(snapshotRoot)) {
  throw new Error(`No extracted Playwright snapshots found at ${snapshotRoot}`);
}

mkdirSync(archiveRoot, { recursive: true });
const snapshots = readdirSync(snapshotRoot, { withFileTypes: true })
  .filter((entry) => entry.isDirectory())
  .map((entry) => entry.name)
  .sort();
const expectedArchives = new Set(snapshots.map((name) => `${name}.tar.gz`));

for (const name of snapshots) {
  // A separate deterministic archive per test keeps PR changes localized.
  const archivePath = path.join(archiveRoot, `${name}.tar.gz`);
  const temporaryPath = `${archivePath}.tmp`;
  execFileSync('tar', [
    '--sort=name',
    '--mtime=@0',
    '--owner=0',
    '--group=0',
    '--numeric-owner',
    '-C',
    snapshotRoot,
    '-czf',
    temporaryPath,
    name,
  ]);
  renameSync(temporaryPath, archivePath);
}

for (const entry of readdirSync(archiveRoot)) {
  if (entry.endsWith('.tar.gz') && !expectedArchives.has(entry)) {
    rmSync(path.join(archiveRoot, entry));
  }
}
