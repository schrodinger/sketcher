import { execFileSync } from 'node:child_process';
import { existsSync, mkdirSync, readdirSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const sketcherDir = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const e2eDir = path.resolve(sketcherDir, '..');
const archiveRoot = path.join(e2eDir, 'reference_archives');
const snapshotRoot = path.join(e2eDir, '__snapshots__');

export default function extractSnapshotArchives() {
  if (!existsSync(archiveRoot)) return;
  mkdirSync(snapshotRoot, { recursive: true });
  for (const archive of readdirSync(archiveRoot).filter((name) => name.endsWith('.tar.gz')).sort()) {
    execFileSync('tar', ['-xzf', path.join(archiveRoot, archive), '-C', snapshotRoot]);
  }
}
