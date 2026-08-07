import { readFileSync, readdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { expect, test } from '@playwright/test';
import { waitForSketcherReady, clearSketcher, tryImport } from './e2e_helpers.js';

// Ported from `suite_2D_sketcher/tst_problematic_structures`: import historically
// problematic structures and confirm each loads without throwing. Multi-structure
// SDF/.mae files are split into individual records (the importer takes one at a time).

const FIXTURES_DIR = join(
  dirname(fileURLToPath(import.meta.url)),
  'fixtures',
  'problematic_structures',
);

/**
 * Split a multi-record SDF into individual mol blocks (re-appending the `$$$$`
 * delimiter to each). The record content is left untrimmed on purpose: a mol
 * block's first three lines are the header (title/program/comment) and any of
 * them may be blank, so trimming would shift the counts line and break format
 * detection. Whitespace-only chunks (e.g. after the final delimiter) are
 * dropped.
 */
function splitSdfRecords(text) {
  return text
    .split(/\$\$\$\$\r?\n/)
    .filter((r) => /\S/.test(r))
    .map((r) => `${r}\n$$$$\n`);
}

/**
 * Split a multi-structure Maestro file into one importable string per CT block,
 * each paired with the leading m2io version header. Top-level `{...}` blocks are
 * found by brace matching; the first is the version header, and every following
 * `f_m_ct` block becomes its own structure.
 */
function splitMaeStructures(text) {
  const blocks = [];
  let depth = 0;
  let start = -1;
  let prevEnd = 0;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (c === '{') {
      if (depth === 0) start = i;
      depth++;
    } else if (c === '}') {
      depth--;
      if (depth === 0) {
        blocks.push({ label: text.slice(prevEnd, start), body: text.slice(start, i + 1) });
        prevEnd = i + 1;
      }
    }
  }
  if (blocks.length < 2) {
    return [text];
  }
  const header = blocks[0].label + blocks[0].body;
  return blocks
    .slice(1)
    .filter((b) => /f_m_ct/.test(b.label))
    .map((b) => `${header}\n${b.label.replace(/^[\s\S]*?(?=f_m_ct)/, '')}${b.body}\n`);
}

/** Build the per-file list of structures to import. */
function structuresFor(fileName) {
  const text = readFileSync(join(FIXTURES_DIR, fileName), 'utf8');
  return fileName.endsWith('.sdf') ? splitSdfRecords(text) : splitMaeStructures(text);
}

const FIXTURE_FILES = readdirSync(FIXTURES_DIR).filter(
  (f) => f.endsWith('.sdf') || f.endsWith('.mae'),
);
const STRUCTURES_PER_TEST = 10;

function chunkStructures(structures) {
  const chunks = [];
  for (let i = 0; i < structures.length; i += STRUCTURES_PER_TEST) {
    chunks.push({
      startIndex: i,
      structures: structures.slice(i, i + STRUCTURES_PER_TEST),
    });
  }
  return chunks;
}

function testName(fileName, startIndex, count, total) {
  if (total === 1) {
    return `import ${fileName}`;
  }
  return `import ${fileName} records ${startIndex + 1}-${startIndex + count}`;
}

test.beforeEach(async ({ page }) => {
  await waitForSketcherReady(page);
});

test.describe('Problematic structures load without crashing', () => {
  for (const fileName of FIXTURE_FILES) {
    const allStructures = structuresFor(fileName);
    for (const { startIndex, structures } of chunkStructures(allStructures)) {
      test(
        testName(fileName, startIndex, structures.length, allStructures.length),
        async ({ page }) => {
          const failures = [];

          for (let i = 0; i < structures.length; i++) {
            await clearSketcher(page);
            const result = await tryImport(page, structures[i]);
            const recordNumber = startIndex + i + 1;
            if (!result.ok) {
              failures.push(`record ${recordNumber}: threw: ${result.error}`);
            } else if (result.empty) {
              failures.push(`record ${recordNumber}: imported nothing (empty canvas)`);
            }
          }

          expect(
            failures,
            `${failures.length}/${structures.length} structures failed in ${fileName}:\n` +
              failures.join('\n'),
          ).toEqual([]);
        },
      );
    }
  }
});
