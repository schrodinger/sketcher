import { expect, test } from '@playwright/test';
import { readFile } from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { Sketcher } from '../wrappers/sketcher.js';

const FIXTURES_DIR = path.resolve(
  path.dirname(fileURLToPath(import.meta.url)),
  '../.runtime_fixtures/export_menu_structures',
);
const REFERENCES_DIR = path.resolve(
  path.dirname(fileURLToPath(import.meta.url)),
  '../.runtime_fixtures/export_menu_references',
);

// The source uses item 10 for image export. The original 2,000+ entry
// enumeration-library corpus is not retained locally, so the five later
// source fragments are reconstructed from the original Squish import SDF
// references, preserving the exact molecular graphs at indices below.
const IMAGE_SOURCE_FIXTURE = '10501119.sdf';
const SOURCE_FRAGMENT_IDS = [1710, 2421, 1406, 2261, 2295];

const TEXT_EXPORTS = [
  ['sdf', 'SDF', 'sdf', 'V3000'],
  ['smi', 'SMILES', 'smi', null],
  ['cxsmi', 'CXSMI', 'cxsmi', null],
  ['inchi', 'InChI', 'inchi', 'InChI=1S/'],
  ['pdb', 'PDB', 'pdb', 'HETATM'],
];

const IMAGE_EXPORTS = [
  ['save_image_height_width_png_1', 'PNG', 420, 690, false],
  ['save_image_height_width_png_2', 'PNG', 690, 420, false],
  ['transparent_background_png_1', 'PNG', 400, 400, true],
  ['transparent_background_png_2', 'PNG', 420, 690, true],
  ['save_image_height_width_svg_1', 'SVG', 420, 690, false],
  ['save_image_height_width_svg_2', 'SVG', 690, 420, false],
  ['transparent_background_svg_1', 'SVG', 400, 400, true],
  ['transparent_background_svg_2', 'SVG', 420, 690, true],
];

// The source suite intentionally performs 33 browser-visible imports/exports.
test.setTimeout(180_000);

function pngDimensions(image) {
  expect(image.subarray(0, 8).toString('hex')).toBe('89504e470d0a1a0a');
  expect(image.subarray(12, 16).toString('ascii')).toBe('IHDR');
  return {
    colorType: image[25],
    height: image.readUInt32BE(20),
    width: image.readUInt32BE(16),
  };
}

function v3000Block(text, name) {
  const lines = String(text).split(/\r?\n/);
  const begin = lines.indexOf(`M  V30 BEGIN ${name}`);
  const end = lines.indexOf(`M  V30 END ${name}`);
  if (begin < 0 || end <= begin) throw new Error(`Missing V3000 ${name} block`);
  return lines.slice(begin + 1, end);
}

/**
 * Compare an SDF's molecular graph independently of V3000 atom numbering or
 * 2-D coordinates. Desktop Squish and Qt/WASM can serialize those differently
 * even when the browser-visible structure and chemistry are the same.
 */
function v3000Graph(text) {
  const atoms = new Map(v3000Block(text, 'ATOM').map((line) => {
    const fields = line.trim().split(/\s+/);
    const charge = fields.find((field) => field.startsWith('CHG=')) || 'CHG=0';
    return [Number(fields[2]), { element: fields[3], charge, edges: [] }];
  }));
  for (const line of v3000Block(text, 'BOND')) {
    const fields = line.trim().split(/\s+/);
    const edge = {
      order: fields[3],
      orientation: fields.find((field) => field.startsWith('CFG=')) || '',
    };
    atoms.get(Number(fields[4])).edges.push({ to: Number(fields[5]), ...edge });
    atoms.get(Number(fields[5])).edges.push({ to: Number(fields[4]), ...edge });
  }
  return atoms;
}

function sameMolecularGraph(actualText, referenceText) {
  const actual = v3000Graph(actualText);
  const reference = v3000Graph(referenceText);
  if (actual.size !== reference.size) return false;
  const label = (atom) => `${atom.element}|${atom.charge}|${atom.edges.length}`;
  const actualIds = [...actual.keys()].sort(
    (left, right) => actual.get(right).edges.length - actual.get(left).edges.length,
  );
  const candidates = new Map(actualIds.map((id) => [
    id,
    [...reference.keys()].filter((candidate) => label(actual.get(id)) === label(reference.get(candidate))),
  ]));
  const mapping = new Map();
  const used = new Set();
  const compatible = (actualId, referenceId) => actual.get(actualId).edges.every((edge) => {
    if (!mapping.has(edge.to)) return true;
    return reference.get(referenceId).edges.some((candidate) =>
      candidate.to === mapping.get(edge.to)
      && candidate.order === edge.order
      && candidate.orientation === edge.orientation,
    );
  });
  const search = (position) => {
    if (position === actualIds.length) return true;
    const actualId = actualIds[position];
    for (const referenceId of candidates.get(actualId)) {
      if (used.has(referenceId) || !compatible(actualId, referenceId)) continue;
      mapping.set(actualId, referenceId);
      used.add(referenceId);
      if (search(position + 1)) return true;
      mapping.delete(actualId);
      used.delete(referenceId);
    }
    return false;
  };
  return search(0);
}

test.describe('tst_export_menu', () => {
  test('main', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    await test.step('import_from_file_10', async () => {
      await sk.import_menu('import_from_file', path.join(FIXTURES_DIR, IMAGE_SOURCE_FIXTURE));
    });
    for (const [name, format, width, height, transparent] of IMAGE_EXPORTS) {
      await test.step(name, async () => {
        const download = await sk.save_image({ filename: name, format, width, height, transparent });
        expect(download.filename).toBe(`${name}.${format.toLowerCase()}`);
        const image = Buffer.from(download.contentBase64, 'base64');
        if (format === 'PNG') {
          expect(pngDimensions(image)).toMatchObject({ width, height, colorType: 6 }); // RGBA
        } else {
          const svg = image.toString('utf8');
          expect(svg).toContain(`<svg width="${width}px" height="${height}px"`);
          expect(svg).toContain(
            transparent ? 'fill-opacity="0"' : 'fill="#ffffff" fill-opacity="1"',
          );
        }
      });
    }

    for (const sourceFragment of SOURCE_FRAGMENT_IDS) {
      await test.step(`import_from_file_${sourceFragment}`, async () => {
        await sk.import_menu(
          'import_from_file',
          path.join(REFERENCES_DIR, `import_from_file_${sourceFragment}_ref.sdf`),
        );
      });
      for (const [extension, format, expectedExtension, expectedText] of TEXT_EXPORTS) {
        await test.step(`save_as_file_${sourceFragment}_${extension}`, async () => {
          const download = await sk.export_menu({
            filename: `fragment_${sourceFragment}_${extension}`,
            format,
          });
          expect(download.filename).toBe(`fragment_${sourceFragment}_${extension}.${expectedExtension}`);
          const contents = Buffer.from(download.contentBase64, 'base64').toString('utf8');
          const referencePath = path.join(
            REFERENCES_DIR,
            `save_as_file_${sourceFragment}_${extension}_ref.${expectedExtension}`,
          );
          // Reference output was produced by the original desktop Squish
          // test. V3000 ordering and coordinates differ in Qt/WASM, so check
          // the chemically meaningful graph rather than byte order.
          const referenceContents = await readFile(referencePath, 'utf8');
          // The original desktop references contain empty InChI files for
          // several enumeration fragments. Preserve that observed behavior
          // instead of incorrectly requiring a payload where none exists.
          if (referenceContents.length === 0) {
            expect(contents).toBe('');
            return;
          }
          expect(contents.length).toBeGreaterThan(0);
          if (expectedText) expect(contents).toContain(expectedText);
          if (format === 'SDF') {
            expect(sameMolecularGraph(contents, referenceContents)).toBe(true);
          } else {
            // PDB's fixed-width coordinates can differ only in a signed zero
            // after the desktop-to-WASM rendering transform.
            const normalizePdbZero = (value) => format === 'PDB'
              ? value.replace(/-0\.000/g, ' 0.000')
              : value;
            expect(normalizePdbZero(contents.trim())).toBe(normalizePdbZero(referenceContents.trim()));
          }
        });
      }
    }
  });
});
