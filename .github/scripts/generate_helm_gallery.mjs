#!/usr/bin/env node
/**
 * Generate HELM gallery using sketcher WASM module via Playwright.
 * Usage: node generate_helm_gallery.mjs test_data/*.csv --wasm-dir path/to/wasm/dir --output gallery.html
 */

import { readFileSync, writeFileSync } from 'fs';
import { resolve, dirname } from 'path';
import { parseArgs } from 'util';
import { chromium } from 'playwright';
import { fileURLToPath } from 'url';
import { createServer } from 'http';
import { createReadStream, statSync } from 'fs';
import { extname } from 'path';

const __dirname = dirname(fileURLToPath(import.meta.url));

// Parse CLI arguments
const { values, positionals } = parseArgs({
  options: {
    'wasm-dir': { type: 'string', required: true },
    'output': { type: 'string', default: 'helm_gallery.html' }
  },
  allowPositionals: true
});

if (positionals.length === 0) {
  console.error('Usage: generate_helm_gallery.mjs <csv_files...> --wasm-dir <path> [--output <file>]');
  process.exit(1);
}

// Parse CSV files
const entries = [];
for (const file of positionals) {
  const lines = readFileSync(file, 'utf-8').trim().split('\n');
  const headers = lines[0].split(',').map(h => h.trim());

  for (let i = 1; i < lines.length; i++) {
    const values = lines[i].split(',').map(v => v.trim());
    const entry = {};
    headers.forEach((h, idx) => entry[h] = values[idx] || '');
    entries.push(entry);
  }
}

console.log(`Processing ${entries.length} entries...`);

// Start HTTP server to serve WASM files
const wasmDir = resolve(values['wasm-dir']);
const server = createServer((req, res) => {
  const filePath = resolve(wasmDir, req.url.slice(1).split('?')[0]); // Strip query params

  // Security: ensure we're serving from wasmDir
  if (!filePath.startsWith(wasmDir)) {
    res.writeHead(403);
    res.end('Forbidden');
    return;
  }

  try {
    const stat = statSync(filePath);
    if (!stat.isFile()) {
      res.writeHead(404);
      res.end('Not found');
      return;
    }

    // Set content type
    const ext = extname(filePath);
    const contentTypes = {
      '.html': 'text/html',
      '.js': 'application/javascript',
      '.wasm': 'application/wasm',
      '.png': 'image/png',
      '.gif': 'image/gif',
      '.svg': 'image/svg+xml'
    };
    res.writeHead(200, { 'Content-Type': contentTypes[ext] || 'application/octet-stream' });
    createReadStream(filePath).pipe(res);
  } catch (err) {
    res.writeHead(404);
    res.end('Not found');
  }
});

await new Promise((resolve) => {
  server.listen(0, 'localhost', () => {
    console.log(`HTTP server started on http://localhost:${server.address().port}`);
    resolve();
  });
});

const baseUrl = `http://localhost:${server.address().port}`;

// Launch browser and load WASM
const browser = await chromium.launch();
const page = await browser.newPage();

// Navigate to the WASM shell
await page.goto(`${baseUrl}/wasm_shell.html`);
await page.waitForFunction(() => typeof window.Module !== 'undefined', { timeout: 30000 });

console.log('WASM module loaded in browser.');

// Generate SVGs
const cards = [];
let success = 0, failed = 0;

for (const entry of entries) {
  try {
    const svg = await page.evaluate((helmString) => {
      Module.sketcher_clear();
      Module.sketcher_allow_monomeric();
      Module.sketcher_import_text(helmString);
      return Module.sketcher_export_image(Module.ImageFormat.SVG);
    }, entry.helm_string);

    // Decode base64 SVG
    const svgContent = Buffer.from(svg, 'base64').toString();

    const searchText = `${entry.helm_string} ${entry.description} ${entry.origin}`.toLowerCase();
    cards.push(`
      <div class="card bg-white rounded-2xl shadow-sm border border-slate-200 overflow-hidden hover:shadow-xl transition-all flex flex-col" data-search="${searchText}">
        <div class="p-8 flex-grow flex items-center justify-center bg-white min-h-[250px]">${svgContent}</div>
        <div class="p-5 flex flex-col gap-3">
          <span class="px-2 py-1 bg-blue-50 text-blue-600 text-[10px] font-bold uppercase rounded tracking-wider w-fit">${entry.origin}</span>
          <h3 class="text-sm font-semibold text-slate-800 line-clamp-2">${entry.description}</h3>
          <div class="mt-2 pt-3 border-t border-slate-100">
            <p class="text-[10px] font-mono text-slate-400 uppercase font-bold mb-1">HELM String</p>
            <p class="text-[11px] font-mono text-slate-600 break-all bg-slate-50 p-2 rounded border border-slate-100">${entry.helm_string}</p>
          </div>
        </div>
      </div>`);
    success++;
  } catch (err) {
    console.error(`Failed: ${entry.helm_string} - ${err.message}`);
    failed++;
  }
}

// Cleanup
await browser.close();
server.close();

// Write HTML
const html = `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>HELM Gallery</title>
  <script src="https://cdn.tailwindcss.com"></script>
  <style>.svg-container svg { max-width: 100%; height: auto; }</style>
</head>
<body class="bg-slate-50 min-h-screen p-6 md:p-12">
  <div class="max-w-7xl mx-auto">
    <header class="mb-12 border-b border-slate-200 pb-8 flex flex-col md:flex-row md:items-end justify-between gap-6">
      <div>
        <h1 class="text-4xl font-extrabold text-slate-900">HELM Gallery</h1>
        <p class="text-slate-500 mt-2 text-lg">Visualizing ${entries.length} entries.</p>
      </div>
      <input type="text" id="search" placeholder="Filter by HELM, description, or origin..."
        class="w-full md:w-80 px-4 py-3 rounded-xl border border-slate-200 shadow-sm focus:ring-2 focus:ring-blue-500 outline-none">
    </header>
    <div id="grid" class="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
      ${cards.join('')}
    </div>
    <footer class="mt-20 py-8 border-t border-slate-200 text-center text-slate-400 text-sm">
      HELM Gallery via Schrodinger Sketcher WASM
    </footer>
  </div>
  <script>
    const search = document.getElementById('search');
    const cards = document.querySelectorAll('.card');
    search.addEventListener('input', (e) => {
      const q = e.target.value.toLowerCase();
      cards.forEach(card => {
        card.style.display = card.getAttribute('data-search').includes(q) ? 'flex' : 'none';
      });
    });
  </script>
</body>
</html>`;

writeFileSync(values.output, html);
console.log(`\nGallery: ${success} successful, ${failed} failed → ${values.output}`);
