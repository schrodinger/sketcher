#!/usr/bin/env node
/**
 * Generate HELM gallery using sketcher WASM module via Playwright.
 * Usage: node generate_helm_gallery.mjs test/testfiles/helm-gallery/*.csv --wasm-dir path/to/wasm/dir --output gallery.html
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

/**
 * Parse a CSV line respecting quoted fields (RFC 4180-compliant)
 */
function parseCSVLine(line) {
  const result = [];
  let current = '';
  let inQuotes = false;

  for (let i = 0; i < line.length; i++) {
    const char = line[i];
    const nextChar = line[i + 1];

    if (char === '"') {
      if (inQuotes && nextChar === '"') {
        // Escaped quote
        current += '"';
        i++; // Skip next quote
      } else {
        // Toggle quote mode
        inQuotes = !inQuotes;
      }
    } else if (char === ',' && !inQuotes) {
      // End of field
      result.push(current.trim());
      current = '';
    } else {
      current += char;
    }
  }

  // Push last field
  result.push(current.trim());
  return result;
}

// Parse CSV files
const entries = [];
for (const file of positionals) {
  // Handle different line endings (CRLF, LF, CR)
  const content = readFileSync(file, 'utf-8').trim();
  const lines = content.split(/\r?\n/);
  const headers = parseCSVLine(lines[0]);

  for (let i = 1; i < lines.length; i++) {
    if (!lines[i].trim()) continue; // Skip empty lines
    const values = parseCSVLine(lines[i]);
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
const failures = [];
let success = 0, failed = 0;

for (let i = 0; i < entries.length; i++) {
  const entry = entries[i];
  const progressiveNumber = i + 1;

  try {
    const result = await page.evaluate((helmString) => {
      try {
        Module.sketcher_clear();
        Module.sketcher_allow_monomeric();
        Module.sketcher_import_text(helmString);
        const svg = Module.sketcher_export_image(Module.ImageFormat.SVG);
        return { success: true, svg };
      } catch (e) {
        let errorMsg = e.message || e.toString();

        // Try to get the detailed error message from the C++ side
        if (typeof Module.get_last_error === 'function') {
          try {
            const lastError = Module.get_last_error();
            if (lastError && lastError.length > 0) {
              errorMsg = lastError;
            }
          } catch (getErr) {
            // If get_last_error fails, fall back to the exception message
          }
        }

        return { success: false, error: errorMsg };
      }
    }, entry.helm_string);

    if (!result.success) {
      throw new Error(result.error);
    }

    // Decode base64 SVG
    const svgContent = Buffer.from(result.svg, 'base64').toString();

    const searchText = `${entry.helm_string} ${entry.description} ${entry.origin}`.toLowerCase();
    cards.push(`
      <div class="card bg-white rounded-2xl shadow-sm border border-slate-200 overflow-hidden hover:shadow-xl transition-all flex flex-col relative" data-search="${searchText}">
        <div class="absolute top-4 left-4 bg-slate-900 text-white text-xs font-bold px-2 py-1 rounded">#${progressiveNumber}</div>
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
    console.error(`Failed #${progressiveNumber}: ${entry.helm_string} - ${err.message}`);
    failures.push({
      number: progressiveNumber,
      helmString: entry.helm_string,
      description: entry.description,
      origin: entry.origin,
      error: err.message
    });
    failed++;
  }
}

// Cleanup
await browser.close();
server.close();

// Build failures section HTML
const failuresSection = failures.length > 0 ? `
  <div class="mb-12 bg-red-50 border-2 border-red-200 rounded-2xl p-8">
    <h2 class="text-2xl font-bold text-red-900 mb-4">⚠️ Failed Structures</h2>
    <p class="text-red-700 mb-6">The following ${failures.length} structure(s) failed to generate:</p>
    <div class="space-y-4">
      ${failures.map(f => `
        <div class="bg-white border border-red-200 rounded-lg p-4">
          <div class="flex items-start gap-3">
            <span class="bg-red-900 text-white text-xs font-bold px-2 py-1 rounded">#${f.number}</span>
            <div class="flex-grow">
              <p class="text-xs font-bold text-red-600 uppercase mb-1">${f.origin} - ${f.description}</p>
              <p class="text-xs font-mono text-slate-600 break-all bg-slate-50 p-2 rounded border border-slate-200 mb-2">${f.helmString}</p>
              <p class="text-xs text-red-700"><span class="font-bold">Error:</span> ${f.error}</p>
            </div>
          </div>
        </div>
      `).join('')}
    </div>
  </div>
` : '';

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
    <header class="mb-8 border-b border-slate-200 pb-8">
      <h1 class="text-4xl font-extrabold text-slate-900">HELM Gallery</h1>
      <div class="mt-4 grid grid-cols-1 md:grid-cols-3 gap-4">
        <div class="bg-blue-50 border border-blue-200 rounded-lg p-4">
          <p class="text-xs font-bold text-blue-600 uppercase mb-1">Total Input Structures</p>
          <p class="text-3xl font-bold text-blue-900">${entries.length}</p>
        </div>
        <div class="bg-green-50 border border-green-200 rounded-lg p-4">
          <p class="text-xs font-bold text-green-600 uppercase mb-1">Successful</p>
          <p class="text-3xl font-bold text-green-900">${success}</p>
        </div>
        <div class="bg-red-50 border border-red-200 rounded-lg p-4">
          <p class="text-xs font-bold text-red-600 uppercase mb-1">Failed</p>
          <p class="text-3xl font-bold text-red-900">${failed}</p>
        </div>
      </div>
    </header>

    ${failuresSection}

    <div class="mb-6">
      <input type="text" id="search" placeholder="Filter by HELM, description, or origin..."
        class="w-full px-4 py-3 rounded-xl border border-slate-200 shadow-sm focus:ring-2 focus:ring-blue-500 outline-none">
    </div>

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
