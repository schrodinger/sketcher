#!/usr/bin/env python3
"""
Generate a visual HELM gallery from CSV files using sketcher's helm_to_svg tool.
This is a modified version that works without $SCHRODINGER/run and bbchem_endpoints.
"""

import argparse
import csv
import subprocess
import sys
from pathlib import Path
from typing import List, NamedTuple

DEFAULT_OUTPUT_HTML = "helm_gallery.html"



class HelmEntry(NamedTuple):
    helm_string: str
    smiles_string: str
    description: str
    origin: str

def parse_helm_csv(file_paths: List[str]) -> List[HelmEntry]:
    """Parse CSV files containing HELM data."""
    helm_entries = []
    for file_path in file_paths:
        with open(file_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                helm_entry = HelmEntry(
                    helm_string=row['helm_string'],
                    smiles_string=row.get('smiles_string', ''),
                    description=row.get('description', ''),
                    origin=row.get('origin', '')
                )
                helm_entries.append(helm_entry)
    return helm_entries


def generate_svg(helm_string: str, helm_to_svg_path: str) -> str:
    """Generate SVG from HELM string using the helm_to_svg tool."""
    try:
        result = subprocess.run(
            [helm_to_svg_path, helm_string],
            capture_output=True,
            text=True,
            timeout=10
        )
        if result.returncode == 0:
            # Clean up SVG for inline embedding
            svg_content = result.stdout.replace('\n', '').replace('\r', '')
            return svg_content
        else:
            print(f"Error generating SVG for HELM: {helm_string}")
            print(f"stderr: {result.stderr}")
            return None
    except subprocess.TimeoutExpired:
        print(f"Timeout generating SVG for HELM: {helm_string}")
        return None
    except Exception as e:
        print(f"Exception generating SVG for HELM: {helm_string}. Error: {e}")
        return None


def create_webpage(entries: List[HelmEntry], output_file: str, helm_to_svg_path: str):
    """Generate a single-file HTML gallery for the provided HelmEntry objects."""
    if not entries:
        print("No entries found to generate gallery.")
        return

    # Header of the HTML using Tailwind CSS
    html_start = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HELM Visualization Gallery</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <style>
        .svg-container svg {{ max-width: 100%; height: auto; }}
    </style>
</head>
<body class="bg-slate-50 min-h-screen p-6 md:p-12 font-sans">
    <div class="max-w-7xl mx-auto">
        <header class="mb-12 border-b border-slate-200 pb-8 flex flex-col md:flex-row md:items-end justify-between gap-6">
            <div>
                <h1 class="text-4xl font-extrabold text-slate-900 tracking-tight">HELM Gallery</h1>
                <p class="text-slate-500 mt-2 text-lg">Visualizing {len(entries)} monomeric entries.</p>
            </div>
            <div class="w-full md:w-80">
                <label class="block text-xs font-bold uppercase text-slate-400 mb-2 ml-1">Search Entries</label>
                <input type="text" id="search" placeholder="Filter by HELM, description, or origin..."
                       class="w-full px-4 py-3 rounded-xl border border-slate-200 shadow-sm focus:ring-2 focus:ring-blue-500 outline-none transition-all">
            </div>
        </header>

        <div id="grid" class="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
    '''

    html_cards = []
    successful = 0
    failed = 0

    for entry in entries:
        svg_content = generate_svg(entry.helm_string, helm_to_svg_path)
        if svg_content is None:
            failed += 1
            continue

        successful += 1

        # We store searchable text in a data attribute
        search_blob = f"{entry.helm_string} {entry.description} {entry.origin}".lower()

        card = f'''
            <div class="card bg-white rounded-2xl shadow-sm border border-slate-200 overflow-hidden hover:shadow-xl transition-all flex flex-col"
                 data-search="{search_blob}">
                <!-- SVG Preview -->
                <div class="p-8 flex-grow flex items-center justify-center bg-white min-h-[250px] svg-container border-b border-slate-50">
                    {svg_content}
                </div>

                <!-- Metadata Area -->
                <div class="p-5 flex flex-col gap-3">
                    <div class="flex items-start justify-between gap-2">
                        <span class="px-2 py-1 bg-blue-50 text-blue-600 text-[10px] font-bold uppercase rounded tracking-wider whitespace-nowrap">
                            {entry.origin}
                        </span>
                    </div>

                    <div>
                        <h3 class="text-sm font-semibold text-slate-800 line-clamp-2" title="{entry.description}">
                            {entry.description}
                        </h3>
                    </div>

                    <div class="mt-2 pt-3 border-t border-slate-100">
                        <p class="text-[10px] font-mono text-slate-400 uppercase font-bold mb-1">HELM String</p>
                        <p class="text-[11px] font-mono text-slate-600 break-all bg-slate-50 p-2 rounded border border-slate-100">
                            {entry.helm_string}
                        </p>
                    </div>
                </div>
            </div>
        '''
        html_cards.append(card)

    html_end = '''
        </div>

        <footer class="mt-20 py-8 border-t border-slate-200 text-center text-slate-400 text-sm">
            HELM Gallery Generated via Schrodinger Sketcher
        </footer>
    </div>

    <script>
        const searchInput = document.getElementById('search');
        const cards = document.querySelectorAll('.card');

        searchInput.addEventListener('input', (e) => {
            const query = e.target.value.toLowerCase();
            cards.forEach(card => {
                const text = card.getAttribute('data-search');
                if (text.includes(query)) {
                    card.style.display = 'flex';
                } else {
                    card.style.display = 'none';
                }
            });
        });
    </script>
</body>
</html>
    '''

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_start + "".join(html_cards) + html_end)

    print(f"\nGallery generated: {successful} successful, {failed} failed")
    print(f"Written to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a visual HELM gallery from CSV files.")
    parser.add_argument(
        "csv_files",
        nargs="+",
        help="Path(s) to the HELM CSV files."
    )
    parser.add_argument(
        "-o", "--output",
        default=DEFAULT_OUTPUT_HTML,
        help=f"Output HTML file path (default: {DEFAULT_OUTPUT_HTML})"
    )
    parser.add_argument(
        "--helm-to-svg",
        required=True,
        help="Path to the helm_to_svg executable"
    )

    args = parser.parse_args()

    # Validate helm_to_svg exists
    helm_to_svg_path = Path(args.helm_to_svg)
    if not helm_to_svg_path.exists():
        print(f"Error: helm_to_svg not found at {helm_to_svg_path}", file=sys.stderr)
        sys.exit(1)

    all_entries = parse_helm_csv(args.csv_files)
    create_webpage(all_entries, args.output, str(helm_to_svg_path))
