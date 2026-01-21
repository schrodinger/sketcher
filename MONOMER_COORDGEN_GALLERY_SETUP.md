# Monomer Coordgen Gallery Setup

This document explains the GitHub Actions workflow that automatically generates a visual gallery of HELM molecules for every build.

## What It Does

**Automatically runs after every successful build** on any branch (including PRs):

1. ✅ Downloads the `helm_to_svg` CLI tool from the build artifacts
2. ✅ Uses test data from `test_data/original_test_data.csv`
3. ✅ Generates SVG images for all 86 HELM molecules
4. ✅ Creates an HTML gallery showing all molecule visualizations
5. ✅ **Deploys to GitHub Pages** for direct browser viewing
6. ✅ Uploads the gallery as a downloadable artifact (backup)
7. ✅ **Posts a comment on PRs** with a browser link to the gallery

## Setup Requirements

**GitHub Pages must be enabled** in repository settings:
1. Go to Settings → Pages
2. Source: **Deploy from a branch**
3. Branch: **gh-pages** / (root)
4. Save

The first workflow run will create the `gh-pages` branch automatically.

## How It Works

The workflow triggers automatically via `workflow_run` after the "Sketcher Build and Test" workflow completes successfully on any branch. No manual intervention needed!

### On Pull Requests

When you open or update a PR:
1. The build workflow runs
2. If successful, the gallery workflow automatically triggers
3. Gallery is deployed to GitHub Pages at `https://<owner>.github.io/<repo>/galleries/<commit>.html`
4. A comment is posted to your PR with a **direct browser link**
5. Click the link to view the gallery instantly - no download needed!

### Manual Trigger

You can also run it manually:

1. Go to Actions → "Monomer Coordgen Gallery"
2. Click "Run workflow"
3. **Leave commit SHA empty** (auto-detects latest build) or enter a specific commit
4. Select your branch
5. Click "Run workflow"



## Viewing the Gallery

### On Pull Requests

**The gallery is automatically generated for ALL PRs and viewable in your browser!**

1. Open or update any PR
2. Wait for builds to complete (~5-10 minutes)
3. A comment will be posted with a **browser link** to the gallery
4. **Click "View Gallery in Browser"** - opens instantly, no download needed
5. _(Alternative: Download from artifacts if needed)_

### On Regular Pushes (non-PR branches)

**Browse all galleries at:** `https://<owner>.github.io/<repo>/`

Or download specific galleries:
1. Go to Actions tab
2. Find the "Monomer Coordgen Gallery" workflow run for your commit
3. Scroll to "Artifacts" section at the bottom
4. Download `helm-gallery-<commit-sha>`
5. Extract and open `helm_gallery.html`

## Gallery Features

The generated HTML gallery includes:

- 🎨 Visual SVG rendering of each HELM molecule
- 🔍 Search/filter by HELM string, description, or origin
- 📱 Responsive grid layout (works on mobile/desktop)
- 🎯 Shows HELM string, description, and origin for each entry
- ⚡ Single-file HTML (no external dependencies except Tailwind CDN)

## Troubleshooting

### "helm_to_svg: command not found"

- The build step failed
- Check the "Build helm_to_svg tool" step logs in the workflow

### Gallery shows no images

- Check the "Generate HELM gallery" step logs
- Look for errors like "Error generating SVG for HELM: ..."
- Invalid HELM strings will be skipped

### Workflow doesn't trigger

- Verify you changed one of the monitored files
- Check the `paths` filter in the workflow file matches your changes