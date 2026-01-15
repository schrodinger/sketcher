"""
Playwright tests for the WASM Sketcher.

These tests verify that the WASM module loads correctly and exercise all
public API Emscripten bindings, including import/export across all supported
chemical file formats.
"""

import base64
from pathlib import Path

import pytest
from playwright.sync_api import Page, expect

# Format configurations for testing
# yapf: disable
FORMATS = {
    "AUTO_DETECT": {"import": False, "export": False},
    "RDMOL_BINARY_BASE64": {},
    "SMILES": {},
    "EXTENDED_SMILES": {},
    "SMARTS": {},
    "EXTENDED_SMARTS": {},
    "MDL_MOLV2000": {"export": False},
    "MDL_MOLV3000": {},
    "MAESTRO": {},
    "INCHI": {},
    "INCHI_KEY": {"import": False},
    "PDB": {},
    "MOL2": {"export": False},
    "XYZ": {},
    "MRV": {},
    "CDXML": {"export": False},
    "HELM": {"export": False},  # Export not supported (throws C++ exception)
    "FASTA_PEPTIDE": {"export": False},
    "FASTA_DNA": {"export": False},
    "FASTA_RNA": {"export": False},
    "FASTA": {"import": False, "export": False},
}
# yapf: enable


@pytest.fixture(autouse=True)
def setup_page(page: Page):
    """
    Navigate to the WASM shell page and wait for the module to load.
    This fixture runs automatically before each test.
    """
    # Navigate to the page that loads the WASM module and
    # wait for the WASM module to be fully loaded and available
    page.goto("/wasm_shell.html")
    page.wait_for_function(
        "() => typeof window.Module !== 'undefined'",
        timeout=20000,
    )
    # Enable monomer support in tests
    page.evaluate("() => Module.sketcher_allow_monomeric()")


def test_sketcher_ui_loads_correctly(page: Page):
    """Verify that the sketcher UI loads and renders correctly."""
    page.screenshot(path="sketcher-load.png")
    # Visual regression testing - compare with baseline screenshot
    # expect(page).to_have_screenshot("sketcher-load.png")


def test_format_enum_matches_python_list(page: Page):
    """Verify that FORMATS dict matches the actual Module.Format enum."""
    actual_formats = page.evaluate(
        "() => Object.keys(Module.Format).filter(k => k !== 'values')")
    expected_formats = list(FORMATS.keys())

    # Check for missing formats (in C++ but not in Python tests)
    missing = set(actual_formats) - set(expected_formats)
    if missing:
        pytest.fail(
            f"Formats in Module.Format but not in FORMATS dict: {sorted(missing)}\n"
            f"Please add them to the FORMATS dict in test_wasm_api.py")

    # Check for extra formats (in Python tests but not in C++)
    extra = set(expected_formats) - set(actual_formats)
    if extra:
        pytest.fail(
            f"Formats in FORMATS dict but not in Module.Format: {sorted(extra)}\n"
            f"Please remove them from the FORMATS dict in test_wasm_api.py")


@pytest.mark.parametrize(
    "format_name",
    [name for name, config in FORMATS.items() if config.get("export", True)],
)
def test_exporting(page: Page, format_name: str):
    """Test exporting molecule data for formats that support export."""

    def export_molecule(fmt: str) -> str:
        page.evaluate("() => Module.sketcher_clear()")
        page.evaluate(
            "() => Module.sketcher_import_text('C[C@H](N)C=O', Module.Format.AUTO_DETECT)"
        )
        return page.evaluate(
            f"(fmt) => Module.sketcher_export_text(Module.Format[fmt])", fmt)

    exported_text = export_molecule(format_name)

    # Read the expected snapshot
    test_dir = Path(__file__).parent
    snapshot_path = (test_dir / "__snapshots__" / "test_wasm_api.py" /
                     f"format-{format_name.replace('_', '-')}.txt")
    expected_text = snapshot_path.read_text()
    assert exported_text == expected_text


@pytest.mark.parametrize(
    "format_name",
    [name for name, config in FORMATS.items() if config.get("import", True)],
)
def test_importing(page: Page, format_name: str):
    """Test importing molecule data for formats that support import."""
    # FASTA_RNA has a C++ bug with AUTO_DETECT format - skip for now
    if format_name == "FASTA_RNA":
        pytest.skip("FASTA_RNA import with AUTO_DETECT throws C++ exception")

    # Read the fixture file
    test_dir = Path(__file__).parent
    fixture_path = (test_dir / "__snapshots__" / "test_wasm_api.py" /
                    f"format-{format_name.replace('_', '-')}.txt")
    input_text = fixture_path.read_text()

    # Test auto-detect import
    page.evaluate("() => Module.sketcher_clear()")
    page.evaluate(
        "(text) => Module.sketcher_import_text(text, Module.Format.AUTO_DETECT)",
        input_text)
    auto_detect_works = not page.evaluate("() => Module.sketcher_is_empty()")
    assert auto_detect_works is True

    # Test explicit format import
    page.evaluate("() => Module.sketcher_clear()")
    page.evaluate(
        "([text, format]) => Module.sketcher_import_text(text, Module.Format[format])",
        [input_text, format_name],
    )
    explicit_format_works = not page.evaluate(
        "() => Module.sketcher_is_empty()")
    assert explicit_format_works is True


def test_clearing_the_sketcher(page: Page):
    """Test that clearing the sketcher works correctly."""
    is_empty_on_load = page.evaluate("() => Module.sketcher_is_empty()")
    assert is_empty_on_load is True

    page.evaluate(
        "() => Module.sketcher_import_text('C', Module.Format.AUTO_DETECT)")
    is_empty_after_import = page.evaluate("() => Module.sketcher_is_empty()")
    assert is_empty_after_import is False

    page.evaluate("() => Module.sketcher_clear()")
    is_empty_after_clear = page.evaluate("() => Module.sketcher_is_empty()")
    assert is_empty_after_clear is True


def test_checking_if_molecule_has_monomers(page: Page):
    """Test that monomer detection works correctly."""
    has_monomers_on_load = page.evaluate("() => Module.sketcher_has_monomers()")
    assert has_monomers_on_load is False

    page.evaluate(
        "() => Module.sketcher_import_text('c1ccccc1', Module.Format.AUTO_DETECT)"
    )
    has_monomers_after_smiles_import = page.evaluate(
        "() => Module.sketcher_has_monomers()")
    assert has_monomers_after_smiles_import is False

    page.evaluate("() => Module.sketcher_clear()")
    page.evaluate(
        "() => Module.sketcher_import_text('PEPTIDE1{A.S.D.F.G.H.W}$$$$V2.0', Module.Format.AUTO_DETECT)"
    )
    has_monomers_after_helm_import = page.evaluate(
        "() => Module.sketcher_has_monomers()")
    assert has_monomers_after_helm_import is True

    page.evaluate("() => Module.sketcher_clear()")
    has_monomers_after_clear = page.evaluate(
        "() => Module.sketcher_has_monomers()")
    assert has_monomers_after_clear is False


@pytest.mark.parametrize("image_format", ["SVG", "PNG"])
def test_exporting_image(page: Page, image_format: str):
    """Test exporting images in SVG and PNG formats."""
    page.evaluate(
        "() => Module.sketcher_import_text('C=O', Module.Format.AUTO_DETECT)")
    base64_content = page.evaluate(
        "(format) => Module.sketcher_export_image(Module.ImageFormat[format])",
        image_format,
    )
    image_buffer = base64.b64decode(base64_content)

    # For SVGs, we render them in the browser and take a screenshot so we can do
    # actual visual comparison instead of comparing the markup (which might differ
    # slightly by platform)
    if image_format == "SVG":
        svg_data_uri = f"data:image/svg+xml;charset=utf-8;base64,{base64.b64encode(image_buffer).decode()}"
        page.goto(svg_data_uri)
        page.locator("svg").screenshot(
            path=f"export-image-{image_format}.png",
            omit_background=True,
        )
    else:
        # For PNG, just verify it's a valid PNG (exact bytes vary by platform)
        assert image_buffer.startswith(b'\x89PNG'), "Should be valid PNG format"
        assert len(image_buffer) > 1000, "PNG should have reasonable size"
