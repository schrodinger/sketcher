import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[2] / "wasm/update_cache_bust.py"
SPEC = importlib.util.spec_from_file_location("update_cache_bust", SCRIPT)
CACHE_BUST = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CACHE_BUST)


class CacheBustTest(unittest.TestCase):

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        for name, content in {
                "Sketcher.wasm": "main wasm",
                "Sketcher.js": "'Sketcher.wasm'; 'custom_monomers.worker.js'",
                "MonomerLoader.wasm": "worker wasm",
                "MonomerLoader.js": "'MonomerLoader.wasm'",
                "custom_monomers.worker.js": "importScripts('MonomerLoader.js')",
                "qtloader.js": "qt loader",
                "wasm_shell.html": "'Sketcher.js'; 'qtloader.js'",
        }.items():
            (self.root / name).write_text(content)

    def update(self):
        with patch.object(sys, "argv", [
                str(SCRIPT), "--wasm-dir",
                str(self.root), "--app-target", "Sketcher"
        ]):
            CACHE_BUST.main()

    def snapshot(self):
        return {file.name: file.read_bytes() for file in self.root.iterdir()}

    def test_repeated_runs_are_stable(self):
        self.update()
        previous = self.snapshot()
        self.update()
        self.assertEqual(self.snapshot(), previous)
        for name in [
                "Sketcher.js", "MonomerLoader.js", "custom_monomers.worker.js",
                "wasm_shell.html"
        ]:
            self.assertIn(b"?cache_bust=", previous[name])

    def test_worker_change_invalidates_main_module_and_shell(self):
        self.update()
        previous = self.snapshot()
        (self.root / "MonomerLoader.wasm").write_text("changed worker wasm")
        self.update()
        updated = self.snapshot()
        for name in [
                "MonomerLoader.js", "custom_monomers.worker.js", "Sketcher.js",
                "wasm_shell.html"
        ]:
            self.assertNotEqual(updated[name], previous[name])
        self.assertEqual(updated["Sketcher.wasm"], previous["Sketcher.wasm"])


if __name__ == "__main__":
    unittest.main()
