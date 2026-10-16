"""
Update cache-busting hashes for the WebAssembly application assets.

Browsers may keep JavaScript and WebAssembly files in cache after a deployment.
This script appends each referenced asset's MD5 digest as a ``cache_bust`` query
parameter, so a changed file is fetched again. It updates the dependency chain
from the monomer loader WASM through its worker and the main Sketcher module,
then updates the shell's references to the main module and Qt loader. The
resulting hashes are deterministic, so running the script repeatedly without
changing assets does not keep modifying the files.
"""

import argparse
import hashlib
import re
from pathlib import Path


def md5(f):
    return hashlib.md5(f.read_bytes()).hexdigest()


def update_reference(path, name, digest):
    content = path.read_text()
    content = re.sub(
        re.escape(name) + r"(?:\?cache_bust=[0-9a-f]+)?",
        f"{name}?cache_bust={digest}", content)
    path.write_text(content)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wasm-dir", required=True)
    parser.add_argument("--app-target", required=True)
    args = parser.parse_args()

    wasm_dir = Path(args.wasm_dir)

    # Propagate worker changes through its script, the main module, and shell.
    worker_module = wasm_dir / "MonomerLoader.js"
    update_reference(worker_module, "MonomerLoader.wasm",
                     md5(wasm_dir / "MonomerLoader.wasm"))
    worker = wasm_dir / "custom_monomers.worker.js"
    update_reference(worker, "MonomerLoader.js", md5(worker_module))

    # Hash the .wasm file and update references in the .js file
    wasm_hash = md5(wasm_dir / f"{args.app_target}.wasm")
    sketcher_js = wasm_dir / f"{args.app_target}.js"
    update_reference(sketcher_js, f"{args.app_target}.wasm", wasm_hash)
    update_reference(sketcher_js, "custom_monomers.worker.js", md5(worker))

    # Hash the modified .js file and qtloader.js
    js_hash = md5(sketcher_js)
    qtloader = wasm_dir / "qtloader.js"
    qt_hash = md5(qtloader)

    # Update wasm_shell.html to cache-bust both JS files
    shell_html = wasm_dir / "wasm_shell.html"
    if shell_html.exists():
        update_reference(shell_html, "qtloader.js", qt_hash)
        update_reference(shell_html, f"{args.app_target}.js", js_hash)


if __name__ == "__main__":
    main()
