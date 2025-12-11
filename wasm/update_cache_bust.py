import argparse
import hashlib
from pathlib import Path


def md5(f):
    return hashlib.md5(f.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wasm-dir", required=True)
    parser.add_argument("--app-target", required=True)
    args = parser.parse_args()

    wasm_dir = Path(args.wasm_dir)

    # Hash the .wasm file and update references in the .js file
    wasm_hash = md5(wasm_dir / f"{args.app_target}.wasm")
    sketcher_js = wasm_dir / f"{args.app_target}.js"
    content = sketcher_js.read_text()
    content = content.replace(f"{args.app_target}.wasm",
                              f"{args.app_target}.wasm?cache_bust={wasm_hash}")
    sketcher_js.write_text(content)

    # Hash the modified .js file and qtloader.js
    js_hash = md5(sketcher_js)
    qtloader = wasm_dir / "qtloader.js"
    qt_hash = md5(qtloader)

    # Update wasm_shell.html to cache-bust both JS files
    shell_html = wasm_dir / "wasm_shell.html"
    if shell_html.exists():
        content = shell_html.read_text()
        content = content.replace("qtloader.js",
                                  f"qtloader.js?cache_bust={qt_hash}")
        content = content.replace(f"{args.app_target}.js",
                                  f"{args.app_target}.js?cache_bust={js_hash}")
        shell_html.write_text(content)


if __name__ == "__main__":
    main()
