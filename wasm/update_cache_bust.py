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
    wasm_hash = md5(wasm_dir / f"{args.app_target}.wasm")
    js_hash = md5(wasm_dir / f"{args.app_target}.js")

    qtloader = wasm_dir / "qtloader.js"
    content = qtloader.read_text()
    content = content.replace(".wasm", f".wasm?cache_bust={wasm_hash}")
    content = content.replace('".js"', f'".js?cache_bust={js_hash}"')
    qtloader.write_text(content)

    shell_html = wasm_dir / "wasm_shell.html"
    if shell_html.exists():
        qt_hash = md5(qtloader)
        content = shell_html.read_text()
        content = content.replace("qtloader.js",
                                  f"qtloader.js?cache_bust={qt_hash}")
        shell_html.write_text(content)


if __name__ == "__main__":
    main()
