"""Pytest configuration for WASM tests."""

import http.server
import os
import socket
import socketserver
import threading
import time
from pathlib import Path

import pytest

SERVER_PORT = 8000
SERVER_DIR = Path(__file__).parent.parent.parent / "build" / "sketcher_app"


def pytest_configure(config):
    """
    Configure Playwright base URL and screenshot comparison tolerance.
    """
    if not config.option.base_url:
        config.option.base_url = f"http://localhost:{SERVER_PORT}"
    config.option.playwright_max_diff_pixel_ratio = 0.1


def pytest_collection_modifyitems(config, items):
    """
    Skip WASM tests when build artifacts are missing.
    """
    if not (SERVER_DIR / "wasm_shell.html").exists():
        skip_wasm = pytest.mark.skip(reason="WASM build artifacts not found")
        wasm_test_dir = Path(__file__).parent
        for item in items:
            if Path(str(item.fspath)).is_relative_to(wasm_test_dir):
                item.add_marker(skip_wasm)


def _wait_for_server(port: int, timeout: float = 10.0) -> bool:
    """
    Wait for server to become available on the given port.
    """
    start_time = time.time()
    while time.time() - start_time < timeout:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.settimeout(1.0)
            try:
                sock.connect(("localhost", port))
                return True
            except (socket.timeout, ConnectionRefusedError, OSError):
                time.sleep(0.1)
    return False


class _QuietHTTPHandler(http.server.SimpleHTTPRequestHandler):
    """
    HTTP handler that only logs errors.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(SERVER_DIR), **kwargs)

    def log_message(self, format, *args):
        if args[1][0] not in ("2", "3"):  # Log non-2xx/3xx status codes
            super().log_message(format, *args)


@pytest.fixture(scope="session", autouse=True)
def web_server():
    """
    Start HTTP server to serve WASM artifacts.

    Handles pytest-xdist: if multiple workers race to start the server,
    the losers will wait and reuse the winner's server.
    """
    # Reuse existing server (from another xdist worker or manual start)
    if _wait_for_server(SERVER_PORT, timeout=0.1):
        if not os.getenv("CI"):
            print(f"Using existing server on port {SERVER_PORT}")
        yield
        return

    # Try to start the server
    httpd = None
    try:
        httpd = socketserver.TCPServer(("", SERVER_PORT), _QuietHTTPHandler)
        server_thread = threading.Thread(target=httpd.serve_forever,
                                         daemon=True)
        server_thread.start()

        if not _wait_for_server(SERVER_PORT):
            raise RuntimeError(f"Server failed to start on port {SERVER_PORT}")

        print(f"Started web server on http://localhost:{SERVER_PORT}")
    except OSError:
        # Another xdist worker won the race - wait for their server
        if not _wait_for_server(SERVER_PORT):
            raise RuntimeError(
                f"Server on port {SERVER_PORT} never became available")
        print(f"Using server started by another worker on port {SERVER_PORT}")

    yield

    if httpd:
        httpd.shutdown()
        httpd.server_close()


@pytest.fixture(scope="session")
def browser_context_args(browser_context_args):
    """
    Configure browser viewport size.
    """
    return {
        **browser_context_args,
        "viewport": {
            "width": 1280,
            "height": 720
        },
    }
