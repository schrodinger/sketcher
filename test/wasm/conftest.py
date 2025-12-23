"""
Pytest configuration for WASM tests.

This file configures the test environment, including starting a web server
to serve the WASM build artifacts.
"""

import http.server
import os
import socket
import socketserver
import threading
import time
from pathlib import Path

import pytest

# Web server configuration
SERVER_PORT = 8000
SERVER_DIR = Path(__file__).parent.parent.parent / "build" / "sketcher_app"


def pytest_configure(config):
    """
    Configure pytest settings and Playwright defaults.
    """
    # Set base URL for playwright tests
    if not config.option.base_url:
        config.option.base_url = f"http://localhost:{SERVER_PORT}"

    # Set default max diff pixel ratio for screenshot comparisons
    # This allows up to 10% pixel difference in snapshot tests
    config.option.playwright_max_diff_pixel_ratio = 0.1


def pytest_collection_modifyitems(config, items):
    """
    Skip WASM tests if build artifacts are not present.
    """
    # Check if WASM build artifacts exist
    wasm_file = SERVER_DIR / "wasm_shell.html"
    if not wasm_file.exists():
        skip_wasm = pytest.mark.skip(reason="WASM build artifacts not found")
        for item in items:
            # Skip all tests in the wasm directory
            if "test/wasm" in str(item.fspath):
                item.add_marker(skip_wasm)


def is_server_running(port: int, timeout: float = 1.0) -> bool:
    """Check if a server is already running on the given port."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.settimeout(timeout)
        try:
            sock.connect(("localhost", port))
            return True
        except (socket.timeout, ConnectionRefusedError):
            return False


@pytest.fixture(scope="session", autouse=True)
def web_server():
    """
    Start a web server to serve the WASM build artifacts.
    This fixture runs once per test session.
    """
    # Check if server is already running (e.g., in development)
    if is_server_running(SERVER_PORT) and not os.getenv("CI"):
        print(f"Using existing server on port {SERVER_PORT}")
        yield
        return

    # Start the server
    class QuietHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
        """HTTP request handler with minimal logging."""

        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=str(SERVER_DIR), **kwargs)

        def log_message(self, format, *args):
            # Only log errors
            if args[1][0] not in ("2", "3"):
                super().log_message(format, *args)

    httpd = socketserver.TCPServer(("", SERVER_PORT), QuietHTTPRequestHandler)

    # Run server in a separate thread
    server_thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    server_thread.start()

    # Wait for server to be ready
    max_wait = 10  # seconds
    start_time = time.time()
    while not is_server_running(SERVER_PORT):
        if time.time() - start_time > max_wait:
            raise RuntimeError(f"Server failed to start on port {SERVER_PORT}")
        time.sleep(0.1)

    print(f"Started web server on http://localhost:{SERVER_PORT}")

    yield

    # Shutdown server
    httpd.shutdown()
    httpd.server_close()


@pytest.fixture(scope="session")
def browser_context_args(browser_context_args):
    """
    Configure browser context arguments.
    """
    return {
        **browser_context_args,
        "viewport": {
            "width": 1280,
            "height": 720
        },
    }
