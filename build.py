#!/usr/bin/env python
"""
Build script for IsoGen - automates pre-release build steps.
Runs up to (but not including) the wheel testing step.

Usage:
    python build.py              Run the full build process
    python build.py --help       Show this help message
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

# Handle Windows encoding issues
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')


def run_command(cmd, description, check=True):
    """Run a shell command and handle errors."""
    print(f"\n[*] {description}...")
    try:
        result = subprocess.run(cmd, shell=True, check=check)
        return result.returncode == 0
    except Exception as e:
        print(f"[!] Error: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="IsoGen build script - automates pre-release build steps"
    )
    parser.add_argument(
        "--skip-tests",
        action="store_true",
        help="Skip running tests (useful for debugging)"
    )
    args = parser.parse_args()

    print("\n" + "=" * 60)
    print("  IsoGen Build Script")
    print("=" * 60)

    # Step 1: Check git status
    print("\n[1/6] Checking git status...")
    try:
        result = subprocess.run(
            "git status --porcelain",
            shell=True,
            capture_output=True,
            text=True,
            check=True
        )
        if result.stdout.strip():
            print("[!] Working directory has uncommitted changes:")
            print(result.stdout)
            print("\nPlease commit your changes before releasing:")
            print("  git add .")
            print("  git commit -m 'Description of changes'")
            sys.exit(1)
        print("[OK] Git working directory is clean")
    except subprocess.CalledProcessError as e:
        print(f"[!] Error checking git status: {e}")
        sys.exit(1)

    # Step 2: Run tests
    if not args.skip_tests:
        print("\n[2/6] Running tests...")
        if not run_command("python -m pytest tests/ -v", "Tests", check=False):
            print("[!] Tests failed")
            sys.exit(1)
        print("[OK] All tests passed")
    else:
        print("\n[2/6] Skipping tests (--skip-tests)")

    # Step 3: Install/upgrade build tools
    step = 3 if not args.skip_tests else 2
    print(f"\n[{step}/6] Installing/upgrading build tools...")
    if not run_command(
        "python -m pip install --upgrade build twine",
        "Installing build tools",
        check=False
    ):
        print("[!] Error installing build tools")
        sys.exit(1)
    print("[OK] Build tools ready")

    # Step 4: Clean previous build artifacts
    step = 4 if not args.skip_tests else 3
    print(f"\n[{step}/6] Cleaning previous build artifacts...")
    paths_to_clean = ["build", "dist", "wheelhouse"]
    for path in paths_to_clean:
        if Path(path).exists():
            shutil.rmtree(path, ignore_errors=True)
            print(f"  - Removed: {path}")
    print("[OK] Build artifacts cleaned")

    # Step 5: Build source distribution and wheel
    step = 5 if not args.skip_tests else 4
    print(f"\n[{step}/6] Building source distribution and wheel...")
    if not run_command("python -m build", "Building", check=False):
        print("[!] Build failed")
        sys.exit(1)
    print("[OK] Build completed successfully")

    # Step 6: Validate with twine
    step = 6 if not args.skip_tests else 5
    print(f"\n[{step}/6] Validating artifacts with twine...")
    if not run_command("python -m twine check dist/*", "Twine validation", check=False):
        print("[!] Twine validation failed")
        sys.exit(1)
    print("[OK] All artifacts validated")

    # Success summary
    print("\n" + "=" * 60)
    print("  [OK] Build Successful")
    print("=" * 60)
    print("\nArtifacts ready in ./dist/")
    print("\nNext steps:")
    print("  1. Test the wheel in a clean virtual environment:")
    print("     python -m venv wheel-test")
    print("     ./wheel-test/bin/python -m pip install dist/isogen-*.whl")
    print("     ./wheel-test/bin/python -c 'import isogen; print(isogen.isodist(1000, isolen=8))'")
    print("")
    print("  2. Upload to TestPyPI (first time) or PyPI:")
    print("     python -m twine upload --repository testpypi dist/*  # TestPyPI")
    print("     python -m twine upload dist/*  # PyPI")
    print("")


if __name__ == "__main__":
    main()
