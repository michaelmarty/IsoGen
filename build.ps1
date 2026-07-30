# Build script for IsoGen - automates pre-release build steps
# Runs up to (but not including) the test step

param(
    [switch]$Help
)

if ($Help) {
    Write-Host @"
IsoGen Build Script
===================

This script automates the build and validation steps for releasing IsoGen.
It runs up to (but not including) the wheel testing step.

Usage:
  .\build.ps1              Run the full build process
  .\build.ps1 -Help        Show this help message

Steps performed:
  1. Check git status
  2. Run tests
  3. Install/upgrade build tools (build, twine)
  4. Clean previous build artifacts (build, dist, wheelhouse)
  5. Build source distribution and wheel
  6. Validate wheel with twine check

The script will stop at any error. Fix the issue and re-run.

Next steps after this script:
  - Test the wheel in a clean virtual environment
  - Upload to TestPyPI or PyPI using twine
"@
    exit 0
}

$ErrorActionPreference = "Stop"

Write-Host "╔════════════════════════════════════════════════════════════╗" -ForegroundColor Cyan
Write-Host "║           IsoGen Build Script                              ║" -ForegroundColor Cyan
Write-Host "╚════════════════════════════════════════════════════════════╝" -ForegroundColor Cyan

# Step 1: Check git status
Write-Host "`n[1/6] Checking git status..." -ForegroundColor Yellow
$gitStatus = & git status --porcelain
if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: Not a git repository or git is not available" -ForegroundColor Red
    exit 1
}

if ($gitStatus) {
    Write-Host "⚠️  Working directory has uncommitted changes:" -ForegroundColor Yellow
    Write-Host $gitStatus
    Write-Host "`nPlease commit your changes before releasing:" -ForegroundColor Yellow
    Write-Host "  git add ." -ForegroundColor Cyan
    Write-Host "  git commit -m 'Description of changes'" -ForegroundColor Cyan
    exit 1
}
Write-Host "✅ Git working directory is clean" -ForegroundColor Green

# Step 2: Run tests
Write-Host "`n[2/6] Running tests..." -ForegroundColor Yellow
try {
    & python -m pytest tests/ -v
    if ($LASTEXITCODE -ne 0) {
        Write-Host "❌ Tests failed" -ForegroundColor Red
        exit 1
    }
    Write-Host "✅ All tests passed" -ForegroundColor Green
} catch {
    Write-Host "❌ Error running tests: $_" -ForegroundColor Red
    exit 1
}

# Step 3: Install/upgrade build tools
Write-Host "`n[3/6] Installing/upgrading build tools..." -ForegroundColor Yellow
try {
    & python -m pip install --upgrade build twine -q
    Write-Host "✅ Build tools ready" -ForegroundColor Green
} catch {
    Write-Host "❌ Error installing build tools: $_" -ForegroundColor Red
    exit 1
}

# Step 4: Clean previous build artifacts
Write-Host "`n[4/6] Cleaning previous build artifacts..." -ForegroundColor Yellow
try {
    $pathsToClean = @('build', 'dist', 'wheelhouse')
    foreach ($path in $pathsToClean) {
        if (Test-Path $path) {
            Remove-Item $path -Recurse -Force -ErrorAction SilentlyContinue
            Write-Host "  Removed: $path" -ForegroundColor Gray
        }
    }
    Write-Host "✅ Build artifacts cleaned" -ForegroundColor Green
} catch {
    Write-Host "❌ Error cleaning artifacts: $_" -ForegroundColor Red
    exit 1
}

# Step 5: Build source distribution and wheel
Write-Host "`n[5/6] Building source distribution and wheel..." -ForegroundColor Yellow
try {
    & python -m build
    if ($LASTEXITCODE -ne 0) {
        Write-Host "❌ Build failed" -ForegroundColor Red
        exit 1
    }
    Write-Host "✅ Build completed successfully" -ForegroundColor Green
} catch {
    Write-Host "❌ Error during build: $_" -ForegroundColor Red
    exit 1
}

# Step 6: Validate with twine
Write-Host "`n[6/6] Validating artifacts with twine..." -ForegroundColor Yellow
try {
    & python -m twine check dist/*
    if ($LASTEXITCODE -ne 0) {
        Write-Host "❌ Twine validation failed" -ForegroundColor Red
        exit 1
    }
    Write-Host "✅ All artifacts validated" -ForegroundColor Green
} catch {
    Write-Host "❌ Error during validation: $_" -ForegroundColor Red
    exit 1
}

# Summary
Write-Host "`n╔════════════════════════════════════════════════════════════╗" -ForegroundColor Green
Write-Host "║           ✅ Build Successful                              ║" -ForegroundColor Green
Write-Host "╚════════════════════════════════════════════════════════════╝" -ForegroundColor Green

Write-Host "`nArtifacts ready in ./dist/" -ForegroundColor Cyan
Write-Host "`nNext steps:" -ForegroundColor Cyan
Write-Host "  1. Test the wheel in a clean virtual environment:" -ForegroundColor Cyan
Write-Host "     python -m venv wheel-test" -ForegroundColor Gray
Write-Host "     wheel-test\Scripts\pip install dist\*.whl" -ForegroundColor Gray
Write-Host "     wheel-test\Scripts\python -c 'import isogen; print(isogen.isodist(1000, isolen=8))'" -ForegroundColor Gray
Write-Host "" -ForegroundColor Cyan
Write-Host "  2. Upload to TestPyPI (first time) or PyPI:" -ForegroundColor Cyan
Write-Host "     python -m twine upload --repository testpypi dist/*  # TestPyPI" -ForegroundColor Gray
Write-Host "     python -m twine upload dist/*  # PyPI" -ForegroundColor Gray
