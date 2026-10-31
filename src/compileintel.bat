@echo off
rem Builds IsoGen with the Intel oneAPI C compiler (icx) and installs into ..\isogen\bin
setlocal

set "INTEL_ROOT=C:/PROGRA~2/Intel/oneAPI/compiler/2026.1"
set "VCVARS=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
set "SRC_DIR=%~dp0"
set "BUILD_DIR=%SRC_DIR%cmake-build-release-intel"

if not exist "%VCVARS%" (
    echo Visual Studio vcvars64.bat not found: %VCVARS%
    exit /b 1
)
if not exist "%INTEL_ROOT%/bin/icx.exe" (
    echo Intel icx.exe not found: %INTEL_ROOT%/bin/icx.exe
    exit /b 1
)

call "%VCVARS%" >nul
if errorlevel 1 exit /b 1

cmake -S "%SRC_DIR%." -B "%BUILD_DIR%" -G Ninja -DCMAKE_BUILD_TYPE=Release ^
    -DCMAKE_C_COMPILER=%INTEL_ROOT%/bin/icx.exe ^
    -DCMAKE_CXX_COMPILER=%INTEL_ROOT%/bin/icx.exe ^
    -DCMAKE_EXE_LINKER_FLAGS=/Qoption,link,/LIBPATH:%INTEL_ROOT%/lib ^
    -DCMAKE_SHARED_LINKER_FLAGS=/Qoption,link,/LIBPATH:%INTEL_ROOT%/lib ^
    -DISOGEN_REQUIRE_INTELLLVM=ON
if errorlevel 1 exit /b 1

cmake --build "%BUILD_DIR%"
if errorlevel 1 exit /b 1

cmake --install "%BUILD_DIR%"
if errorlevel 1 exit /b 1

echo Intel build installed to %SRC_DIR%..\isogen\bin
endlocal
