# ======================== Setup ========================
cmake_minimum_required(VERSION 3.22.1)
include(GenerateExportHeader)

project(isogen C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_SHARED_LIBRARY_PREFIX "")
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(
        CMAKE_INSTALL_PREFIX
        "${CMAKE_CURRENT_SOURCE_DIR}/../isogen/bin"
        CACHE PATH "Install path prefix, prepended onto install directories." FORCE
    )
endif()
option(
    ISOGEN_STATIC_GNU_RUNTIME
    "Statically link libstdc++ and libgcc into Linux binaries for portability"
    ON
)
option(
    ISOGEN_REQUIRE_INTELLLVM
    "Fail configuration unless the Intel LLVM compiler is active"
    OFF
)
set(FFTW3_ROOT "" CACHE PATH "Root directory of an FFTW3 installation")
set(FFTW3_DIR "" CACHE PATH "FFTW3 installation or library directory")

if(ISOGEN_REQUIRE_INTELLLVM AND NOT CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
    message(FATAL_ERROR
        "ISOGEN_REQUIRE_INTELLLVM is ON, but the C compiler is ${CMAKE_C_COMPILER_ID}"
    )
endif()

# Set SRC_DIR to the path of the IsoDec source code
get_filename_component(SRC_DIR ./ ABSOLUTE)
message(STATUS "SRC_DIR: ${SRC_DIR}")

# FFTW directory is src_dir/fftw3
set(FFTW_DIR "${SRC_DIR}/fftw")
message(STATUS "FFTW_DIR: ${FFTW_DIR}")

include_directories("${SRC_DIR}")
include_directories("${FFTW_DIR}")

#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath='$ORIGIN'")
# ======================== IsoGen Builds ========================

set(ISOGEN_SOURCES
                isogendep.c isogendep.h
                isogen_models.c isogen_models.h
                ${FFTW_DIR}/fftw3.h
                isogenpep.c isogenpep.h
                isogenrna.c isogenrna.h
                isogenatom.c isogenatom.h
)

add_library(isogen SHARED ${ISOGEN_SOURCES})
set_target_properties(isogen PROPERTIES LINKER_LANGUAGE C)
target_include_directories(isogen PUBLIC ${FFTW_DIR})
target_compile_definitions(isogen PUBLIC ISOGEN_BUILD_DLL)


set(ISOGEN_X64 OFF)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$")
    set(ISOGEN_X64 ON)
endif()

# The bundled FFTW archive contains AMD64 COFF objects and is valid only for
# native x64 MSVC builds. All other targets discover an architecture-matched
# FFTW installation supplied by the system, Homebrew, vcpkg, or an explicit
# FFTW3_ROOT/FFTW3_DIR hint.
if(MSVC AND ISOGEN_X64)
    add_library(fftw3 STATIC IMPORTED)
    set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION ${FFTW_DIR}/fftw3.lib)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
else()
    set(FFTW3_HINTS)
    foreach(FFTW3_HINT IN ITEMS
            "${FFTW3_ROOT}"
            "${FFTW3_DIR}"
            "$ENV{FFTW3_ROOT}"
            "$ENV{FFTW3_DIR}")
        if(FFTW3_HINT)
            list(APPEND FFTW3_HINTS "${FFTW3_HINT}")
        endif()
    endforeach()
    find_library(
        FFTW3_LIBRARY
        NAMES fftw3 libfftw3-3 libfftw3.so.3 fftw3.so.3
        HINTS ${FFTW3_HINTS}
        PATH_SUFFIXES lib lib64
    )
    if(NOT FFTW3_LIBRARY)
        message(FATAL_ERROR
            "FFTW3 was not found for ${CMAKE_SYSTEM_NAME}/${CMAKE_SYSTEM_PROCESSOR}. "
            "Install libfftw3-dev (Debian/Ubuntu), brew install fftw (macOS), "
            "or vcpkg install fftw3:<triplet> (Windows), then set FFTW3_ROOT, "
            "FFTW3_DIR, FFTW3_LIBRARY, or CMAKE_PREFIX_PATH."
        )
    endif()
    message(STATUS "FFTW3_LIBRARY: ${FFTW3_LIBRARY}")
    add_library(fftw3 UNKNOWN IMPORTED)
    set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION "${FFTW3_LIBRARY}")
    if(NOT WIN32)
        target_link_libraries(isogen PUBLIC m)
    endif()
    if(MSVC)
        add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
    endif()
endif()

target_link_libraries(isogen PUBLIC fftw3)

# Dynamic vcpkg installations need their FFTW DLL beside isogen.dll. Static
# installations have no matching runtime file and simply skip this step.
if(WIN32 AND NOT (MSVC AND ISOGEN_X64))
    find_file(
        FFTW3_RUNTIME_LIBRARY
        NAMES fftw3.dll libfftw3-3.dll
        HINTS ${FFTW3_HINTS}
        PATH_SUFFIXES bin
    )
    if(FFTW3_RUNTIME_LIBRARY)
        message(STATUS "FFTW3_RUNTIME_LIBRARY: ${FFTW3_RUNTIME_LIBRARY}")
        install(FILES "${FFTW3_RUNTIME_LIBRARY}" DESTINATION . COMPONENT Runtime)
    endif()
endif()
install(
    TARGETS isogen
    RUNTIME DESTINATION . COMPONENT Runtime
    LIBRARY DESTINATION . COMPONENT Runtime
    ARCHIVE DESTINATION . COMPONENT Development
)

# IntelLLVM can emit calls into Intel's optimized math runtime on Windows.
# Bundle those redistributable DLLs beside isogen.dll so the wheel remains
# self-contained on machines that do not have oneAPI installed.
if(WIN32 AND CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
    get_filename_component(ISOGEN_INTEL_COMPILER_BIN "${CMAKE_C_COMPILER}" DIRECTORY)
    foreach(ISOGEN_INTEL_RUNTIME IN ITEMS libmmd.dll svml_dispmd.dll)
        set(ISOGEN_INTEL_RUNTIME_PATH
            "${ISOGEN_INTEL_COMPILER_BIN}/${ISOGEN_INTEL_RUNTIME}"
        )
        if(NOT EXISTS "${ISOGEN_INTEL_RUNTIME_PATH}")
            message(FATAL_ERROR
                "Required Intel runtime was not found: ${ISOGEN_INTEL_RUNTIME_PATH}"
            )
        endif()
        install(
            FILES "${ISOGEN_INTEL_RUNTIME_PATH}"
            DESTINATION .
            COMPONENT Runtime
        )
    endforeach()
endif()

add_executable(isogen_test isogen_test.cpp)
target_include_directories(isogen_test PUBLIC ${SRC_DIR})
target_link_libraries(isogen_test PUBLIC isogen)

# Prebuilt Linux artifacts must not inherit the build host's libstdc++ symbol
# requirements. Otherwise a build made with a recent GCC cannot be loaded on
# older supported distributions or in many conda environments.
if(ISOGEN_STATIC_GNU_RUNTIME AND UNIX AND NOT APPLE AND
   CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|IntelLLVM")
    target_link_options(isogen PRIVATE -static-libstdc++ -static-libgcc)
    target_link_options(isogen_test PRIVATE -static-libstdc++ -static-libgcc)
endif()

install(
    TARGETS isogen_test
    RUNTIME DESTINATION bin
    COMPONENT Development
)







