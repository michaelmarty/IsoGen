# ======================== Setup ========================
cmake_minimum_required(VERSION 3.22.1)
include(GenerateExportHeader)

project(isogen C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_SHARED_LIBRARY_PREFIX "")
option(
    ISOGEN_STATIC_GNU_RUNTIME
    "Statically link libstdc++ and libgcc into Linux binaries for portability"
    ON
)
option(
    ISOGEN_ENABLE_AVX2
    "Build the runtime-dispatched AVX2/FMA neural-network implementation"
    ON
)
option(
    ISOGEN_REQUIRE_INTELLLVM
    "Fail configuration unless the Intel LLVM compiler is active"
    OFF
)

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

set(ISOGEN_AVX2_DISPATCH_AVAILABLE OFF)
if(ISOGEN_ENABLE_AVX2 AND
   CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|i[3-6]86)$")
    if(MSVC)
        set_source_files_properties(isogen_avx2.c PROPERTIES COMPILE_OPTIONS "/arch:AVX2")
        set(ISOGEN_AVX2_DISPATCH_AVAILABLE ON)
    elseif(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|IntelLLVM")
        set_source_files_properties(isogen_avx2.c PROPERTIES COMPILE_OPTIONS "-mavx2;-mfma")
        set(ISOGEN_AVX2_DISPATCH_AVAILABLE ON)
    else()
        message(WARNING "AVX2 runtime dispatch is not supported by ${CMAKE_C_COMPILER_ID}; using scalar code")
    endif()
endif()

if(ISOGEN_AVX2_DISPATCH_AVAILABLE)
    list(APPEND ISOGEN_SOURCES isogen_avx2.c)
    message(STATUS "Building runtime-dispatched AVX2/FMA support")
else()
    message(STATUS "Building scalar neural-network implementation")
endif()

add_library(isogen SHARED ${ISOGEN_SOURCES})
set_target_properties(isogen PROPERTIES LINKER_LANGUAGE C)
target_include_directories(isogen PUBLIC ${FFTW_DIR})
target_compile_definitions(isogen PUBLIC ISOGEN_BUILD_DLL)
if(ISOGEN_AVX2_DISPATCH_AVAILABLE)
    target_compile_definitions(isogen PRIVATE ISOGEN_AVX2_DISPATCH=1)
endif()


if(MSVC)
    add_library(fftw3 STATIC IMPORTED)
    set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION ${FFTW_DIR}/fftw3.lib)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
else()
    message(STATUS "Not MSVC")
    message(STATUS "CMAKE_C_COMPILER_ID: ${CMAKE_C_COMPILER_ID}")
    find_library(
        FFTW3_LIBRARY
        NAMES fftw3 libfftw3.so.3 fftw3.so.3
        REQUIRED
    )
    add_library(fftw3 SHARED IMPORTED)
    set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION "${FFTW3_LIBRARY}")
    target_link_libraries(isogen PUBLIC m)
endif()

target_link_libraries(isogen PUBLIC fftw3)
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








