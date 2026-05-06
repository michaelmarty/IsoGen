# ======================== Setup ========================
cmake_minimum_required(VERSION 3.22.1)
include(GenerateExportHeader)

project(isogen C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_SHARED_LIBRARY_PREFIX "")

get_filename_component(INSTALL_DIR ../bin ABSOLUTE)
message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")

# Set SRC_DIR to the path of the IsoDec source code
get_filename_component(SRC_DIR ./ ABSOLUTE)
message(STATUS "SRC_DIR: ${SRC_DIR}")

# FFTW directory is src_dir/fftw3
set(FFTW_DIR "${SRC_DIR}/fftw")
message(STATUS "FFTW_DIR: ${FFTW_DIR}")

include_directories("${SRC_DIR}")
include_directories("${FFTW_DIR}")

#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath='$ORIGIN'")
set(CMAKE_INSTALL_RPATH ${INSTALL_DIR})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Set _CRT_SECURE_NO_WARNINGS option
if(MSVC)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()


# ======================== IsoGen Builds ========================

add_library(isogen SHARED
                isogendep.c isogendep.h
                ${FFTW_DIR}/fftw3.h
                isogenpep.c isogenpep.h
                isogenrna.c isogenrna.h
        isogen_test.cpp
)
set_target_properties(isogen PROPERTIES LINKER_LANGUAGE C)
target_include_directories(isogen PUBLIC ${FFTW_DIR})
target_compile_definitions(isogen PUBLIC ISOGEN_BUILD_DLL)


if(MSVC)
    add_library(fftw3 STATIC IMPORTED)
    set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION ${FFTW_DIR}/fftw3.lib)
else()
    message(STATUS "Not MSVC")
    message(STATUS "CMAKE_C_COMPILER_ID: ${CMAKE_C_COMPILER_ID}")
    #    set(FFTW_LIB "/usr/lib/x86_64-linux-gnu/libfftw3.so")
    #target_link_libraries(isogendep PUBLIC ${FFTW_LIB})
    add_library(fftw3 SHARED IMPORTED)
    set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION /usr/lib/x86_64-linux-gnu/libfftw3.so)
    target_link_libraries(isogen PUBLIC m)
endif()

target_link_libraries(isogen PUBLIC fftw3)
install(TARGETS isogen DESTINATION ${INSTALL_DIR})

add_executable(isogen_test isogen_test.cpp)
target_include_directories(isogen_test PUBLIC ${SRC_DIR})
target_link_libraries(isogen_test PUBLIC isogen isogen)

install(TARGETS isogen_test DESTINATION ${INSTALL_DIR})









