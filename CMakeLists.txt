# CMake minimum version and project name
cmake_minimum_required(VERSION 3.14)
project(tagger)

# ---- Configuration Options ----

# Option to enable tag array subsampling (default: OFF)
option(TAG_ARRAY_SUBSAMPLING "Enable tag array subsampling" OFF)
if(TAG_ARRAY_SUBSAMPLING)
    add_definitions(-DTAG_ARRAY_SUBSAMPLING)
    message(STATUS "Tag array subsampling enabled.")
else()
    message(STATUS "Tag array subsampling disabled.")
endif()

# Option to enable ground truth checks (default: OFF)
option(GROUND_TRUTH_CHECKS "Enable ground truth checks" OFF)
if(GROUND_TRUTH_CHECKS)
    add_definitions(-DGROUND_TRUTH_CHECKS)
    message(STATUS "Ground truth checks enabled.")
else()
    message(STATUS "Ground truth checks disabled.")
endif()

# Option to enable developer mode (default: OFF)
option(DEVELOPER_MODE "Enable developer mode" OFF)
if(DEVELOPER_MODE)
    add_definitions(-DDEVELOPER_MODE)
    message(STATUS "Developer mode enabled.")
else()
    message(STATUS "Developer mode disabled.")
endif()

# Compiler flags for warnings and C++ standard
set(COMMON_FLAGS "-Wall -pedantic -mpopcnt -std=gnu++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")

# Set default 32-bit or 64-bit type
set(THIRTY_TWO_DEFAULT OFF)

option(THIRTY_TWO "Use 32-bit types" ${THIRTY_TWO_DEFAULT})
if(THIRTY_TWO)
    add_definitions(-DTHIRTY_TWO)
    message(STATUS "Using 32-bit types.")
else()
    message(STATUS "Using 64-bit types.")
endif()

# Always define ALPHABET as 6
add_definitions(-DALPHABET=6)

# ---- Source Files ----

set(SOURCE_FILES
    src/indexhelpers.cpp
    src/nucleotide.cpp
    src/logger.cpp
    src/parameters/alignparameters.cpp
    src/indexinterface.cpp
)

# Add compression-specific source files
list(APPEND SOURCE_FILES src/bmove/bmove.cpp src/bmove/moverepr.cpp)

# Create the auxiliary library
add_library(AuxiliaryClassesTagging ${SOURCE_FILES})


# ---- Update Git Submodules ----

# Update the Big-BWT submodule
# Find Git
find_package(Git QUIET)

set(BIG_BWT_USABLE TRUE)


# Find Python 3.8 or greater
find_package(Python3 3.8 COMPONENTS Interpreter)

if (Python3_FOUND)

message(STATUS "Python executable detected by CMake: ${Python3_EXECUTABLE}")

# Check Python version
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import sys; print(f'{sys.version_info[0]} {sys.version_info[1]}')"
        OUTPUT_VARIABLE PYTHON_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Split the major and minor version
    string(REPLACE " " ";" PYTHON_VERSION_LIST "${PYTHON_VERSION}")
    list(GET PYTHON_VERSION_LIST 0 PYTHON_MAJOR_VERSION)
    list(GET PYTHON_VERSION_LIST 1 PYTHON_MINOR_VERSION)

    # Check if the Python version is >= 3.8
    if (PYTHON_MAJOR_VERSION GREATER 3 OR (PYTHON_MAJOR_VERSION EQUAL 3 AND PYTHON_MINOR_VERSION GREATER_EQUAL 8))
        message(STATUS "Found Python: ${Python3_EXECUTABLE} (Version ${PYTHON_MAJOR_VERSION}.${PYTHON_MINOR_VERSION})")

        # Check if psutil is installed
        execute_process(
                COMMAND ${Python3_EXECUTABLE} -c "import psutil"
                RESULT_VARIABLE PSUTIL_INSTALLED
                ERROR_QUIET
                OUTPUT_QUIET
        )

        if (NOT PSUTIL_INSTALLED EQUAL "0")
                message(WARNING "Python package 'psutil' is required for prefix-free parsing, but not installed. Please install it using pip (pip install psutil) if you want to use prefix-free parsing.")
                set(BIG_BWT_USABLE FALSE)
        else()
                message(STATUS "Python 3.8 or greater and psutil are installed.")
                set(BIG_BWT_USABLE TRUE)
        endif()

    else()
        message(WARNING "Found Python version ${PYTHON_MAJOR_VERSION}.${PYTHON_MINOR_VERSION}, but Python 3.8 or greater is required for prefix-free parsing.")
        set(BIG_BWT_USABLE FALSE)
    endif()
else()
    message(WARNING "Python not found, but is required for prefix-free parsing")
    set(BIG_BWT_USABLE FALSE)
endif()


# Proceed with submodule update if usable
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    if(BIG_BWT_USABLE)
        message(STATUS "Updating Big-BWT submodule")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init -q external/Big-BWT
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    else()
        message(WARNING "Prefix-free parsing is not usable due to missing dependencies (Python 3.8 or psutil).")
    endif()
endif()

# ---- Executables and Linking ----

add_executable(tagger src/taggingMEM.cpp src/fastq.cpp src/seqfile.cpp src/util.cpp)
target_link_libraries(tagger AuxiliaryClassesTagging)

add_executable(tagger_build src/logger.cpp src/buildindex.cpp src/bmove/moverepr.cpp src/parameters/buildparameters.cpp src/seqfile.cpp src/util.cpp)
if(BIG_BWT_USABLE)
    add_definitions(-DBIG_BWT_USABLE)
    # Build the Big-BWT submodule if compression is enabled
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/external/Big-BWT/bwtcheck
        COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_SOURCE_DIR}/external/Big-BWT make
        COMMENT "Building Big-BWT executable"
    )
    add_custom_target(bigbwt_target ALL DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/external/Big-BWT/bwtcheck)
    add_dependencies(tagger_build bigbwt_target)
endif()
# Copy the script to the build directory after tagger_build is built
add_custom_command(
    TARGET tagger_build POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/src/bmove/tagger_build_pfp.sh ${CMAKE_BINARY_DIR}/tagger_build_pfp.sh
    COMMENT "Copying tagger_build_pfp.sh to the build directory"
)

# ---- External Dependencies ----

# overwrite the message function to suppress messages if MESSAGE_QUIET is set
function(message)
    if (NOT MESSAGE_QUIET)
        _message(${ARGN})
    endif()
endfunction()
set(MESSAGE_QUIET OFF)

set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "Suppress developer warnings")
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)

# FetchContent module for including external projects
include(FetchContent)

# libsais library (added to tagger_build)

if(THIRTY_TWO)
    # Set specific build options for libsais
    message(STATUS "Fetching libsais library")
    set(LIBSAIS_USE_OPENMP OFF CACHE BOOL "Disable OpenMP support in libsais")
    set(LIBSAIS_BUILD_SHARED_LIB OFF CACHE BOOL "Build libsais as a static library")

    FetchContent_Declare(
    libsais
    GIT_REPOSITORY https://github.com/IlyaGrebnov/libsais.git
    GIT_TAG v2.8.5
    )
    FetchContent_MakeAvailable(libsais)
    target_link_libraries(tagger_build libsais)
else()
    find_package(Git QUIET)
    # update git submodule libdivsufsort
    message(STATUS "Updating libdivsufsort submodule")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init  external/libdivsufsort
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()


    set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build divsufsort64")
    set(BUILD_EXAMPLES OFF CACHE BOOL "Build examples")
   

    
    set(MESSAGE_QUIET ON)
    add_subdirectory(external/libdivsufsort)
    set(MESSAGE_QUIET OFF)

     # If needed, add dependencies to ensure libdivsufsort is built first
    add_dependencies(tagger_build divsufsort64)


    # Link against the libdivsufsort library
    target_link_libraries(tagger_build  divsufsort64)


    # Specify include directories for the headers that are generated during build
    # Set the include path to the appropriate build directory where the headers are located
    target_include_directories(tagger_build PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/external/libdivsufsort/include
    )

   

endif()

# parallel-hashmap library
message(STATUS "Fetching parallel-hashmap library")
FetchContent_Declare(
    parallel-hashmap
    GIT_REPOSITORY https://github.com/greg7mdp/parallel-hashmap.git
    GIT_TAG e5b892baed478513adcb6425773cae1eda033057
)
FetchContent_MakeAvailable(parallel-hashmap)
include_directories(${parallel-hashmap_SOURCE_DIR})

# fmt library
message(STATUS "Fetching fmt library")
FetchContent_Declare(
  fmtlib
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 10.2.1
)

# FetchContent module for including external projects
include(FetchContent)

# Fetch fmtlib from GitHub
FetchContent_Declare(
  fmtlib
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 10.2.1 # Adjust the version as needed
)

FetchContent_MakeAvailable(fmtlib)

target_link_libraries(AuxiliaryClassesTagging fmt::fmt)
target_link_libraries(tagger_build fmt::fmt)
target_link_libraries(tagger fmt::fmt)



# ---- Build Configuration ----

# Set default build type to Release if not specified
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# ----Installed Dependencies and Linking ----

# Check for zlib
find_package(ZLIB)
if(ZLIB_FOUND)
    add_definitions(-DHAVE_ZLIB)
    include_directories(${ZLIB_INCLUDE_DIR})
    target_link_libraries(tagger ${ZLIB_LIBRARY})
    target_link_libraries(tagger_build ${ZLIB_LIBRARY})
endif()

# Link Threads package
find_package(Threads REQUIRED)
target_link_libraries(tagger Threads::Threads)

# SDLS library 
# Set the path to your custom FindSDSL.cmake module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
find_package(SDSL)

if(SDSL_FOUND)
    message(STATUS "Found SDSL library at ${SDSL_LIBRARY}")
    message(STATUS "SDSL include directory: ${SDSL_INCLUDE_DIR}")

    # Include directories
    include_directories(${SDSL_INCLUDE_DIR})
    # Link the SDSL library
    target_link_libraries(tagger ${SDSL_LIBRARY})
    target_link_libraries(tagger_build ${SDSL_LIBRARY})
else()
    message(FATAL_ERROR "\nFATAL ERROR: SDSL library not found. Please install it and try again. If "
        "you installed SDSL in a non-standard location, you can point cmake to the "
        "installation location as follows: \ncmake -DSDSL_INCLUDE_DIR=<path-to-sdsl>/include .")
endif()

# ---- Installation ----

# Install targets to the bin directory
install(TARGETS tagger DESTINATION bin)

# # ---- Unit Tests ----

# # Add the unittest subdirectory
# add_subdirectory(unittest)
