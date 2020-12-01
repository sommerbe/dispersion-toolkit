## CMAKE file for
##
## dispersion toolkit
##
## Authors: 
## Benjamin Sommer: 2020 - current
## ========================================

project(dispersion-toolkit)

cmake_minimum_required(VERSION 3.4)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Set additional project information
set(COMPANY "Benjamin Sommer")
set(COPYRIGHT "Copyright (c) 2020 Benjamin Sommer. All rights reserved.")
set(IDENTIFIER "benjaminsommer@mailbox.org")

##############################################################
# Options
##############################################################

# set(CMAKE_INSTALL_PATH_APPS "/usr/share/applications/" CACHE PATH "Application Desktop File install directory")

# Find includes in corresponding build directories
# set(CMAKE_INCLUDE_CURRENT_DIR ON)

# install location
set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE BOOL "Installation Prefix")

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE CACHE BOOL "Add lib paths outside of build tree to install RPATH")

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

set(CMAKE_CXX_STANDARD 11)

set( VERSION_MAJOR 1)
set( VERSION_MINOR 0)
set( VERSION_PATCH 0)
set( VERSION_STRING ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

set(CONFIGURE_MANPAGE ON CACHE BOOL "Generates manpages and their embedding into the main.cpp files")

##############################################################
# ISA configuration
##############################################################

set(BUILD_ISA_DEVICE_TARGET_SSE2  OFF CACHE BOOL "Generate code path for SSE2")
set(BUILD_ISA_DEVICE_TARGET_SSE41 OFF CACHE BOOL "Generate code path for SSE41")
set(BUILD_ISA_DEVICE_TARGET_SSE42 OFF CACHE BOOL "Generate code path for SSE42")
set(BUILD_ISA_DEVICE_TARGET_AVX   OFF CACHE BOOL "Generate code path for AVX")
set(BUILD_ISA_DEVICE_TARGET_AVX2  OFF CACHE BOOL "Generate code path for AVX2")
set(BUILD_ISA_DEVICE_TARGET_NATIVE  OFF CACHE BOOL "Generate code path for NATIVE")

set(ISA_TARGETS "")
if(BUILD_ISA_DEVICE_TARGET_SSE2 STREQUAL "ON")
  set(ISA_TARGETS ${ISA_TARGETS} ${FLAGS_SSE2})
endif()

if(BUILD_ISA_DEVICE_TARGET_SSE41 STREQUAL "ON")
  set(ISA_TARGETS ${ISA_TARGETS} ${FLAGS_SSE41})
endif()

if(BUILD_ISA_DEVICE_TARGET_SSE42 STREQUAL "ON")
  set(ISA_TARGETS ${ISA_TARGETS} ${FLAGS_SSE42})
endif()

if(BUILD_ISA_DEVICE_TARGET_AVX STREQUAL "ON")
  set(ISA_TARGETS ${ISA_TARGETS} ${FLAGS_AVX})
endif()

if(BUILD_ISA_DEVICE_TARGET_AVX2 STREQUAL "ON")
  set(ISA_TARGETS ${ISA_TARGETS} ${FLAGS_AVX2})
endif()

if(BUILD_ISA_DEVICE_TARGET_NATIVE STREQUAL "ON")
  set(ISA_TARGETS ${ISA_TARGETS} ${FLAGS_NATIVE})
endif()

set(ISA_TARGETS "${ISA_TARGETS}")
string(REGEX REPLACE ";" "," ISA_TARGETS "${ISA_TARGETS}")

##############################################################
# Output paths
##############################################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_MAN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/man")

##############################################################
# Compiler Selection
##############################################################
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  include(gcc)
# elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#   include(clang)
# elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
#   include(icc)
# elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
#   include(msvc)
else()
  message(FATAL_ERROR "Unknown compiler specified: " ${COMPILER})
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the build type: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_VERBOSE_MAKEFILE true CACHE BOOL "Enables verbose mode.")

##############################################################
# Dependencies
##############################################################
# find_package(glfw REQUIRED)
find_package(OpenMP)

##############################################################
# Helper scripts and functions
##############################################################
function(configure_manpage_embedding NAME SECTION)
  set(MANPAGE "" CACHE STRING "manual preprint page")
  file(READ ${NAME}.1.pre MANPAGE)
  configure_file(${NAME}.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/${NAME}.hpp ESCAPE_QUOTES @ONLY)
endfunction()

##############################################################
# Compiler configuration
##############################################################
message(STATUS "Current CXX Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Current C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "Current CXX FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "Current LINKER FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

##############################################################
# Modules
##############################################################
set(src "${CMAKE_SOURCE_DIR}/src")

add_subdirectory(src/set/fibonacci)
add_subdirectory(src/set/kritzinger)
add_subdirectory(src/set/cswap)
add_subdirectory(src/measure/dispcombinatorial)
add_subdirectory(src/measure/dispgs)
add_subdirectory(src/opt/dispoptgs)
add_subdirectory(src/stat/confidence)
add_subdirectory(src/vis/psspy)
add_subdirectory(src/adapter/utk)