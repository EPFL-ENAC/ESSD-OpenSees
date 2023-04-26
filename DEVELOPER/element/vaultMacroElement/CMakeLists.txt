cmake_minimum_required(VERSION 3.0)

# Set project name and version
project(VaultMacroElement VERSION 1.0)

# Set the source files for the DLL
set(SRC_FILES VaultMacroElement.cpp VaultMacroElement.h)

# Set the DLL output directory
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build)

# Headers and sources
set(OPENSEES_INDLUDE_DIRS ../../core)
file(GLOB OPENSEES_SOURCES ../../core/*.cpp)

# Set the target to be a shared library
add_library(VaultMacroElement SHARED ${SRC_FILES})

# Set the include directories
target_include_directories(VaultMacroElement PRIVATE ${OPENSEES_INDLUDE_DIRS})
target_sources(VaultMacroElement PRIVATE ${OPENSEES_SOURCES})
