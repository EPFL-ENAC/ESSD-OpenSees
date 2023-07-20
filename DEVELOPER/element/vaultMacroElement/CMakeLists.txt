cmake_minimum_required(VERSION 3.0)

# Set project name and version
project(VaultMacroElement_2D VERSION 1.0)

# Set the source files for the DLL
set(SRC_FILES VaultMacroElement_2D.cpp VaultMacroElement_2D.h)

# Set the DLL output directory
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build)

# Headers and sources
set(OPENSEES_INCLUDE_DIRS ../../core)
file(GLOB OPENSEES_SOURCES ../../core/*.cpp)

# Find Python
find_package(Python COMPONENTS Interpreter Development)
if(NOT Python_FOUND)
	message(FATAL_ERROR "Python not found")
endif()
include_directories(${Python_INCLUDE_DIRS})

# Set the target to be a shared library
add_library(VaultMacroElement_2D SHARED ${SRC_FILES})

# Set the include directories
target_include_directories(VaultMacroElement_2D PRIVATE ${OPENSEES_INCLUDE_DIRS})
target_sources(VaultMacroElement_2D PRIVATE ${OPENSEES_SOURCES})
target_link_libraries(VaultMacroElement_2D PRIVATE ${Python_LIBRARIES})
