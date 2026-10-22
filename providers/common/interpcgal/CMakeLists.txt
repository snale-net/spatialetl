cmake_minimum_required(VERSION 3.10.2)
project(interpcgal)

include(CheckIncludeFileCXX)
check_include_file_cxx("cstdint" HAVE_CSTDINT)
if(NOT HAVE_CSTDINT)
    message(FATAL_ERROR "C++ compiler must support <cstdint>")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_BUILD_TYPE Release)

# Try to import all Python components potentially needed by nanobind
find_package(Python 3.9
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)
# Import nanobind through CMake's find_package mechanism
find_package(nanobind CONFIG REQUIRED)

set(SOURCE_DIR "src/")
include_directories(${SOURCE_DIR})

# We are now ready to compile the actual extension module
nanobind_add_module(
  # Name of the extension
  nninterpol

  # Target the stable ABI for Python 3.12+, which reduces
  # the number of binary wheels that must be built. This
  # does nothing on older Python versions
  STABLE_ABI

  # Source code goes here
 ${SOURCE_DIR}/nninterpol.cpp
)

# Install directive for scikit-build-core
install(TARGETS nninterpol LIBRARY DESTINATION spatialetl/operator/interpolator/interpcgal)

find_package(CGAL CONFIG REQUIRED)
target_link_libraries(nninterpol PRIVATE CGAL::CGAL)

target_compile_definitions(nninterpol PRIVATE
    -D__STDC_LIMIT_MACROS
    -D__STDC_CONSTANT_MACROS
)