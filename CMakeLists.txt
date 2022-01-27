cmake_minimum_required (VERSION 3.0)
project(sz_cpp VERSION 1.0.11)

include(CTest)
include(GNUInstallDirs)

#correct was to set a default build type
# https://blog.kitware.com/cmake-and-the-default-build-type/
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type was set. Setting build type to ${default_build_type}.")
  set(CMAKE_BUILD_TYPE ${default_build_type} CACHE 
    STRING "Choose the type to build" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
option(BUILD_SHARED_LIBS "build libpressio as a shared library" ON)

find_package(PkgConfig REQUIRED)
pkg_search_module(ZSTD IMPORTED_TARGET GLOBAL libzstd)

add_library (sz_cpp  
	src/sz_compression_utils.cpp
	src/sz_compress_cp_preserve_2d.cpp
	src/sz_compress_cp_preserve_3d.cpp
	src/sz_compress_3d.cpp
	src/sz_compress_pwr.cpp
	src/sz_decompression_utils.cpp
	src/sz_decompress_cp_preserve_2d.cpp
	src/sz_decompress_cp_preserve_3d.cpp
	src/sz_decompress_3d.cpp
	src/sz_decompress_pwr.cpp
	src/sz_huffman.cpp
	src/sz_lossless.cpp)
target_compile_features(sz_cpp PUBLIC cxx_std_11)
target_include_directories(sz_cpp PUBLIC 
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
target_link_libraries(sz_cpp PRIVATE PkgConfig::ZSTD)

file(GLOB header_files ${PROJECT_SOURCE_DIR}/include/*.hpp)
install(FILES ${header_files} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/sz_cpp)
install (TARGETS sz_cpp EXPORT SZCppTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(EXPORT SZCppTargets NAMESPACE SZCpp:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZCpp/)
include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/SZCppConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/SZCppConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZCpp
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/SZCppConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY AnyNewerVersion
)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SZCppConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SZCppConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SZCpp
)

if(BUILD_TESTING)
  add_subdirectory (test)
endif()
