cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(MEUtils)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Core Eve EG Geom Ged RGL)
include(${ROOT_USE_FILE})
file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_library(MEUtils SHARED ${SOURCES})
target_include_directories(MEUtils
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/MEUtils>)
target_link_libraries(MEUtils PRIVATE ${ROOT_LIBRARIES})
set_property(TARGET MEUtils PROPERTY CXX_STANDARD 11)
set(PUBLIC_HEADERS ${HEADERS})
set_target_properties(MEUtils
  PROPERTIES
  PUBLIC_HEADER "${PUBLIC_HEADERS}")
install(TARGETS MEUtils
  EXPORT MEUtilsTargets
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/MEUtils)
install(EXPORT MEUtilsTargets
  FILE MEUtilsTargets.cmake
  DESTINATION lib/cmake)
include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/MEUtilsConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/MEUtilsConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/example"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/MEUtilsConfigVersion.cmake"
  VERSION "${MEUtils_VERSION_MAJOR}.${MEUtils_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/MEUtilsConfig.cmake
  DESTINATION lib/cmake)
export(EXPORT MEUtilsTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/MEUtilsTargets.cmake")
export(PACKAGE MEUtils)
