cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(MEUtils)
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "~/install/MatrixElementUtils")
endif()
message("-- INSTALL PREFIX: " ${CMAKE_INSTALL_PREFIX})
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})
file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/inc/*.h)
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/me-utils.sh.in me-utils.sh)
add_library(MEUtils SHARED ${SOURCES})
target_include_directories(MEUtils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/inc)
target_link_libraries(MEUtils ${ROOT_LIBRARIES})
set(PUBLIC_HEADERS ${HEADERS})
set_target_properties(MEUtils
  PROPERTIES
  PUBLIC_HEADER "${PUBLIC_HEADERS}"
  RESOURCE "${PROJECT_BINARY_DIR}/me-utils.sh")
install(TARGETS MEUtils
		LIBRARY DESTINATION lib
		PUBLIC_HEADER DESTINATION inc
		RESOURCE DESTINATION env COMPONENT base)
