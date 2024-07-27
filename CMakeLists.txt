cmake_minimum_required( VERSION 3.6 )
project( SIDIS )
set(CMAKE_CXX_FLAGS "-std=c++17")

find_package(ROOT REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(hipo4 REQUIRED IMPORTED_TARGET hipo4)
#set clas12 env
set (CLAS12ROOT $ENV{CLAS12ROOT})

SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_SOURCE_DIR}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_SOURCE_DIR}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_SOURCE_DIR}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

#set path for include directories and libraries
set (INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${HIPO}/include)
set (LIBRARIES ${HIPO}/lib)

list (APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
include(${ROOT_USE_FILE})


ADD_DEFINITIONS( -D_DATA="${CMAKE_SOURCE_DIR}/data" )


add_subdirectory(classes)

add_subdirectory(skimmers)
add_subdirectory(corrections)
add_subdirectory(analysis)
add_subdirectory(rhoAnalysis)
add_subdirectory(kaonAnalysis)
