#
#	This is the top level CMakeLists.txt file
#	It should be located within master/percept
#	But not in the BUILD folder itself
#
###########################################################################################################

cmake_minimum_required (VERSION 2.8.11)

###########################################################################################################
#
#  THESE PATHS MAY NEED TO BE CHANGED DEPENDING ON WHERE THEY ARE STORED!!!
#

SET(CMAKE_FIND_LIBRARY_PREFIXES lib)
SET(CMAKE_FIND_LIBRARY_SUFFIXES a)

MESSAGE("-- CMAKE_SOURCE_DIR = ${CMAKE_SOURCE_DIR}")

MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
SET(CMAKE_PREFIX_PATH ${Trilinos_DIR} ${CMAKE_PREFIX_PATH})
SET(CMAKE_PREFIX_PATH ${Trilinos_PATH} ${CMAKE_PREFIX_PATH})
find_package(Trilinos REQUIRED)

# This needs to be cleaned up
FIND_LIBRARY(YAML_LIBRARY NAMES libyaml-cpp.a PATHS ${YAML_DIR}/lib )
FIND_PATH(YAML_INCLUDES yaml.h PATHS ${YAML_DIR}/include/yaml-cpp )
INCLUDE_DIRECTORIES(${YAML_INCLUDES}/..)

# OpenNURBS
FIND_LIBRARY(OPENNURBS_LIBRARY NAMES libopenNURBS.a PATHS ${OPENNURBS_DIR} )
FIND_PATH(OPENNURBS_INCLUDES opennurbs.h PATHS ${OPENNURBS_DIR} )
INCLUDE_DIRECTORIES(${OPENNURBS_INCLUDES})

MESSAGE("-- OPENNURBS_DIR = ${OPENNURBS_DIR}")
MESSAGE("-- OPENNURBS_INCLUDES = ${OPENNURBS_INCLUDES}")
MESSAGE("-- OPENNURBS_LIBRARY = ${OPENNURBS_LIBRARY}")



# gtest
FIND_LIBRARY(GTEST_LIBRARY NAMES libgtest.a PATHS ${PERCEPT_TPLS_INSTALL_DIR}/lib NO_DEFAULT_PATH)
FIND_PATH(GTEST_INCLUDES gtest.h PATHS ${PERCEPT_TPLS_INSTALL_DIR}/include/gtest )
MESSAGE("-- GTEST_LIBRARY = ${GTEST_LIBRARY}")
MESSAGE("-- GTEST_INCLUDES = ${GTEST_INCLUDES}")
INCLUDE_DIRECTORIES(${GTEST_INCLUDES})

INCLUDE_DIRECTORIES(${STK_EXTRA_DIR})

SET(BOOST_LIBRARIES
   ${BOOST_DIR}/lib/libboost_exception.a
   ${BOOST_DIR}/lib/libboost_filesystem.a
   ${BOOST_DIR}/lib/libboost_graph.a
   ${BOOST_DIR}/lib/libboost_graph_parallel.a
   ${BOOST_DIR}/lib/libboost_mpi.a
   ${BOOST_DIR}/lib/libboost_program_options.a
   ${BOOST_DIR}/lib/libboost_regex.a
   ${BOOST_DIR}/lib/libboost_serialization.a
   ${BOOST_DIR}/lib/libboost_signals.a
   ${BOOST_DIR}/lib/libboost_system.a
   ${BOOST_DIR}/lib/libboost_thread.a
   ${BOOST_DIR}/lib/libboost_wserialization.a
)

SET(BOOST_INCLUDES ${BOOST_DIR}/include )
INCLUDE_DIRECTORIES(${BOOST_INCLUDES})
MESSAGE("-- BOOST_INCLUDES = ${BOOST_INCLUDES}")
MESSAGE("-- BOOST_LIBRARIES = ${BOOST_LIBRARIES}")

#include(FindPackageHandleStandardArgs)

INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Trilinos_TPL_INCLUDE_DIRS})

#######################################     TRILINOS	 ##################################################

MESSAGE("\nFound Trilinos!  Here are the details: ")
# Uncomment lines to view more details about Trilinos
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_BIN_DIRS = ${Trilinos_BIN_DIRS}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
MESSAGE("End of Trilinos details\n")

# Optional Installation helpers
SET (INSTALL_PERCEPT FALSE)
IF (ENABLE_INSTALL)
    SET (INSTALL_PERCEPT TRUE)
ENDIF()

IF (INSTALL_PERCEPT)
  set(BINARY_INSTALL_DIR bin)
  set(INCLUDE_INSTALL_DIR include)
  set(LIB_INSTALL_DIR lib)
  INCLUDE(CMakePackageConfigHelpers)
ENDIF ()

MESSAGE("Setting and checking of compilers:")
SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

IF (CMAKE_Fortran_COMPILER) # Enable Fortran if it is enabled in Trilinos.
  ENABLE_LANGUAGE(Fortran)
ENDIF()

# Build Percept as shared libraries if Trilinos was compiled that way

IF(Trilinos_BUILD_SHARED_LIBS)
  SET(Percept_LIBRARY_TYPE SHARED)
  SET(BUILD_SHARED_LIBS ON)
  MESSAGE("-- Building Percept with shared libraries")
ELSE(Trilinos_BUILD_SHARED_LIBS)
  SET(Percept_LIBRARY_TYPE STATIC)
ENDIF(Trilinos_BUILD_SHARED_LIBS)

#######################################     PERCEPT     ##########################

PROJECT(Percept)

SET(GCC_NO_MPI_WARNINGS_FLAGS "-Wno-literal-suffix")   # C++ 11 warnings originating with MPI

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_NO_MPI_WARNINGS_FLAGS} -Wno-deprecated-declarations")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Trilinos_CXX_COMPILER_FLAGS}")
SET(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} ${Trilinos_Fortran_COMPILER_FLAGS})
MESSAGE("-- CMAKE_CXX_FLAGS     = ${CMAKE_CXX_FLAGS}")
MESSAGE("-- CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")

IF (STK_ADAPT_HAVE_YAML_CPP)
  ADD_DEFINITIONS(-DSTK_ADAPT_HAVE_YAML_CPP=1)
  MESSAGE("-- using YAML in Percept")
ENDIF()

SET(HEADERS "")
SET(SOURCES "")

FILE (GLOB ADAPT_MAIN

     src/adapt/main/*.cpp

)

FILE (GLOB SOURCE

     src/percept/*.cpp
     src/percept/function/*.cpp
     src/percept/function/internal/*.cpp
     src/percept/norm/*.cpp
     src/percept/mesh/gen/*.cpp

     src/percept/element/intrepid/*.cpp
     src/percept/mesh/mod/smoother/*.cpp
     src/percept/mesh/geometry/kernel/*.cpp
     src/percept/mesh/geometry/kernel/search/*.cpp
     src/percept/mesh/geometry/kernel/xfer/*.cpp

     src/percept/mesh/geometry/recovery/*.cpp
     src/percept/mesh/geometry/volume/*.cpp
     src/percept/mesh/geometry/volume/sierra_only/*.cpp
     src/percept/mesh/geometry/stk_geom/*.cpp
     src/percept/mesh/geometry/stk_geom/3D/*.cpp
     src/percept/structured/*.cpp
     src/percept/xfer/*.cpp
     src/percept/util/*.cpp
     src/percept/stk_rebalance/*.cpp
     src/percept/stk_rebalance_utils/*.cpp
     src/percept/fixtures/*.cpp

     src/adapt/*.cpp
     src/adapt/sierra_element/*.cpp
     src/adapt/geometry/*.cpp
     src/adapt/markers/*.cpp

)

FILE (GLOB HEADER

     src/percept/*.hpp
     src/percept/function/*.hpp
     src/percept/function/internal/*.hpp
     src/percept/norm/*.hpp
     src/percept/mesh/gen/*.hpp

     src/percept/element/intrepid/*.hpp
     src/percept/mesh/mod/smoother/*.hpp
     src/percept/mesh/geometry/kernel/*.hpp
     src/percept/mesh/geometry/kernel/search/*.hpp
     src/percept/mesh/geometry/kernel/xfer/*.hpp

     src/percept/mesh/geometry/recovery/*.hpp
     src/percept/mesh/geometry/volume/*.hpp
     src/percept/mesh/geometry/volume/sierra_only/*.hpp
     src/percept/mesh/geometry/stk_geom/*.hpp
     src/percept/mesh/geometry/stk_geom/3D/*.hpp
     src/percept/structured/*.hpp
     src/percept/xfer/*.hpp
     src/percept/util/*.hpp
     src/percept/stk_rebalance/*.hpp
     src/percept/stk_rebalance_utils/*.hpp
     src/percept/fixtures/*.hpp

     src/adapt/*.hpp
     src/adapt/sierra_element/*.hpp
     src/adapt/geometry/*.hpp
     src/adapt/markers/*.hpp
)

MESSAGE("-- SOURCE = ${SOURCE}")
MESSAGE("-- HEADER = ${HEADER}")

FILE (GLOB ADAPT_SOURCE src/adapt/main/*.cpp)
MESSAGE("-- ADAPT_SOURCE = ${ADAPT_SOURCE}")

FILE (GLOB PERCEPT_UTEST_SOURCE test/unit_tests/*.cpp)
MESSAGE("-- PERCEPT_UTEST_SOURCE = ${PERCEPT_UTEST_SOURCE}")

FILE (GLOB PERCEPT_HTEST_SOURCE test/heavy_tests/*.cpp)
MESSAGE("-- PERCEPT_HTEST_SOURCE = ${PERCEPT_HTEST_SOURCE}")

FILE (GLOB PERCEPT_RTEST_SOURCE test/regression_tests/*.cpp)
MESSAGE("-- PERCEPT_RTEST_SOURCE = ${PERCEPT_RTEST_SOURCE}")

FILE (GLOB PERCEPT_UTEST_HEADER test/unit_tests/*.hpp)
MESSAGE("-- PERCEPT_UTEST_HEADER = ${PERCEPT_UTEST_HEADER}")

FILE (GLOB PERCEPT_RTEST_HEADER test/regression_tests/*.hpp)
MESSAGE("-- PERCEPT_RTEST_HEADER = ${PERCEPT_RTEST_HEADER}")

INCLUDE_DIRECTORIES (${CMAKE_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES (${CMAKE_SOURCE_DIR}/test)

SET(percept_lib_name "percept")

#SET(percept_ex_name "perceptX")
SET(adapt_ex_name "mesh_adapt")
SET(percept_utest_name "percept_utest")
SET(percept_rtest_name "percept_rtest")
SET(percept_htest_name "percept_htest")

MESSAGE("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
IF (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
   #SET (percept_ex_name "perceptXd")
#   ADD_DEFINITIONS("-Wall" "-Werror")
   SET(CMAKE_CXX_FLAGS "-g  -std=c++11")
   MESSAGE("-- Debug Build, CMAKE_CXX_FLAGS= ${CMAKE_CXX_FLAGS}")
ENDIF()

IF (NOT STK_PERCEPT_LITE)

  ADD_EXECUTABLE(${adapt_ex_name} ${ADAPT_MAIN})
  TARGET_LINK_LIBRARIES(${adapt_ex_name} ${percept_lib_name})

  ADD_EXECUTABLE(${percept_utest_name} ${PERCEPT_UTEST_SOURCE} ${PERCEPT_UTEST_HEADER})
  TARGET_LINK_LIBRARIES(${percept_utest_name} ${percept_lib_name} ${GTEST_LIBRARY})

  ADD_EXECUTABLE(${percept_htest_name} ${PERCEPT_HTEST_SOURCE} ${PERCEPT_RTEST_HEADER})
  TARGET_LINK_LIBRARIES(${percept_htest_name} ${percept_lib_name} ${GTEST_LIBRARY})

  ADD_EXECUTABLE(${percept_rtest_name} ${PERCEPT_RTEST_SOURCE} ${PERCEPT_RTEST_HEADER})
  TARGET_LINK_LIBRARIES(${percept_rtest_name} ${percept_lib_name} ${GTEST_LIBRARY})

  ADD_LIBRARY (${percept_lib_name} ${SOURCE} ${HEADER})
  TARGET_LINK_LIBRARIES(${percept_lib_name} ${Trilinos_LIBRARIES} ${YAML_LIBRARY} ${OPENNURBS_LIBRARY} ${BOOST_LIBRARIES})

  SET_PROPERTY(
    TARGET ${percept_lib_name} ${adapt_ex_name} ${percept_utest_name} ${percept_rtest_name}
    PROPERTY COMPILE_DEFINITIONS STK_PERCEPT_LITE=0 STK_PERCEPT_HAS_GEOMETRY STK_PERCEPT_USE_INTREPID
  )

  SET_PROPERTY(
    TARGET ${percept_htest_name} ${percept_rtest_name} ${percept_utest_name}
    PROPERTY COMPILE_DEFINITIONS STK_PERCEPT_LITE=0
  )

ENDIF()


MESSAGE("\nAnd CMake says...:")
