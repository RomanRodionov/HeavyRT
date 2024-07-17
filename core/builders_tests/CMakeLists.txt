cmake_minimum_required (VERSION 3.8)
project (cbvh_test)

# get rid of annoying MSVC warnings.
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

IF(CMAKE_BUILD_TYPE MATCHES Debug)
  message("Debug build.")
ELSEIF(CMAKE_BUILD_TYPE MATCHES Release)
  message("Release build.")
  add_definitions(-DNDEBUG)
ELSE()
  message("Other build type.")
ENDIF()

set (CMAKE_CXX_STANDARD 17)

option(USE_EMBREE "Switch implementation to Embree" OFF)

include_directories(".")
include_directories("../common")
include_directories("../builders")
include_directories("../../examples/loader")

IF (USE_EMBREE)
  include_directories("../embree/include")
  include_directories("../embree")
  link_directories("../embree/build")
ENDIF()

if(WIN32)
  add_definitions(-D_ENABLE_EXTENDED_ALIGNED_STORAGE)
  link_directories(${ADDITIONAL_LIBRARY_DIRS})
  IF (USE_EMBREE)
    link_directories("embreewin64")
  ENDIF()
endif()


if (CMAKE_SYSTEM_PROCESSOR MATCHES "arm" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpu=neon")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
endif()

#uncomment this to detect broken memory problems via gcc sanitizers
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize-address-use-after-scope -fno-omit-frame-pointer -fsanitize=leak -fsanitize=undefined -fsanitize=bounds-strict")

IF (USE_EMBREE)
  add_definitions(-DCBVH_USE_EMBREE)
  set(SOURCES
    ../builders/cbvh_embree.cpp
    ../builders/cbvh_core.cpp
    ../builders/cbvh_esc.cpp
    ../builders/cbvh_esc2.cpp
    ../builders/cbvh_split.cpp
    ../builders/EXT_TriBoxOverlap.cpp
    ../builders/cbvh_hlbvh.cpp
    ../../examples/loader/cmesh.cpp
    cbvh_metrics.cpp
    cbvh_metrics_ray.cpp
    tests_traversal.cpp)
ELSE()
  set(SOURCES
    ../builders/cbvh_core.cpp
    ../builders/cbvh_esc.cpp
    ../builders/cbvh_esc2.cpp
    ../builders/cbvh_split.cpp
    ../builders/EXT_TriBoxOverlap.cpp
    ../builders/cbvh_hlbvh.cpp
    ../../examples/loader/cmesh.cpp
    cbvh_metrics.cpp
    cbvh_metrics_ray.cpp
    tests_traversal.cpp)
ENDIF()

set(SOURCES "${SOURCES}" main_tests.cpp)

add_executable(cbvh_test "${SOURCES}")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(cbvh_test PUBLIC OpenMP::OpenMP_CXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

set_target_properties(cbvh_test PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
IF (USE_EMBREE)
  IF (WIN32)
    target_link_libraries(cbvh_test LINK_PUBLIC embree3 embree_sse42 sys tasking simd lexers)
  ELSE()
    target_link_libraries(cbvh_test LINK_PUBLIC embree3 embree_sse42 sys tasking simd lexers pthread)
  ENDIF()
ENDIF()