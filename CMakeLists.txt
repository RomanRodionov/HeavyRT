cmake_minimum_required(VERSION 3.15)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

project(loader_app CXX)
add_compile_definitions(HALFFLOAT)

add_subdirectory(external/LiteMath)

add_library(project_warnings INTERFACE)
include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

add_library(HeavyRT)
add_executable(render_app)

if(CMAKE_SYSTEM_NAME STREQUAL Windows)
    link_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/embree/lib_win64)
    set(RAYTRACING_EMBREE_LIBS embree3)
else()
    find_package(Threads REQUIRED)
    target_link_directories(HeavyRT PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/external/embree/lib)
    target_link_directories(render_app PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/embree/lib)
    set(RAYTRACING_EMBREE_LIBS embree3 embree_sse42 embree_avx embree_avx2 lexers simd sys tasking) # embree_avx embree_avx2
endif()


set(HeavyRT_LIB_INCLUDES 
        external
        external/embree
        external/nanort
        external/tiny_gltf
        examples)

target_include_directories(HeavyRT PUBLIC ${HeavyRT_LIB_INCLUDES})
target_include_directories(HeavyRT PUBLIC core)
target_include_directories(HeavyRT PUBLIC core/common)
target_include_directories(HeavyRT PUBLIC core/builders)
target_include_directories(render_app PRIVATE ${HeavyRT_LIB_INCLUDES})
target_include_directories(render_app PRIVATE examples/loader)


find_package(OpenMP)

set(SCENE_LOADER_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/examples/loader/pugixml.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/examples/loader/hydraxml.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/examples/loader/cmesh.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/examples/loader/gltf_loader.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/examples/loader/gltf_loader_def.cpp)


set(CORE_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/core/common/Image2d.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/common/Timer.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/common/qmc_sobol_niederreiter.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/FactoryRT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/EmbreeRT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BruteForceRT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2CommonRT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2CommonRT_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2CommonLoftRT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2CommonLoftRT_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2CommonLoft16RT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2CommonLoft16RT_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2CommonLRFT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2CommonLRFT_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2CommonRTStacklessLBVH_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2CommonRTStacklessLBVH.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH4CommonRT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH4CommonRT_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH4HalfRT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH4HalfRT_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2Stackless.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2Stackless_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2FatRT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2FatRT_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2Fat16RT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2Fat16RT_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2FatCompactRT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/BVH2FatCompactRT_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/NanoRT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/NonunifNodeStorage.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/ShtRT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/ShtRT_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/ShtRT64.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/ShtRT64_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/Sht4NodesRT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/Sht4NodesRT_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/Sht4NodesHalfRT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/Sht4NodesHalfRT_host.cpp
        )

set(EXAMPLES_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/examples/01_eye_rays/eye_ray.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/examples/01_eye_rays/eye_ray_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/examples/02_rtao/rtao.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/examples/02_rtao/rtao_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/examples/RenderFactory.cpp
        )

set(BUILDERS_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/core/builders/lbvh.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/builders/lbvh_host.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/builders/cbvh.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/builders/cbvh_fat.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/builders/cbvh_core.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/builders/cbvh_embree2.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/builders/cbvh_split.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/builders/EXT_TriBoxOverlap.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/builders/cbvh_esc.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/builders/cbvh_esc2.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/builders/cbvh_hlbvh.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/builders/FatBVHReorder.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/builders/FatBVHClusterize.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/builders/FatBVHPrint.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/core/builders/refitter.cpp)


message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Debug is ON!")
  add_definitions(-D_DEBUG)
  #uncomment this to detect broken memory problems via gcc sanitizers
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize-address-use-after-scope -fno-omit-frame-pointer -fsanitize=leak -fsanitize=undefined -fsanitize=bounds-strict")
else()
  add_definitions(-DNDEBUG)
endif()

#add_library(HeavyRT)
target_compile_features(HeavyRT PUBLIC cxx_std_17)
target_compile_features(render_app PUBLIC cxx_std_17)

target_sources(HeavyRT PRIVATE
           ${CORE_SRC}
           ${EXAMPLES_SRC}
           ${BUILDERS_SRC}
           )

target_sources(render_app PRIVATE main.cpp
                ${SCENE_LOADER_SRC}
)
target_link_libraries(HeavyRT PUBLIC LiteMath)
target_link_libraries(render_app PRIVATE HeavyRT)



if(CMAKE_SYSTEM_NAME STREQUAL Windows)
    set_target_properties(HeavyRT PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
    set_target_properties(render_app PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

    target_link_libraries(HeavyRT PRIVATE project_warnings ${RAYTRACING_EMBREE_LIBS})
    target_link_libraries(render_app PRIVATE project_warnings ${RAYTRACING_EMBREE_LIBS})

    add_custom_command(TARGET render_app POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory "${PROJECT_SOURCE_DIR}/external/embree/bin_win64" $<TARGET_FILE_DIR:render_app>)
else()
    target_link_libraries(render_app PRIVATE project_warnings Threads::Threads dl ${RAYTRACING_EMBREE_LIBS})
    target_link_libraries(HeavyRT PRIVATE project_warnings Threads::Threads dl ${RAYTRACING_EMBREE_LIBS})
endif()

if(OpenMP_CXX_FOUND)
  target_link_libraries(HeavyRT PUBLIC OpenMP::OpenMP_CXX)
  target_link_libraries(render_app PUBLIC OpenMP::OpenMP_CXX) ##
endif()

