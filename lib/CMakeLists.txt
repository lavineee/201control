#cmake_minimum_required(VERSION 3.16.3)
#
#SET(LIBRARY_OUTPUT_PATH ./lib)
#
#SET(LIBUITOOLS_SRC uitools.c)
#add_library(uitools STATIC ${LIBUITOOLS_SRC})
#target_include_directories(uitools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

cmake_minimum_required (VERSION 3.8)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib)

set(UITOOLS_SOURCES uitools.c)
#file(GLOB STATIC_LIB_SOURCES ./src/*.cpp)
add_library( uitools STATIC ${UITOOLS_SOURCES})
target_include_directories( uitools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_directories(uitools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib)
target_link_libraries(uitools PUBLIC mujoco glfw)

set(LODEPNG_SOURCES lodepng.cpp)
add_library(lodepng STATIC ${LODEPNG_SOURCES})
target_include_directories(lodepng PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(GLFWEX_SOURCES glfw_dispatch.cc glfw_adapter.cc)
add_library(glfwex STATIC ${GLFWEX_SOURCES})
target_include_directories(glfwex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_directories(glfwex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib)
target_link_libraries(glfwex PUBLIC mujoco glfw)

set(PLATFORMUI_SOURCES platform_ui_adapter.cc)
add_library(platformui STATIC ${PLATFORMUI_SOURCES})
target_include_directories(platformui PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(SIMULATE_SOURCES simulate.cc)
add_library(simulate STATIC ${SIMULATE_SOURCES})
target_include_directories(simulate PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_directories(glfwex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib)
target_link_libraries(glfwex PUBLIC mujoco lodepng)