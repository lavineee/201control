cmake_minimum_required (VERSION 3.8)
project(mujocoTest3D)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib)


set(TOOLS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tools/tools.cpp)
add_library(tools STATIC ${TOOLS_SOURCES})
target_include_directories(tools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../lib/include)
target_link_directories(tools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../lib/lib)
target_link_libraries(tools PUBLIC mujoco)


set(ROBBIE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/robbie/Robbie.cpp)
add_library(robbie STATIC ${ROBBIE_SOURCES})
target_include_directories(robbie PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../lib/include)
target_link_directories(robbie PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib)
target_link_libraries(robbie PUBLIC tools)