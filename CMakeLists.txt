cmake_minimum_required(VERSION 3.20)
project(RobBIE_V2_SIMULATE)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pthread")
add_subdirectory(lib)
add_subdirectory(src)
#
find_package(OpenCV 4 REQUIRED)
find_package(Python3 COMPONENTS Development NumPy)

include_directories(lib/include)
include_directories(src)
include_directories("usr/include/eigen3")

add_executable(main main.cpp src/controller.cpp src/log/logFile.cpp src/robbie/calculateRobBIE.cpp src/parameters/parameters.cpp )
target_link_libraries(main -lmujoco -lglfw uitools lodepng glfwex platformui simulate ${OpenCV_LIBS} Python3::Python Python3::NumPy robbie tools)


add_executable(test_newEst test/test_newEst.cpp)

add_executable(test_standing test/test_standing.cpp src/parameters/parameters.cpp)
add_executable(test_k test/test_k.cpp src/parameters/parameters.cpp test/test_k.cpp)