# CMakeLists
cmake_minimum_required(VERSION 3.10)
project(MyBtraj VERSION 1.0 LANGUAGES C CXX)
set (CMAKE_CXX_STANDARD 14)

find_package(OsqpEigen)
find_package(Eigen3)
find_package(OpenCV REQUIRED )

set(CMAKE_BUILD_TYPE DEBUG)
#find_package(osqp REQUIRED)
#get_target_property(osqp_INCLUDE_DIR osqp::osqpstatic INTERFACE_INCLUDE_DIRECTORIES)
#include_directories(${osqp_INCLUDE_DIR})              
##src/grid_map_2d.cc 
##src/signed_distance_field_2d.cc

set(SOURCES  src/timeallocation.cpp  src/A_star.cpp src/corridor.cpp
             src/trajectory_generator.cpp src/bezier_base.cpp 
             #src/grid_map_2d.cpp 
            #src/signed_distance_field_2d.cpp
  src/main.cpp)

include_directories (${CMAKE_SOURCE_DIR}/include)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
include_directories ("/usr/include/opencv4/")
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
add_executable(mybtraj ${SOURCES})
target_link_libraries(mybtraj OsqpEigen::OsqpEigen)
##target_link_libraries(mybtraj PUBLIC ${OpenCV_LIBS})
