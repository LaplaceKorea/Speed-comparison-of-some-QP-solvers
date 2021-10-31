cmake_minimum_required(VERSION 3.5)
project(basic_utils LANGUAGES CXX)
add_definitions(-std=c++14 -w)

SET(CMAKE_BUILD_TYPE "Release")  
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -g2 -ggdb")  
# SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

find_package(Eigen3 REQUIRED)

include_directories ( ${EIGEN3_INCLUDE_DIRS} )
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/MPC/qpOASES/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

FILE (GLOB SRC MPC/qpOASES/src/*.cpp)
add_library(qpOASES ${SRC})

find_package(blasfeo REQUIRED)
find_package(hpipm REQUIRED)
find_package(casadi REQUIRED)
find_package(OsqpEigen REQUIRED)
find_package(acados REQUIRED)

if(blasfeo_FOUND)
    message(STATUS "Found BLASFEO")
else()
    message(WARNING "Could not find BLASFEO")
endif()
if(hpipm_FOUND)
    message(STATUS "Found HPIPM")
else()
    message(WARNING "Could not find HPIPM")
endif()

include_directories(${blasfeo_INCLUDE_DIRS})
include_directories(${hpipm_INCLUDE_DIRS})
list(APPEND HPIPM_LIBS hpipm blasfeo)
add_library(mpc_control MPC/MpcControl.cpp)
target_link_libraries(mpc_control qpOASES ${HPIPM_LIBS} casadi OsqpEigen::OsqpEigen acados)

add_executable(mpc_test examples/mpcTest.cpp)
target_link_libraries(mpc_test mpc_control)

add_executable(hpipm_test examples/hpipmTest.cpp)
target_link_libraries(hpipm_test mpc_control)

add_executable(acados_test examples/acadosTest.cpp)
target_link_libraries(acados_test mpc_control)
