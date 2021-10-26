cmake_minimum_required(VERSION 3.10.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  message(STATUS "Build type not set, building in Release configuration")
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -fopenmp") 

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CUDA_TOOLKIT_ROOT_DIR "/usr/local/cuda-11.3") 

project(ROBUSfTt VERSION 1.0.0)

find_package(OpenCV REQUIRED)
find_package (Eigen3 REQUIRED)
find_package(PopSift CONFIG REQUIRED)

if(PopSift_USE_NVTX_PROFILING)
  list(APPEND PD_LINK_LIBS ${CUDA_NVTX_LIBRARY})
endif(PopSift_USE_NVTX_PROFILING)

include_directories(${OpenCV_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})
add_definitions( ${EIGEN3_DEFINITIONS} )

add_definitions(-DROOT_FOLDER="${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(ROBUSfT)
add_subdirectory(examples)






