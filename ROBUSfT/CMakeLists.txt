find_package(OpenCV REQUIRED)
find_package (Eigen3 REQUIRED)
find_package(PopSift CONFIG REQUIRED)

if(PopSift_USE_NVTX_PROFILING)
  list(APPEND PD_LINK_LIBS ${CUDA_NVTX_LIBRARY})
endif(PopSift_USE_NVTX_PROFILING)

include_directories(${OpenCV_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})
add_definitions( ${EIGEN3_DEFINITIONS} )

add_library(ROBUSfT SHARED
            src/ROBUSfT/ROBUSfT.cpp 
            src/ROBUSfT/ROBUSfT.h 
            src/ROBUSfT/delaunator_h.cpp 
            src/ROBUSfT/delaunator_h.hpp 
            src/ROBUSfT/functions.cpp 
            src/ROBUSfT/functions.h 
            src/ROBUSfT/bbs_eigen.cpp 
            src/ROBUSfT/bbs_eigen.h 
            src/ROBUSfT/headers.h 
            ) 

target_include_directories(ROBUSfT PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_compile_options(ROBUSfT PRIVATE -std=gnu++17)
target_link_libraries(ROBUSfT ${OpenCV_LIBS} Eigen3::Eigen PopSift::popsift) 