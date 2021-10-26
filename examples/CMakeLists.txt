add_executable(CPU_GPU_RegularTemplate CPU_GPU_RegularTemplate.cpp )

target_compile_options(CPU_GPU_RegularTemplate PRIVATE -std=gnu++17)
target_link_libraries(CPU_GPU_RegularTemplate ${OpenCV_LIBS} Eigen3::Eigen PopSift::popsift) 
target_link_libraries(CPU_GPU_RegularTemplate ROBUSfT) 


add_executable(CPU_GPU_IrregularTemplate CPU_GPU_IrregularTemplate.cpp )

target_compile_options(CPU_GPU_IrregularTemplate PRIVATE -std=gnu++17)
target_link_libraries(CPU_GPU_IrregularTemplate ${OpenCV_LIBS} Eigen3::Eigen PopSift::popsift) 
target_link_libraries(CPU_GPU_IrregularTemplate ROBUSfT) 





