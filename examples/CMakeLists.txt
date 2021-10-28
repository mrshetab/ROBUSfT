add_executable(Example1_RegularTemplate Example1_RegularTemplate.cpp )

target_compile_options(Example1_RegularTemplate PRIVATE -std=gnu++17)
target_link_libraries(Example1_RegularTemplate ${OpenCV_LIBS} Eigen3::Eigen PopSift::popsift) 
target_link_libraries(Example1_RegularTemplate ROBUSfT) 


add_executable(Example2_IrregularTemplate Example2_IrregularTemplate.cpp )

target_compile_options(Example2_IrregularTemplate PRIVATE -std=gnu++17)
target_link_libraries(Example2_IrregularTemplate ${OpenCV_LIBS} Eigen3::Eigen PopSift::popsift) 
target_link_libraries(Example2_IrregularTemplate ROBUSfT) 





