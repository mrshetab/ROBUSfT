# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/reza/Codes/Git/ROBUSfT

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/reza/Codes/Git/ROBUSfT/build

# Include any dependencies generated for this target.
include examples/CMakeFiles/Example2_IrregularTemplate.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/Example2_IrregularTemplate.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/Example2_IrregularTemplate.dir/flags.make

examples/CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.o: examples/CMakeFiles/Example2_IrregularTemplate.dir/flags.make
examples/CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.o: ../examples/Example2_IrregularTemplate.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/reza/Codes/Git/ROBUSfT/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.o"
	cd /home/reza/Codes/Git/ROBUSfT/build/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.o -c /home/reza/Codes/Git/ROBUSfT/examples/Example2_IrregularTemplate.cpp

examples/CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.i"
	cd /home/reza/Codes/Git/ROBUSfT/build/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/reza/Codes/Git/ROBUSfT/examples/Example2_IrregularTemplate.cpp > CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.i

examples/CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.s"
	cd /home/reza/Codes/Git/ROBUSfT/build/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/reza/Codes/Git/ROBUSfT/examples/Example2_IrregularTemplate.cpp -o CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.s

examples/CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.o.requires:

.PHONY : examples/CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.o.requires

examples/CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.o.provides: examples/CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.o.requires
	$(MAKE) -f examples/CMakeFiles/Example2_IrregularTemplate.dir/build.make examples/CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.o.provides.build
.PHONY : examples/CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.o.provides

examples/CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.o.provides.build: examples/CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.o


# Object files for target Example2_IrregularTemplate
Example2_IrregularTemplate_OBJECTS = \
"CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.o"

# External object files for target Example2_IrregularTemplate
Example2_IrregularTemplate_EXTERNAL_OBJECTS =

examples/Example2_IrregularTemplate: examples/CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.o
examples/Example2_IrregularTemplate: examples/CMakeFiles/Example2_IrregularTemplate.dir/build.make
examples/Example2_IrregularTemplate: ROBUSfT/libROBUSfT.so
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_gapi.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_stitching.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_aruco.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_bgsegm.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_bioinspired.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_ccalib.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_dnn_objdetect.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_dnn_superres.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_dpm.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_highgui.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_face.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_freetype.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_fuzzy.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_hfs.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_img_hash.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_intensity_transform.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_line_descriptor.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_mcc.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_quality.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_rapid.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_reg.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_rgbd.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_saliency.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_stereo.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_structured_light.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_phase_unwrapping.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_superres.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_optflow.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_surface_matching.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_tracking.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_datasets.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_plot.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_text.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_dnn.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_videostab.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_videoio.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_xfeatures2d.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_ml.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_shape.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_ximgproc.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_video.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_xobjdetect.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_imgcodecs.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_objdetect.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_calib3d.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_features2d.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_flann.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_xphoto.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_photo.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_imgproc.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libopencv_core.so.4.5.0
examples/Example2_IrregularTemplate: /usr/local/lib/libpopsift.so.1.0.0
examples/Example2_IrregularTemplate: /usr/local/cuda-11.3/lib64/libcudart.so
examples/Example2_IrregularTemplate: /usr/local/cuda-11.3/lib64/libcudadevrt.a
examples/Example2_IrregularTemplate: /usr/local/cuda-11.3/lib64/libcublas.so
examples/Example2_IrregularTemplate: examples/CMakeFiles/Example2_IrregularTemplate.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/reza/Codes/Git/ROBUSfT/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Example2_IrregularTemplate"
	cd /home/reza/Codes/Git/ROBUSfT/build/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Example2_IrregularTemplate.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/Example2_IrregularTemplate.dir/build: examples/Example2_IrregularTemplate

.PHONY : examples/CMakeFiles/Example2_IrregularTemplate.dir/build

examples/CMakeFiles/Example2_IrregularTemplate.dir/requires: examples/CMakeFiles/Example2_IrregularTemplate.dir/Example2_IrregularTemplate.cpp.o.requires

.PHONY : examples/CMakeFiles/Example2_IrregularTemplate.dir/requires

examples/CMakeFiles/Example2_IrregularTemplate.dir/clean:
	cd /home/reza/Codes/Git/ROBUSfT/build/examples && $(CMAKE_COMMAND) -P CMakeFiles/Example2_IrregularTemplate.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/Example2_IrregularTemplate.dir/clean

examples/CMakeFiles/Example2_IrregularTemplate.dir/depend:
	cd /home/reza/Codes/Git/ROBUSfT/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/reza/Codes/Git/ROBUSfT /home/reza/Codes/Git/ROBUSfT/examples /home/reza/Codes/Git/ROBUSfT/build /home/reza/Codes/Git/ROBUSfT/build/examples /home/reza/Codes/Git/ROBUSfT/build/examples/CMakeFiles/Example2_IrregularTemplate.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/Example2_IrregularTemplate.dir/depend

