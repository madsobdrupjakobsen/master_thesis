# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/switching_times_2nd.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/switching_times_2nd.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/switching_times_2nd.dir/flags.make

CMakeFiles/switching_times_2nd.dir/main_2nd_order.cpp.o: CMakeFiles/switching_times_2nd.dir/flags.make
CMakeFiles/switching_times_2nd.dir/main_2nd_order.cpp.o: ../main_2nd_order.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/switching_times_2nd.dir/main_2nd_order.cpp.o"
	/usr/local/bin/g++-9  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/switching_times_2nd.dir/main_2nd_order.cpp.o -c /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/main_2nd_order.cpp

CMakeFiles/switching_times_2nd.dir/main_2nd_order.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/switching_times_2nd.dir/main_2nd_order.cpp.i"
	/usr/local/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/main_2nd_order.cpp > CMakeFiles/switching_times_2nd.dir/main_2nd_order.cpp.i

CMakeFiles/switching_times_2nd.dir/main_2nd_order.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/switching_times_2nd.dir/main_2nd_order.cpp.s"
	/usr/local/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/main_2nd_order.cpp -o CMakeFiles/switching_times_2nd.dir/main_2nd_order.cpp.s

CMakeFiles/switching_times_2nd.dir/src/switching-times.cpp.o: CMakeFiles/switching_times_2nd.dir/flags.make
CMakeFiles/switching_times_2nd.dir/src/switching-times.cpp.o: ../src/switching-times.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/switching_times_2nd.dir/src/switching-times.cpp.o"
	/usr/local/bin/g++-9  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/switching_times_2nd.dir/src/switching-times.cpp.o -c /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/src/switching-times.cpp

CMakeFiles/switching_times_2nd.dir/src/switching-times.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/switching_times_2nd.dir/src/switching-times.cpp.i"
	/usr/local/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/src/switching-times.cpp > CMakeFiles/switching_times_2nd.dir/src/switching-times.cpp.i

CMakeFiles/switching_times_2nd.dir/src/switching-times.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/switching_times_2nd.dir/src/switching-times.cpp.s"
	/usr/local/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/src/switching-times.cpp -o CMakeFiles/switching_times_2nd.dir/src/switching-times.cpp.s

CMakeFiles/switching_times_2nd.dir/src/switching-times-ice-tank-2ndorder.cpp.o: CMakeFiles/switching_times_2nd.dir/flags.make
CMakeFiles/switching_times_2nd.dir/src/switching-times-ice-tank-2ndorder.cpp.o: ../src/switching-times-ice-tank-2ndorder.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/switching_times_2nd.dir/src/switching-times-ice-tank-2ndorder.cpp.o"
	/usr/local/bin/g++-9  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/switching_times_2nd.dir/src/switching-times-ice-tank-2ndorder.cpp.o -c /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/src/switching-times-ice-tank-2ndorder.cpp

CMakeFiles/switching_times_2nd.dir/src/switching-times-ice-tank-2ndorder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/switching_times_2nd.dir/src/switching-times-ice-tank-2ndorder.cpp.i"
	/usr/local/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/src/switching-times-ice-tank-2ndorder.cpp > CMakeFiles/switching_times_2nd.dir/src/switching-times-ice-tank-2ndorder.cpp.i

CMakeFiles/switching_times_2nd.dir/src/switching-times-ice-tank-2ndorder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/switching_times_2nd.dir/src/switching-times-ice-tank-2ndorder.cpp.s"
	/usr/local/bin/g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/src/switching-times-ice-tank-2ndorder.cpp -o CMakeFiles/switching_times_2nd.dir/src/switching-times-ice-tank-2ndorder.cpp.s

# Object files for target switching_times_2nd
switching_times_2nd_OBJECTS = \
"CMakeFiles/switching_times_2nd.dir/main_2nd_order.cpp.o" \
"CMakeFiles/switching_times_2nd.dir/src/switching-times.cpp.o" \
"CMakeFiles/switching_times_2nd.dir/src/switching-times-ice-tank-2ndorder.cpp.o"

# External object files for target switching_times_2nd
switching_times_2nd_EXTERNAL_OBJECTS =

switching_times_2nd.cpython-37m-darwin.so: CMakeFiles/switching_times_2nd.dir/main_2nd_order.cpp.o
switching_times_2nd.cpython-37m-darwin.so: CMakeFiles/switching_times_2nd.dir/src/switching-times.cpp.o
switching_times_2nd.cpython-37m-darwin.so: CMakeFiles/switching_times_2nd.dir/src/switching-times-ice-tank-2ndorder.cpp.o
switching_times_2nd.cpython-37m-darwin.so: CMakeFiles/switching_times_2nd.dir/build.make
switching_times_2nd.cpython-37m-darwin.so: CMakeFiles/switching_times_2nd.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX shared module switching_times_2nd.cpython-37m-darwin.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/switching_times_2nd.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/switching_times_2nd.dir/build: switching_times_2nd.cpython-37m-darwin.so

.PHONY : CMakeFiles/switching_times_2nd.dir/build

CMakeFiles/switching_times_2nd.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/switching_times_2nd.dir/cmake_clean.cmake
.PHONY : CMakeFiles/switching_times_2nd.dir/clean

CMakeFiles/switching_times_2nd.dir/depend:
	cd /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/cmake-build-debug /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/cmake-build-debug /Users/madsobdrup/Dropbox/Skole/DTU/Studie/MASTER/CODE/switching_time_optimization/control/switching-times/cmake-build-debug/CMakeFiles/switching_times_2nd.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/switching_times_2nd.dir/depend
