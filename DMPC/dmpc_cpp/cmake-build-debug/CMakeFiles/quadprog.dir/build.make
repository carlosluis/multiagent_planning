# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

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
CMAKE_COMMAND = /home/carlos/clion-2017.3.3/bin/cmake/bin/cmake

# The command to remove a file.
RM = /home/carlos/clion-2017.3.3/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/quadprog.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/quadprog.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/quadprog.dir/flags.make

CMakeFiles/quadprog.dir/dmpc.cpp.o: CMakeFiles/quadprog.dir/flags.make
CMakeFiles/quadprog.dir/dmpc.cpp.o: ../dmpc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/quadprog.dir/dmpc.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/quadprog.dir/dmpc.cpp.o -c "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/dmpc.cpp"

CMakeFiles/quadprog.dir/dmpc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/quadprog.dir/dmpc.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/dmpc.cpp" > CMakeFiles/quadprog.dir/dmpc.cpp.i

CMakeFiles/quadprog.dir/dmpc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/quadprog.dir/dmpc.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/dmpc.cpp" -o CMakeFiles/quadprog.dir/dmpc.cpp.s

CMakeFiles/quadprog.dir/dmpc.cpp.o.requires:

.PHONY : CMakeFiles/quadprog.dir/dmpc.cpp.o.requires

CMakeFiles/quadprog.dir/dmpc.cpp.o.provides: CMakeFiles/quadprog.dir/dmpc.cpp.o.requires
	$(MAKE) -f CMakeFiles/quadprog.dir/build.make CMakeFiles/quadprog.dir/dmpc.cpp.o.provides.build
.PHONY : CMakeFiles/quadprog.dir/dmpc.cpp.o.provides

CMakeFiles/quadprog.dir/dmpc.cpp.o.provides.build: CMakeFiles/quadprog.dir/dmpc.cpp.o


# Object files for target quadprog
quadprog_OBJECTS = \
"CMakeFiles/quadprog.dir/dmpc.cpp.o"

# External object files for target quadprog
quadprog_EXTERNAL_OBJECTS =

libquadprog.a: CMakeFiles/quadprog.dir/dmpc.cpp.o
libquadprog.a: CMakeFiles/quadprog.dir/build.make
libquadprog.a: CMakeFiles/quadprog.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libquadprog.a"
	$(CMAKE_COMMAND) -P CMakeFiles/quadprog.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/quadprog.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/quadprog.dir/build: libquadprog.a

.PHONY : CMakeFiles/quadprog.dir/build

CMakeFiles/quadprog.dir/requires: CMakeFiles/quadprog.dir/dmpc.cpp.o.requires

.PHONY : CMakeFiles/quadprog.dir/requires

CMakeFiles/quadprog.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/quadprog.dir/cmake_clean.cmake
.PHONY : CMakeFiles/quadprog.dir/clean

CMakeFiles/quadprog.dir/depend:
	cd "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/cmake-build-debug" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp" "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp" "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/cmake-build-debug" "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/cmake-build-debug" "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/cmake-build-debug/CMakeFiles/quadprog.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/quadprog.dir/depend

