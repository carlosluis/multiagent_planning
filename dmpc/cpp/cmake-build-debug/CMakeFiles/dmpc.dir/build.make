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
include CMakeFiles/dmpc.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/dmpc.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/dmpc.dir/flags.make

CMakeFiles/dmpc.dir/dmpc.cpp.o: CMakeFiles/dmpc.dir/flags.make
CMakeFiles/dmpc.dir/dmpc.cpp.o: ../dmpc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/dmpc.dir/dmpc.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dmpc.dir/dmpc.cpp.o -c "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/dmpc.cpp"

CMakeFiles/dmpc.dir/dmpc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dmpc.dir/dmpc.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/dmpc.cpp" > CMakeFiles/dmpc.dir/dmpc.cpp.i

CMakeFiles/dmpc.dir/dmpc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dmpc.dir/dmpc.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/dmpc.cpp" -o CMakeFiles/dmpc.dir/dmpc.cpp.s

CMakeFiles/dmpc.dir/dmpc.cpp.o.requires:

.PHONY : CMakeFiles/dmpc.dir/dmpc.cpp.o.requires

CMakeFiles/dmpc.dir/dmpc.cpp.o.provides: CMakeFiles/dmpc.dir/dmpc.cpp.o.requires
	$(MAKE) -f CMakeFiles/dmpc.dir/build.make CMakeFiles/dmpc.dir/dmpc.cpp.o.provides.build
.PHONY : CMakeFiles/dmpc.dir/dmpc.cpp.o.provides

CMakeFiles/dmpc.dir/dmpc.cpp.o.provides.build: CMakeFiles/dmpc.dir/dmpc.cpp.o


# Object files for target dmpc
dmpc_OBJECTS = \
"CMakeFiles/dmpc.dir/dmpc.cpp.o"

# External object files for target dmpc
dmpc_EXTERNAL_OBJECTS =

libdmpc.a: CMakeFiles/dmpc.dir/dmpc.cpp.o
libdmpc.a: CMakeFiles/dmpc.dir/build.make
libdmpc.a: CMakeFiles/dmpc.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/cmake-build-debug/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libdmpc.a"
	$(CMAKE_COMMAND) -P CMakeFiles/dmpc.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dmpc.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/dmpc.dir/build: libdmpc.a

.PHONY : CMakeFiles/dmpc.dir/build

CMakeFiles/dmpc.dir/requires: CMakeFiles/dmpc.dir/dmpc.cpp.o.requires

.PHONY : CMakeFiles/dmpc.dir/requires

CMakeFiles/dmpc.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/dmpc.dir/cmake_clean.cmake
.PHONY : CMakeFiles/dmpc.dir/clean

CMakeFiles/dmpc.dir/depend:
	cd "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/cmake-build-debug" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp" "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp" "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/cmake-build-debug" "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/cmake-build-debug" "/home/carlos/Documents/UTIAS/First Year/Winter 2018/ECE1505/Project/dec_SQP/DMPC/dmpc_cpp/cmake-build-debug/CMakeFiles/dmpc.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/dmpc.dir/depend

