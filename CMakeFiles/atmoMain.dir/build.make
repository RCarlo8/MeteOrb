# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.25.1/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.25.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/carlocolli/Desktop/projects/Tesi

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/carlocolli/Desktop/projects/Tesi

# Include any dependencies generated for this target.
include CMakeFiles/atmoMain.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/atmoMain.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/atmoMain.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/atmoMain.dir/flags.make

CMakeFiles/atmoMain.dir/scripts/main.cpp.o: CMakeFiles/atmoMain.dir/flags.make
CMakeFiles/atmoMain.dir/scripts/main.cpp.o: scripts/main.cpp
CMakeFiles/atmoMain.dir/scripts/main.cpp.o: CMakeFiles/atmoMain.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/carlocolli/Desktop/projects/Tesi/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/atmoMain.dir/scripts/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/atmoMain.dir/scripts/main.cpp.o -MF CMakeFiles/atmoMain.dir/scripts/main.cpp.o.d -o CMakeFiles/atmoMain.dir/scripts/main.cpp.o -c /Users/carlocolli/Desktop/projects/Tesi/scripts/main.cpp

CMakeFiles/atmoMain.dir/scripts/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/atmoMain.dir/scripts/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/carlocolli/Desktop/projects/Tesi/scripts/main.cpp > CMakeFiles/atmoMain.dir/scripts/main.cpp.i

CMakeFiles/atmoMain.dir/scripts/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/atmoMain.dir/scripts/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/carlocolli/Desktop/projects/Tesi/scripts/main.cpp -o CMakeFiles/atmoMain.dir/scripts/main.cpp.s

CMakeFiles/atmoMain.dir/src/GNSSlibrary/utility.cpp.o: CMakeFiles/atmoMain.dir/flags.make
CMakeFiles/atmoMain.dir/src/GNSSlibrary/utility.cpp.o: src/GNSSlibrary/utility.cpp
CMakeFiles/atmoMain.dir/src/GNSSlibrary/utility.cpp.o: CMakeFiles/atmoMain.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/carlocolli/Desktop/projects/Tesi/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/atmoMain.dir/src/GNSSlibrary/utility.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/atmoMain.dir/src/GNSSlibrary/utility.cpp.o -MF CMakeFiles/atmoMain.dir/src/GNSSlibrary/utility.cpp.o.d -o CMakeFiles/atmoMain.dir/src/GNSSlibrary/utility.cpp.o -c /Users/carlocolli/Desktop/projects/Tesi/src/GNSSlibrary/utility.cpp

CMakeFiles/atmoMain.dir/src/GNSSlibrary/utility.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/atmoMain.dir/src/GNSSlibrary/utility.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/carlocolli/Desktop/projects/Tesi/src/GNSSlibrary/utility.cpp > CMakeFiles/atmoMain.dir/src/GNSSlibrary/utility.cpp.i

CMakeFiles/atmoMain.dir/src/GNSSlibrary/utility.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/atmoMain.dir/src/GNSSlibrary/utility.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/carlocolli/Desktop/projects/Tesi/src/GNSSlibrary/utility.cpp -o CMakeFiles/atmoMain.dir/src/GNSSlibrary/utility.cpp.s

CMakeFiles/atmoMain.dir/src/GNSSlibrary/atmo.cpp.o: CMakeFiles/atmoMain.dir/flags.make
CMakeFiles/atmoMain.dir/src/GNSSlibrary/atmo.cpp.o: src/GNSSlibrary/atmo.cpp
CMakeFiles/atmoMain.dir/src/GNSSlibrary/atmo.cpp.o: CMakeFiles/atmoMain.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/carlocolli/Desktop/projects/Tesi/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/atmoMain.dir/src/GNSSlibrary/atmo.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/atmoMain.dir/src/GNSSlibrary/atmo.cpp.o -MF CMakeFiles/atmoMain.dir/src/GNSSlibrary/atmo.cpp.o.d -o CMakeFiles/atmoMain.dir/src/GNSSlibrary/atmo.cpp.o -c /Users/carlocolli/Desktop/projects/Tesi/src/GNSSlibrary/atmo.cpp

CMakeFiles/atmoMain.dir/src/GNSSlibrary/atmo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/atmoMain.dir/src/GNSSlibrary/atmo.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/carlocolli/Desktop/projects/Tesi/src/GNSSlibrary/atmo.cpp > CMakeFiles/atmoMain.dir/src/GNSSlibrary/atmo.cpp.i

CMakeFiles/atmoMain.dir/src/GNSSlibrary/atmo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/atmoMain.dir/src/GNSSlibrary/atmo.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/carlocolli/Desktop/projects/Tesi/src/GNSSlibrary/atmo.cpp -o CMakeFiles/atmoMain.dir/src/GNSSlibrary/atmo.cpp.s

# Object files for target atmoMain
atmoMain_OBJECTS = \
"CMakeFiles/atmoMain.dir/scripts/main.cpp.o" \
"CMakeFiles/atmoMain.dir/src/GNSSlibrary/utility.cpp.o" \
"CMakeFiles/atmoMain.dir/src/GNSSlibrary/atmo.cpp.o"

# External object files for target atmoMain
atmoMain_EXTERNAL_OBJECTS =

atmoMain: CMakeFiles/atmoMain.dir/scripts/main.cpp.o
atmoMain: CMakeFiles/atmoMain.dir/src/GNSSlibrary/utility.cpp.o
atmoMain: CMakeFiles/atmoMain.dir/src/GNSSlibrary/atmo.cpp.o
atmoMain: CMakeFiles/atmoMain.dir/build.make
atmoMain: libGNSSlibrary.a
atmoMain: CMakeFiles/atmoMain.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/carlocolli/Desktop/projects/Tesi/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable atmoMain"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/atmoMain.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/atmoMain.dir/build: atmoMain
.PHONY : CMakeFiles/atmoMain.dir/build

CMakeFiles/atmoMain.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/atmoMain.dir/cmake_clean.cmake
.PHONY : CMakeFiles/atmoMain.dir/clean

CMakeFiles/atmoMain.dir/depend:
	cd /Users/carlocolli/Desktop/projects/Tesi && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/carlocolli/Desktop/projects/Tesi /Users/carlocolli/Desktop/projects/Tesi /Users/carlocolli/Desktop/projects/Tesi /Users/carlocolli/Desktop/projects/Tesi /Users/carlocolli/Desktop/projects/Tesi/CMakeFiles/atmoMain.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/atmoMain.dir/depend

