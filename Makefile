# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.31

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

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
CMAKE_COMMAND = /opt/homebrew/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/xiongzheng/software/B4

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/xiongzheng/software/B4/B4e

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Running CMake cache editor..."
	/opt/homebrew/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --cyan "Running CMake to regenerate build system..."
	/opt/homebrew/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /Users/xiongzheng/software/B4/B4e/CMakeFiles /Users/xiongzheng/software/B4/B4e//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /Users/xiongzheng/software/B4/B4e/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named B4

# Build rule for target.
B4: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 B4
.PHONY : B4

# fast build rule for target.
B4/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/B4.dir/build.make CMakeFiles/B4.dir/build
.PHONY : B4/fast

#=============================================================================
# Target rules for targets named exampleB4a

# Build rule for target.
exampleB4a: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 exampleB4a
.PHONY : exampleB4a

# fast build rule for target.
exampleB4a/fast:
	$(MAKE) $(MAKESILENT) -f B4a/CMakeFiles/exampleB4a.dir/build.make B4a/CMakeFiles/exampleB4a.dir/build
.PHONY : exampleB4a/fast

#=============================================================================
# Target rules for targets named exampleB4b

# Build rule for target.
exampleB4b: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 exampleB4b
.PHONY : exampleB4b

# fast build rule for target.
exampleB4b/fast:
	$(MAKE) $(MAKESILENT) -f B4b/CMakeFiles/exampleB4b.dir/build.make B4b/CMakeFiles/exampleB4b.dir/build
.PHONY : exampleB4b/fast

#=============================================================================
# Target rules for targets named exampleB4e

# Build rule for target.
exampleB4e: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 exampleB4e
.PHONY : exampleB4e

# fast build rule for target.
exampleB4e/fast:
	$(MAKE) $(MAKESILENT) -f B4e/CMakeFiles/exampleB4e.dir/build.make B4e/CMakeFiles/exampleB4e.dir/build
.PHONY : exampleB4e/fast

#=============================================================================
# Target rules for targets named exampleB4d

# Build rule for target.
exampleB4d: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 exampleB4d
.PHONY : exampleB4d

# fast build rule for target.
exampleB4d/fast:
	$(MAKE) $(MAKESILENT) -f B4d/CMakeFiles/exampleB4d.dir/build.make B4d/CMakeFiles/exampleB4d.dir/build
.PHONY : exampleB4d/fast

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... B4"
	@echo "... exampleB4a"
	@echo "... exampleB4b"
	@echo "... exampleB4e"
	@echo "... exampleB4d"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

