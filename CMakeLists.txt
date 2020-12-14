project(TubeSystem)

cmake_minimum_required(VERSION 3.0)

add_compile_options(-fno-omit-frame-pointer) 
#add_compile_options(-fsanitize=address) 
#add_compile_options(-lasan)

# What we want: LDFLAGS="-fsanitize=address"
#link_libraries("-fsanitize=address") 

# set output paths
set(PROJECT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

if(WIN32)
  MESSAGE("CMAKE_PREFIX_PATH: " $CMAKE_PREFIX_PATH)
  set(CMAKE_LIBRARY_PATH ${CMAKE_SOURCE_DIR}/src/ ${CMAKE_SOURCE_DIR}/lib/win7)
endif()

# default to Release builds
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set (CMAKE_CXX_STANDARD 11)

# attempt to find the system's GLEW; build the included one if unsuccessful
find_package(GLEW QUIET)
if (GLEW_FOUND)
  MESSAGE("Using system GLEW")
  set(GLEW_INCLUDE_DIRS ${GLEW_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS}/GL) # At least on Ubuntu, system glew.h is actually in GL/
  link_libraries(${GLEW_LIBRARIES})
else()
  MESSAGE("System GLEW not found... falling back to local GLEW")
  add_subdirectory(lib/glew)
  set(GLEW_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/lib/glew)
  set(GLEW_LIBRARIES glew)
endif()

# Attempt to find the system's GLFW; build the included one if unsuccessful
find_package(glfw3 QUIET)
if (glfw3_FOUND)
  MESSAGE("Using system GLFW")
  # Note: target_link_libraries(glfw) performed in src/CMakeLists.txt should
  # actually bring in the necessary header files.
else()
  MESSAGE("System GLFW not found... falling back to local GLFW")
  include_directories(${CMAKE_SOURCE_DIR}/lib/glfw/include/)
  add_subdirectory(lib/glfw)
endif()

# add source directory to include path
include_directories(${CMAKE_SOURCE_DIR}/src/)
include_directories(${CMAKE_SOURCE_DIR}/lib/lodePNG/)
#include_directories(${CMAKE_SOURCE_DIR}/lib/eigen-3.3.7)


# build source directory
add_subdirectory(lib/lodePNG)
#add_subdirectory(lib/eigen-3.3.7)
add_subdirectory(src)

# documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
	# prepare doxygen configuration file
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
	# add doc as target
	add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
endif(DOXYGEN_FOUND)


# google test

# Download and unpack googletest at configure time
#configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
#execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
#  RESULT_VARIABLE result
 # WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
#if(result)
 # message(FATAL_ERROR "CMake step for googletest failed: ${result}")
#endif()
#execute_process(COMMAND ${CMAKE_COMMAND} --build .
 # RESULT_VARIABLE result
 # WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
#if(result)
 # message(FATAL_ERROR "Build step for googletest failed: ${result}")
#endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)



# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

# Profiler gprof
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

