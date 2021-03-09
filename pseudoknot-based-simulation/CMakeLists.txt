# Check that the build does not happen in the source tree
IF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
   MESSAGE(FATAL_ERROR "Prevented in-tree built. This is bad practice.")
ENDIF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})

cmake_minimum_required (VERSION 2.6)
project (../)
cmake_policy(SET CMP0012 NEW)
set(CMAKE_CXX_STANDARD 17)
include_directories("src")
# Set the name of the executable to be produced
set(EXECUTABLE_NAME testsim )
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
# Add all .cpp - files under src/ directory
file(GLOB_RECURSE SOURCES LIST_DIRECTORIES true src/*.hpp src/*.cpp)

find_package(OpenMP)

add_executable(${EXECUTABLE_NAME} ${SOURCES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(${EXECUTABLE_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()
