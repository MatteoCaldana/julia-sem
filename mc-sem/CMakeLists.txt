cmake_minimum_required(VERSION 3.12.0 FATAL_ERROR)
project(mc-sem)

set(CMAKE_CXX_STANDARD "17")
set(CMAKE_CXX_STANDARD_REQUIRED "ON")

# Set default build type to Release.
if (NOT CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  add_compile_options(-O2 -DNDEBUG) # -fopt-info-vec-optimized -unroll-loops -omit-frame-pointer -Winline -march=native -mno-vzeroupper -mavx)
endif()
message(STATUS)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS)

string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPERCASE)
add_definitions("-DBUILD_TYPE_${CMAKE_BUILD_TYPE_UPPERCASE}")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfloat-conversion -Wmissing-braces -Wnon-virtual-dtor -Wall -Wextra")

file(GLOB_RECURSE TARGET_H "${CMAKE_SOURCE_DIR}/src/*.h")
file(GLOB_RECURSE TARGET_CPP "${CMAKE_SOURCE_DIR}/src/*.cpp")
set(TARGET_SRC ${TARGET_CPP} ${TARGET_H})
add_executable(${CMAKE_PROJECT_NAME} ${TARGET_SRC})