cmake_minimum_required(VERSION 3.16)
project(discrete_hj_refactored CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -fno-math-errno")

find_package(Threads REQUIRED)

add_executable(dhj src/main.cpp)
target_include_directories(dhj PRIVATE include)
target_link_libraries(dhj PRIVATE Threads::Threads z)
target_compile_options(dhj PRIVATE -Wall -Wextra)

add_executable(dhj_compare_checkpoints tools/compare_checkpoints.cpp)
target_include_directories(dhj_compare_checkpoints PRIVATE include)
