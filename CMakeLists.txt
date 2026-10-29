cmake_minimum_required(VERSION 3.16)
project(tolstack CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

add_library(stackup src/stackup.cpp)
target_include_directories(stackup PUBLIC src)

add_executable(tolstack src/main.cpp)
target_link_libraries(tolstack PRIVATE stackup)

enable_testing()
add_executable(test_stackup tests/test_stackup.cpp)
target_link_libraries(test_stackup PRIVATE stackup)
add_test(NAME stackup_tests COMMAND test_stackup)
