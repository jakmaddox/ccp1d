cmake_minimum_required(VERSION 3.16)

project(ccp1d
    VERSION 0.1.0
    DESCRIPTION "1D drift-diffusion model of a capacitively-coupled RF discharge"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "build type" FORCE)
endif()

add_library(ccp1d_core STATIC
    src/discharge.cpp)

target_include_directories(ccp1d_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_compile_options(ccp1d_core PRIVATE
    -Wall -Wextra -Wpedantic -Wshadow -Wconversion)

add_executable(ccp1d src/main.cpp)
target_link_libraries(ccp1d PRIVATE ccp1d_core)

enable_testing()

add_executable(test_numerics tests/test_numerics.cpp)
target_link_libraries(test_numerics PRIVATE ccp1d_core)

add_test(NAME numerics COMMAND test_numerics)
