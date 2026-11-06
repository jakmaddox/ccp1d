cmake_minimum_required(VERSION 3.16)

project(ccp1d
    VERSION 0.2.0
    DESCRIPTION "1D drift-diffusion model of a capacitively-coupled RF discharge"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "build type" FORCE)
endif()

# Warning flags live on an interface target so that every executable linking
# the core gets them, not only the core library itself.
add_library(ccp1d_warnings INTERFACE)
target_compile_options(ccp1d_warnings INTERFACE
    -Wall -Wextra -Wpedantic -Wshadow -Wconversion)

add_library(ccp1d_core STATIC
    src/discharge.cpp)
target_include_directories(ccp1d_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ccp1d_core PUBLIC ccp1d_warnings)

add_executable(ccp1d src/main.cpp)
target_link_libraries(ccp1d PRIVATE ccp1d_core)

add_executable(ccp1d_audit tools/audit.cpp)
target_link_libraries(ccp1d_audit PRIVATE ccp1d_core)

enable_testing()

add_executable(test_numerics tests/test_numerics.cpp)
target_link_libraries(test_numerics PRIVATE ccp1d_core)
add_test(NAME numerics COMMAND test_numerics)

add_executable(test_conservation tests/test_conservation.cpp)
target_link_libraries(test_conservation PRIVATE ccp1d_core)
add_test(NAME conservation COMMAND test_conservation)
set_tests_properties(conservation PROPERTIES TIMEOUT 300)
