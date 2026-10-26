cmake_minimum_required(VERSION 3.16)

project(hex_tiles VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Boost CONFIG REQUIRED)
find_package(Eigen3 3.3 REQUIRED)

add_executable(hex_tiles
    src/main.cpp
    src/geometry.cpp
    src/tile.cpp
    src/hex.cpp
)

target_link_libraries(hex_tiles
    PRIVATE
        Eigen3::Eigen
)

# The exact core and renderer have no Eigen/Boost or geometric dependencies.
add_library(triangle_grid src/triangle_grid.cpp src/combinatorial_hex.cpp src/combinatorial_svg.cpp
    src/hex_grid.cpp src/hex_grid_svg.cpp src/finite_patch.cpp
    src/collars.cpp src/collar_certificates.cpp src/collar_coverage.cpp)
target_include_directories(triangle_grid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(hex_tiles PRIVATE triangle_grid)
add_executable(collar_example examples/collars.cpp)
target_link_libraries(collar_example PRIVATE triangle_grid)

include(CTest)
if(BUILD_TESTING)
    add_executable(collar_partial_tests tests/collar_partial_tests.cpp)
    target_link_libraries(collar_partial_tests PRIVATE triangle_grid)
    add_test(NAME collar_partial_reconstruction COMMAND collar_partial_tests)
    add_executable(collars_tests tests/collars_tests.cpp)
    target_link_libraries(collars_tests PRIVATE triangle_grid)
    target_compile_definitions(collars_tests PRIVATE COLLAR_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(NAME collar_catalog_certificates COMMAND collars_tests)
    add_executable(collars_svg_tests tests/collars_svg_tests.cpp)
    target_link_libraries(collars_svg_tests PRIVATE triangle_grid)
    add_test(NAME collar_svg_generation COMMAND collars_svg_tests)
    add_executable(triangle_grid_tests tests/triangle_grid_tests.cpp)
    target_link_libraries(triangle_grid_tests PRIVATE triangle_grid)
    add_test(NAME triangle_grid_transformations COMMAND triangle_grid_tests)

    add_executable(finite_patch_tests tests/finite_patch_tests.cpp)
    target_link_libraries(finite_patch_tests PRIVATE triangle_grid)
    add_test(NAME finite_patch_substitution COMMAND finite_patch_tests)

    add_executable(hex_grid_tests tests/hex_grid_tests.cpp)
    target_link_libraries(hex_grid_tests PRIVATE triangle_grid)
    add_test(NAME hex_grid_reconstruction COMMAND hex_grid_tests)

    add_executable(combinatorial_hex_tests tests/combinatorial_hex_tests.cpp
        src/hex.cpp src/tile.cpp src/geometry.cpp)
    target_link_libraries(combinatorial_hex_tests PRIVATE triangle_grid Eigen3::Eigen)
    add_test(NAME combinatorial_hex_validation COMMAND combinatorial_hex_tests)
    add_test(NAME combinatorial_cli_validation COMMAND ${CMAKE_COMMAND}
        -DHEX_TILES=$<TARGET_FILE:hex_tiles>
        -DOUTPUT_ROOT=${CMAKE_CURRENT_BINARY_DIR}/grid-test-output
        -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_validation.cmake)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        add_test(NAME collar_svg_xml COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/collars_svg_validation.py
            ${CMAKE_CURRENT_BINARY_DIR}/collar-test-output)
        set_tests_properties(collar_svg_xml PROPERTIES DEPENDS collar_svg_generation)
        add_test(NAME combinatorial_svg_xml COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/svg_validation.py
            ${CMAKE_CURRENT_BINARY_DIR}/grid-test-output)
        set_tests_properties(combinatorial_svg_xml PROPERTIES DEPENDS combinatorial_cli_validation)
    endif()
endif()
