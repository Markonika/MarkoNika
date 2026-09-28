cmake_minimum_required(VERSION 3.15)
project(multiphase_pipe_flow CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra)

add_library(mfs
    src/PipeGeometry.cpp
    src/Closures.cpp
    src/FourFieldSolver.cpp
    src/FlowRegimeClassifier.cpp
)
target_include_directories(mfs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(mfs_demo src/main.cpp)
target_link_libraries(mfs_demo PRIVATE mfs)

# Validation driver comparing predicted flow regimes against the Shoham
# (1982) experimental dataset -- see VALIDATION.md.
add_executable(validate_shoham validation/validate_shoham.cpp)
target_link_libraries(validate_shoham PRIVATE mfs)

enable_testing()
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt)
    add_subdirectory(tests)
endif()
