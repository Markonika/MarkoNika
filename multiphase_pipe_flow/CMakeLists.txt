cmake_minimum_required(VERSION 3.15)
project(multiphase_pipe_flow CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra)

# Optional OpenMP support (SolverOptions has no dependency on this --
# every parallel region is guarded by #ifdef _OPENMP, so the solver and
# its drivers build and behave identically with or without it). Used two
# ways: (1) validate_shoham/validate_quantitative parallelize their outer
# per-case loop, since each case is an independent FourFieldSolver
# instance -- embarrassingly parallel, and where the real speedup is for
# the batch sizes (182-349 cases) this codebase's own validation sweeps
# actually run; (2) FourFieldSolver's own per-cell/per-face loops
# (updateLayerMomentum, updateContinuity, computeClosures) are also
# annotated, for users running a single large-N simulation, though see
# README.md, "Parallelization" for honest scaling results at this
# project's typical N (60-300 cells) vs. larger N.
find_package(OpenMP)

add_library(mfs
    src/PipeGeometry.cpp
    src/Closures.cpp
    src/FourFieldSolver.cpp
    src/FlowRegimeClassifier.cpp
)
target_include_directories(mfs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(OpenMP_CXX_FOUND)
    target_link_libraries(mfs PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(mfs_demo src/main.cpp)
target_link_libraries(mfs_demo PRIVATE mfs)

# Validation driver comparing predicted flow regimes against the Shoham
# (1982) experimental dataset -- see VALIDATION.md.
add_executable(validate_shoham validation/validate_shoham.cpp)
target_link_libraries(validate_shoham PRIVATE mfs)
if(OpenMP_CXX_FOUND)
    target_link_libraries(validate_shoham PRIVATE OpenMP::OpenMP_CXX)
endif()

# Quantitative validation driver: predicted liquid holdup and frictional
# pressure gradient vs. Kokal (1987), Newton (1997), and a multi-source
# compiled holdup database -- see the paper's quantitative validation
# section for methodology and results.
add_executable(validate_quantitative validation/validate_quantitative.cpp)
target_link_libraries(validate_quantitative PRIVATE mfs)
if(OpenMP_CXX_FOUND)
    target_link_libraries(validate_quantitative PRIVATE OpenMP::OpenMP_CXX)
endif()

enable_testing()
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt)
    add_subdirectory(tests)
endif()
