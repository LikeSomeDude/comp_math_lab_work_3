cmake_minimum_required(VERSION 3.16)

project(comp_math_lab_work_3 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(main
    main.cpp
    src/scalar_newton.cpp
    src/nonlinear_system.cpp
    src/newton.cpp
    src/linear_solver.cpp
)

target_include_directories(main PRIVATE include)
