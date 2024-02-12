cmake_minimum_required(VERSION 3.15)
project(pyqpp)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include("../cmake/qpp_eigen3.cmake")
include("../cmake/pybind11.cmake")
