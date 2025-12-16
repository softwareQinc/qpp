# LSP and CMake support for pyqpp

# Select the target
include(${CMAKE_CURRENT_LIST_DIR}/qpp_select_target.cmake)
qpp_select_target(QPP_TARGET "qpp_types")

# pybind11
include_directories(SYSTEM ${PYBIND11_INCLUDE_DIRS})
target_include_directories(
  ${QPP_TARGET} INTERFACE $<BUILD_INTERFACE:${PYBIND11_INCLUDE_DIRS}/include/>)

# Python development
target_include_directories(${QPP_TARGET}
                           INTERFACE $<BUILD_INTERFACE:${Python3_INCLUDE_DIRS}>)

# pyqpp
target_include_directories(
  ${QPP_TARGET}
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pyqpp/include/>)
