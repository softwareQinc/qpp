# Python development
target_include_directories(
  libqpp INTERFACE $<BUILD_INTERFACE:${Python3_INCLUDE_DIRS}>
                   $<INSTALL_INTERFACE:include/>)
# pybind11 library
target_include_directories(
  libqpp INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libs/>
                   $<INSTALL_INTERFACE:include/>)

# pyqpp
target_include_directories(
  libqpp INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pyqpp/include/>
                   $<INSTALL_INTERFACE:include/>)
