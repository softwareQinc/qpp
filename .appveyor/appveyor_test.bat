@echo off
IF %COMPILER%==msvc2019 (
    cd build
    ctest -E qpp_Timer
)
IF %COMPILER%==msys2 (
    cd build
    ctest -E qpp_Timer
)
