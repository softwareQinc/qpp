@echo off
IF %COMPILER%==msvc2019 (
    ctest -E qpp_Timer
)
IF %COMPILER%==msys2 (
    ctest -E qpp_Timer
)
