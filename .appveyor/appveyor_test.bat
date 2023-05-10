@echo off
IF %COMPILER%==msvc2019 (
    ctest --test-dir build -E qpp_Timer
)
IF %COMPILER%==msys2 (
    ctest --test-dir build -E qpp_Timer
)
