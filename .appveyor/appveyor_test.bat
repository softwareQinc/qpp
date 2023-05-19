@echo off
IF %COMPILER%==msvc2022 (
    ctest --test-dir build -E qpp_Timer
)
IF %COMPILER%==msys2 (
    ctest --test-dir build -E qpp_Timer
)
