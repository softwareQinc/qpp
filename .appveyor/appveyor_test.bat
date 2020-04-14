@echo off
IF %COMPILER%==msvc2019 (
    %APPVEYOR_BUILD_FOLDER%\unit_tests\build\tests\Debug\qpp_testing.exe --gtest_filter=-qpp_Timer*
)
IF %COMPILER%==msys2 (
    %APPVEYOR_BUILD_FOLDER%\unit_tests\build\tests\qpp_testing.exe --gtest_filter=-qpp_Timer*
)
