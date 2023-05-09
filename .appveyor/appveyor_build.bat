@echo off
IF %COMPILER%==msvc2019 (
    @echo on
    CALL "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"

    cmake -B build -DWITH_EXAMPLES=ON -DWITH_UNIT_TESTS=ON -DEIGEN3_INSTALL_DIR="c:\%EIGEN3_INSTALL_PATH%"
    cmake --build build --target INSTALL

    cd %APPVEYOR_BUILD_FOLDER%\examples\standalone
    cmake -B build -DEIGEN3_INSTALL_DIR="c:\%EIGEN3_INSTALL_PATH%"
    cmake --build build
    build\Debug\standalone.exe

    cd %APPVEYOR_BUILD_FOLDER%
    cmake --build build --target examples --target unit_tests --parallel 8
)
IF %COMPILER%==msys2 (
    @echo on
    SET "PATH=C:\msys64\mingw64\bin;%PATH%"

    cd %APPVEYOR_BUILD_FOLDER%
    bash -lc "cmake -B build -DWITH_EXAMPLES=ON -DWITH_UNIT_TESTS=ON -DEIGEN3_INSTALL_DIR=/c/%EIGEN3_INSTALL_PATH%"
    bash -lc "cmake --build build --target install"

    cd %APPVEYOR_BUILD_FOLDER%\examples\standalone
    bash -lc "cmake -B build -DEIGEN3_INSTALL_DIR=/c/%EIGEN3_INSTALL_PATH%"
    bash -lc "cmake --build build"
    bash -lc "./build/standalone"

    cd %APPVEYOR_BUILD_FOLDER%
    bash -lc "cmake --build build --target examples --target unit_tests --parallel 8"
)
