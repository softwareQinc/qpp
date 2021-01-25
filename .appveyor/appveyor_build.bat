@echo off
IF %COMPILER%==msvc2019 (
    @echo on
    CALL "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"
    mkdir build
    cd build
    cmake .. -DEIGEN3_INSTALL_DIR="c:\%EIGEN3_INSTALL_PATH%"
    cmake --build . --target INSTALL
    cd %APPVEYOR_BUILD_FOLDER%\examples\standalone
    mkdir build
    cd build
    cmake .. -DEIGEN3_INSTALL_DIR="c:\%EIGEN3_INSTALL_PATH%"
    msbuild -verbosity:minimal standalone.sln
    Debug\standalone.exe
    cd %APPVEYOR_BUILD_FOLDER%/build
    msbuild -verbosity:minimal qpp.sln
)
IF %COMPILER%==msys2 (
    @echo on
    SET "PATH=C:\msys64\mingw64\bin;%PATH%"
    cd %APPVEYOR_BUILD_FOLDER%
    mkdir build
    cd build
    bash -lc "cmake .. -DEIGEN3_INSTALL_DIR=/c/%EIGEN3_INSTALL_PATH% -GNinja"
    bash -lc "sudo ninja install"
    cd %APPVEYOR_BUILD_FOLDER%\examples\standalone
    mkdir build
    cd build
    bash -lc "cmake .. -DEIGEN3_INSTALL_DIR=/c/%EIGEN3_INSTALL_PATH% -GNinja && ninja && ./standalone"
    cd %APPVEYOR_BUILD_FOLDER%\build
    bash -lc "ninja examples"
    bash -lc "ninja unit_tests"
)
