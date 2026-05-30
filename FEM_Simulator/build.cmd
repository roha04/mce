@echo off
set "CMAKE=C:\Program Files\CMake\bin\cmake.exe"
if not exist build mkdir build
cd /d "%~dp0build"
"%CMAKE%" ..
"%CMAKE%" --build . --config Release
if %ERRORLEVEL% equ 0 (
    echo.
    echo Build OK: build\Release\FEM_App.exe
)
