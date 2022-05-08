@echo off
set include_folder="../../Include/"
set libraries_folder="../../Libraries/"

@echo on
g++ *.cpp glad.c -o out -I %include_folder% -L %libraries_folder% -lglfw3 -lopengl32 -mwindows
pause