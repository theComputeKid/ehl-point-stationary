@echo off

@REM Output Name
SET OUT=pentasolver.mexw64

@REM Include Directories
SET INC=-I%MATLABROOT%\extern\include -I%MATLABROOT%\toolbox\parallel\gpu\extern\include

@REM Link Line
SET MATLIBDIR=%MATLABROOT%\extern\lib\win64\microsoft
SET LDFLAGS=-shared -L%MATLIBDIR% -llibmx -llibmex -llibmat -lgpu -lcusparse -lcublas -Xlinker -EXPORT:mexFunction -Xlinker -noimplib

@REM Defines
SET DEFINES=-DMATLAB_MEXCMD_RELEASE=R2018a -DMX_COMPAT_64 -DMATLAB_MEX_FILE

@REM Compiler Flags
SET CXXFLAGS=-x cu -std=c++17
SET CXXDEBUGFLAGS=-g -G %CXXFLAGS%
SET CXXRELEASEFLAGS=-O2 %CXXFLAGS% -Xcompiler -O2

IF "%~1"=="debug" (
	nvcc %CXXDEBUGFLAGS% %DEFINES% %INC% pentasolver.cpp -o %OUT% %LDFLAGS%
) ELSE (
	nvcc %CXXRELEASEFLAGS% %DEFINES% %INC% pentasolver.cpp -o %OUT% %LDFLAGS%
)
