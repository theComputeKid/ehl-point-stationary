#!/bin/bash 

# Output Name
OUT="pentasolver.mexa64"

# MATLAB ROOT
MATLABROOT="/usr/local/MATLAB/R2021b"

# NVCC Path
NVCC="/usr/local/cuda/bin/nvcc"

# Include Directories
INC="-I$MATLABROOT/extern/include -I$MATLABROOT/toolbox/parallel/gpu/extern/include"

# Link Line
MATLIBDIR="$MATLABROOT/bin/glnxa64"
LDFLAGS="-shared -L$MATLIBDIR -lmwgpu -lmx -lmex -lmat -lcusparse -lcublas"

# Defines
DEFINES="-DMATLAB_MEXCMD_RELEASE=R2018a -DMX_COMPAT_64 -DMATLAB_MEX_FILE"

# Compiler Flags
CXXFLAGS="-x cu -std=c++17 -Xcompiler -fPIC"
CXXDEBUGFLAGS="-g -G $CXXFLAGS"
CXXRELEASEFLAGS="-O2 $CXXFLAGS -Xcompiler -O2"

if [ "$1" = "debug" ]; then
	$NVCC $CXXDEBUGFLAGS $DEFINES $INC pentasolver.cpp -o $OUT $LDFLAGS
else
	$NVCC $CXXRELEASEFLAGS $DEFINES $INC pentasolver.cpp -o $OUT $LDFLAGS
fi
