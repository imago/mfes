bin_PROGRAMS=mfes 
mfes_SOURCES= mfes.cpp lib/LSMS/grid.cpp lib/LSMS/sceneBuilder.cpp  
AM_CPPFLAGS =  -std=c++0x -Wall -O2 -g -ftemplate-depth-999 -funroll-loops -finline-limit=10000 -fopenmp -DPARALLEL -DUSE_TIMEOFDAY -DLAPACK -I${NETGENDIR}/../include    -Wno-unused-result  -I./lib/vcglib -I${NETGEN_SOURCE}/libsrc/include/
AM_LDFLAGS =   -fopenmp -DPARALLEL -DUSE_TIMEOFDAY -DLAPACK -L/usr/local/lib -lboost_system -lboost_filesystem -lboost_program_options -lngstd -lngsolve -lngbla -lngcomp -lngfem -lngla -lnglib -lngmg  -llapack -lblas
