bin_PROGRAMS=mfes
mfes_SOURCES= mfes.cpp lib/LSMS/grid.cpp lib/LSMS/sceneBuilder.cpp 
AM_CPPFLAGS =  -Wall -O2 -g -ftemplate-depth-999 -funroll-loops -finline-limit=10000 -fopenmp -DPARALLEL -DUSE_TIMEOFDAY -DLAPACK -I/opt/netgen/include    -Wno-unused-result -I/opt/netgen/include/ -I./lib/vcglib -I/home/parallels/ngsolve/sources/MUMPS_4.10.0/include/
AM_LDFLAGS =   -fopenmp -DPARALLEL -DUSE_TIMEOFDAY -DLAPACK -L/usr/local/lib -lboost_system -lboost_filesystem -lboost_program_options -lngstd -lngsolve -lngbla -lngcomp -lngfem -lngla -lnglib -lngmg  -llapack -lblas
