OBJS1 = writePotatAscii.o
OBJS2 = energydiff.o
OBJS3 = pointcharges.o
OBJS4 = precalculation.o
CC = ngscxx
DEBUG = 
CFLAGS = -Wall -c $(DEBUG) -I/scratch/scratch/sakalli/ngsolve/sources/MUMPS_4.10.0/include/
LFLAGS = -Wall $(DEBUG)

all: writePotatAscii.so energydiff.so pointcharges.so precalculation.so

writePotatAscii.so: $(OBJS1)
	$(CC) -shared -fopenmp -fpic $(LFLAGS) $(OBJS1) -lngsolve -o $@

energydiff.so: $(OBJS2)
	$(CC) -shared -fopenmp -fpic $(LFLAGS) $(OBJS2) -lngsolve -o $@

pointcharges.so: $(OBJS3)
	$(CC) -shared -fopenmp -fpic $(LFLAGS) $(OBJS3) -lngsolve -o $@

precalculation.so: ${OBJS4}
	$(CC) -shared -fopenmp -fpic $(LFLAGS) $(OBJS4) -lngsolve -o $@

writePotatAscii.o: writePotatAscii.cpp
	$(CC) $(CFLAGS) -O2 -fopenmp -fpic -DNETGEN_ELTRANS -DUSE_TIMEOFDAY -DLAPACK -I. -I$(NETGENDIR)/../include -I/usr/include/tcl8.5/ -c writePotatAscii.cpp

energydiff.o: energydiff.cpp 
	$(CC) $(CFLAGS) -O2 -fopenmp -fpic  -DNETGEN_ELTRANS -DUSE_TIMEOFDAY -DLAPACK -I. -I$(NETGENDIR)/../include -I/usr/include/tcl8.5/ -c energydiff.cpp

pointcharges.o: pointcharges.cpp 
	$(CC) $(CFLAGS) -O2 -fopenmp -fpic  -DNETGEN_ELTRANS -DUSE_TIMEOFDAY -DLAPACK -I. -I$(NETGENDIR)/../include -I/usr/include/tcl8.5/ -c pointcharges.cpp

precalculation.o: precalculation.cpp 
	$(CC) $(CFLAGS) -O2 -fopenmp -fpic  -DNETGEN_ELTRANS -DUSE_TIMEOFDAY -DLAPACK -I. -I$(NETGENDIR)/../include -I/usr/include/tcl8.5/ -c precalculation.cpp


clean:
	\rm -f *.o *.so


