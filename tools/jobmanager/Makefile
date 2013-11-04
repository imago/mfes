OBJS1 = analyze.o

CC = g++ 
DEBUG = 
CFLAGS = -Wall -c $(DEBUG) 
LFLAGS = -Wall $(DEBUG) 
all: analyze

analyze: $(OBJS1)
	$(CC) $(LFLAGS) $(OBJS1) -o $@

analyze.o: analyze.cpp
	$(CC) $(CFLAGS) -O2 -c analyze.cpp

clean:
	\rm -f *.o *.so


