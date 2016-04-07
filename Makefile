CCFLAGS = -O2

all: Needleman Hirschberg

Needleman: Needleman.o
		g++ $(CCFLAGS) -o Needleman Needleman.o

Hirschberg: Hirschberg.o
		g++ $(CCFLAGS) -o Hirschberg Hirschberg.o

Needleman.o: Needleman.cpp
		g++ $(CCFLAGS) -c Needleman.cpp

Hirschberg.o: Hirschberg.cpp
		g++ $(CCFLAGS) -c Hirschberg.cpp


clean:
		rm -rf *.o *~ *.x
