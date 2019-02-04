CC=gcc
CFLAGS=-Wall -O2 -I/scratch/Cuba-4.2/ -I/scratch/gsl-2.5/

LDGSL=-L/scratch/gsl-2.5/.libs/ -L/scratch/gsl-2.5/cblas/.libs -lgsl -lgslcblas
LDCUBA=-L/scratch/Cuba-4.2/ -lcuba
LDFLAGS=$(LDGSL) -lm

EXECUTABLE=main.prog

all: main.o utilities.o
	$(CC) -o $(EXECUTABLE) $^ $(LDFLAGS)

run: all
	./$(EXECUTABLE)

main.o: main.c utilities.h

utilities.o: utilities.c utilities.h

clean:
	rm -rf $(OBJECT) $(EXECUTABLE)
