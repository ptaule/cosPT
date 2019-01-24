CC=gcc
CFLAGS=-I/scratch/Cuba-4.2/
LDCUBA=-L/scratch/Cuba-4.2/ -lcuba
LDFLAGS=$(LDCUBA) -lm

SOURCE=main.c
OBJECT=main.o
EXECUTABLE=main.prog


all: $(OBJECT)
	$(CC) -o $(EXECUTABLE) $^ $(LDFLAGS)

run: all
	./$(EXECUTABLE)

$(OBJECT): $(SOURCE)

.PHONY: clean
clean:
	rm -rf $(OBJECT) $(EXECUTABLE)
