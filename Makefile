EXE=cosPT

SRC_DIR = src
OBJ_DIR = obj
INC_DIR = include
TEST_DIR = tests

SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
HEADERS = $(wildcard $(INC_DIR)/*.h)


CFLAGS += -Wall -Wextra -Wpedantic
CPPFLAGS += -I/scratch/Cuba-4.2/ -I/scratch/gsl-2.5/ $(OPTIONS)

# Information about git revision and compile time
GIT_HASH = `git rev-parse --short HEAD`
CPPFLAGS += -DGIT_HASH="\"$(GIT_HASH)\""

LDFLAGS_GSL  = -L/scratch/gsl-2.5/.libs/ -L/scratch/gsl-2.5/cblas/.libs
LDLIBS_GSL   = -lgsl -lgslcblas
LDFLAGS_CUBA = -L/scratch/Cuba-4.2/
LDLIBS_CUBA  = -lcuba

LDFLAGS += $(LDFLAGS_GSL) $(LDFLAGS_CUBA)
LDLIBS += $(LDLIBS_GSL) $(LDLIBS_CUBA) -lm

all:   CFLAGS   += -O3
1loop: CPPFLAGS += -DDEBUG=0 -DLOOPS=1
1loop: CFLAGS   += -O3
2loop: CPPFLAGS += -DDEBUG=0 -DLOOPS=2
2loop: CFLAGS   += -O3
debug: CPPFLAGS += -DDEBUG=2
debug: CFLAGS   += -O0 -g

.PHONY: all clean run 1loop 2loop

all: $(EXE)
1loop: $(EXE)
2loop: $(EXE)
debug: $(EXE)

run: all
	./$(EXE)

$(EXE): main.o $(OBJ)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

main.o: main.c $(HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJ) main.o
