EXE=main.prog
TEST_EXE=tests.prog

SRC_DIR = src
OBJ_DIR = obj
INC_DIR = include
TEST_DIR = tests

SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
HEADERS = $(wildcard $(INC_DIR)/*.h)

CFLAGS   += -Wall -O3
CPPFLAGS += -I/scratch/Cuba-4.2/ -I/scratch/gsl-2.5/

LDFLAGS_GSL  = -L/scratch/gsl-2.5/.libs/ -L/scratch/gsl-2.5/cblas/.libs
LDLIBS_GSL   = -lgsl -lgslcblas
LDFLAGS_CUBA = -L/scratch/Cuba-4.2/
LDLIBS_CUBA  = -lcuba

LDLIBS += -lm

.PHONY: all clean run

all: $(EXE)

run: all
	./$(EXE)

$(EXE): main.o $(OBJ)
	$(CC) $(LDFLAGS) $(LDFLAGS_GSL) $^ $(LDLIBS) $(LDLIBS_GSL) -o $@

main.o: main.c $(HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@




clean:
	$(RM) $(OBJ) main.o
