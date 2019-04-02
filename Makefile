EXE=main.prog

SRC_DIR = src
OBJ_DIR = obj
INC_DIR = include
TEST_DIR = tests

SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
HEADERS = $(wildcard $(INC_DIR)/*.h)

TEST_INTERFACE = $(TEST_DIR)/test_interface

CFLAGS += -Wall -Wextra -Wpedantic
CPPFLAGS += -I/scratch/Cuba-4.2/ -I/scratch/gsl-2.5/

LDFLAGS_GSL  = -L/scratch/gsl-2.5/.libs/ -L/scratch/gsl-2.5/cblas/.libs
LDLIBS_GSL   = -lgsl -lgslcblas
LDFLAGS_CUBA = -L/scratch/Cuba-4.2/
LDLIBS_CUBA  = -lcuba

LDLIBS += -lm

all:            CFLAGS += -O3
1loop:          CFLAGS += -O3 -DDEBUG=0 -DLOOPS=1
2loop:          CFLAGS += -O3 -DDEBUG=0 -DLOOPS=2
debug:          CFLAGS += -O0 -DDEBUG=2 -g
test_interface: CFLAGS += -fPIC

.PHONY: all clean run 1loop 2loop test_interface

all: $(EXE)
1loop: $(EXE)
2loop: $(EXE)
debug: $(EXE)

run: all
	./$(EXE)

$(EXE): main.o $(OBJ)
	$(CC) $(LDFLAGS) $(LDFLAGS_GSL) $(LDFLAGS_CUBA) $^ $(LDLIBS_GSL) $(LDLIBS_CUBA) $(LDLIBS) -o $@

main.o: main.c $(HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

test_interface: $(TEST_INTERFACE).so

$(TEST_INTERFACE).so: $(TEST_INTERFACE).o $(OBJ)
	$(CC) $(LDFLAGS) $(LDFLAGS_GSL) $^ $(LDLIBS) $(LDLIBS_GSL) -shared -o $@

$(TEST_INTERFACE).o: $(TEST_INTERFACE).c $(HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

diff_eq:
	$(CC) $(CPPFLAGS) $(CFLAGS) gsl_diff_eq.c \
		$(LDFLAGS) $(LDFLAGS_GSL) $^ $(LDLIBS_GSL) $(LDLIBS) -o diff_eq.prog
	./diff_eq.prog

clean:
	$(RM) $(OBJ) main.o $(TEST_INTERFACE).o
