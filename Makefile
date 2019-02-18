EXE=main.prog
TEST_EXE=test.prog

SRC_DIR = src
OBJ_DIR = obj

SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

CFLAGS   += -Wall -O2
CXXFLAGS += -Wall -O2
CPPFLAGS += -I/scratch/Cuba-4.2/ -I/scratch/gsl-2.5/

LDFLAGS_GSL  = -L/scratch/gsl-2.5/.libs/ -L/scratch/gsl-2.5/cblas/.libs
LDLIBS_GSL   = -lgsl -lgslcblas
LDFLAGS_CUBA = -L/scratch/Cuba-4.2/
LDLIBS_CUBA  = -lcuba

LDLIBS += -lm

.PHONY: all clean tests

all: $(EXE)

run: all
	./$(EXE)

$(EXE): $(OBJ)
	$(CC) $(LDFLAGS) $(LDFLAGS_GSL) $^ $(LDLIBS) $(LDLIBS_GSL) -o $@

$(OBJ_DIR)/main.o: $(SRC_DIR)/main.c include/constants.h \
	include/utilities.h include/kernels.h include/spt_kernels.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/utilities.o: $(SRC_DIR)/utilities.c include/constants.h include/utilities.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/kernels.o: $(SRC_DIR)/kernels.c include/constants.h \
	include/utilities.h include/kernels.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/spt_kernels.o: $(SRC_DIR)/spt_kernels.c include/constants.h \
	include/utilities.h include/kernels.h include/spt_kernels.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@


clean:
	$(RM) $(OBJ) $(EXE)
