EXE=cosPT

SRC_DIR = src
OBJ_DIR = obj
INC_DIR = include

SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
HEADERS = $(wildcard $(INC_DIR)/*.h)

GIT = git

CFLAGS += -Wall -Wextra -Wpedantic
CPPFLAGS += -I/scratch/Cuba-4.2/ -I/scratch/gsl-2.5/ $(OPTIONS)

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

.PHONY: all clean run 1loop 2loop force

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

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS) compiler_flags
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(SRC_DIR)/version.c: force
	$(GIT) rev-parse HEAD | awk ' BEGIN {print "#include \"../include/version.h\""} {print "const char * build_git_sha = \"" $$0"\";"} END {}' > $@
	date '+%d.%m.%Y %R %Z' | awk 'BEGIN {} {print "const char * build_git_time = \""$$0"\";"} END {} ' >> $@

compiler_flags: force
	echo '$(CPPFLAGS) $(CFLAGS)' | cmp -s - $@ || echo '$(CPPFLAGS) $(CFLAGS)' > $@

clean:
	$(RM) $(OBJ) main.o
