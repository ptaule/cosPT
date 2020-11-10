EXE=build/cosPT
DEBUG_EXE=build/debug
BENCHMARK_EXE=build/bench

SRC_DIR = src
OBJ_DIR = obj
INC_DIR = include

SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
HEADERS = $(wildcard $(INC_DIR)/*.h)

GIT = git

CFLAGS += -Wall -Wextra -Wpedantic
CPPFLAGS += $(OPTIONS)

all: CPPFLAGS += -DDEBUG=0 -I/space/ge52sir/local/include/
all: CFLAGS += -O3
all: LDFLAGS  += -L/space/ge52sir/local/lib/

debug: CPPFLAGS += -DDEBUG=2 -DN_CORES=0 -I/space/ge52sir/local/include/
debug: CFLAGS   += -g -O0
debug: LDFLAGS  += -L/space/ge52sir/local/lib/

benchmark: CXXFLAGS += -O3
benchmark: LDLIBS   += -lbenchmark -pthread

LDLIBS  += -lcuba -lgsl -lgslcblas -lm

.PHONY: all clean run 1loop 2loop benchmark force

all: $(EXE)
1loop: $(EXE)
2loop: $(EXE)
debug: $(EXE)
benchmark: $(BENCHMARK_EXE)

run: all
	./$(EXE)

$(EXE): main.o $(OBJ)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(BENCHMARK_EXE): benchmark.o $(OBJ)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

benchmark.o: benchmark.cpp $(HEADERS) compiler_flags
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

main.o: main.c $(HEADERS) compiler_flags
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
