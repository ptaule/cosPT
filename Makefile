EXE=build/cosPT
DEBUG_EXE=build/debug
BENCHMARK_EXE=build/bench

SRC_DIR = src
OBJ_DIR = obj
INC_DIR = include

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
HEADERS = $(wildcard $(INC_DIR)/*.hpp)

GIT = git
CXX ?= g++

CXXFLAGS += -Wall -Wextra -Wpedantic -std=c++14
CPPFLAGS += -I/scratch/Cuba-4.2/ -I/scratch/gsl-2.5/

LDFLAGS_GSL  = -L/scratch/gsl-2.5/.libs/ -L/scratch/gsl-2.5/cblas/.libs
LDLIBS_GSL   = -lgsl -lgslcblas
LDFLAGS_CUBA = -L/scratch/Cuba-4.2/
LDLIBS_CUBA  = -lcuba

LDFLAGS += $(LDFLAGS_GSL) $(LDFLAGS_CUBA)
LDLIBS += $(LDLIBS_GSL) $(LDLIBS_CUBA)

all:       CXXFLAGS += -O3
debug:     CPPFLAGS += -D_GLIBCXX_DEBUG -DDEBUG=2
debug:     CXXFLAGS += -g3 -Og
benchmark: CXXFLAGS += -O3
benchmark: LDLIBS   += -lbenchmark -pthread
profile:   CXXFLAGS += -O3 -fno-omit-frame-pointer
profile:   LDLIBS   += -lbenchmark -pthread

.PHONY: all clean run benchmark profile force

all: $(EXE)
debug: $(DEBUG_EXE)
benchmark: $(BENCHMARK_EXE)
profile: $(BENCHMARK_EXE)

run: all
	./$(EXE)

$(EXE): main.o $(OBJ)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(DEBUG_EXE): main.o $(OBJ)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(BENCHMARK_EXE): benchmark.o $(OBJ)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

main.o: main.cpp $(HEADERS) compiler_flags
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

benchmark.o: benchmark.cpp $(HEADERS) compiler_flags
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS) compiler_flags
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(SRC_DIR)/version.cpp: force
	$(GIT) rev-parse HEAD | awk ' BEGIN {print "#include \"../include/version.hpp\""} {print "std::string build_git_sha = \"" $$0"\";"} END {}' > $@
	date '+%d.%m.%Y %R %Z' | awk 'BEGIN {} {print "std::string build_git_time = \""$$0"\";"} END {} ' >> $@

compiler_flags: force
	echo '$(CPPFLAGS) $(CXXFLAGS)' | cmp -s - $@ || echo '$(CPPFLAGS) $(CXXFLAGS)' > $@

clean:
	$(RM) $(OBJ) main.o benchmark.o $(EXE) $(DEBUG_EXE) $(BENCHMARK_EXE)
