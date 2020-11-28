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

CXXFLAGS += -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion \
			-Wcast-align -Wunused -Wlogical-op -Wnull-dereference \
			-std=c++17

all: CPPFLAGS += -DDEBUG=0 -DHAVE_INLINE -I/space/ge52sir/local/include/
all: CXXFLAGS += -O3
all: LDFLAGS  += -L/space/ge52sir/local/lib/

cluster: CPPFLAGS += -DDEBUG=0 -DHAVE_INLINE -I./local/include/
cluster: CXXFLAGS += -O3
cluster: LDFLAGS  += -L./local/lib/

debug: CPPFLAGS   += -D_GLIBCXX_DEBUG -DDEBUG=2 -I/space/ge52sir/local/include/
debug: CXXFLAGS   += -g -O0
debug: LDFLAGS    += -L/space/ge52sir/local/lib/

benchmark: CPPFLAGS += -DHAVE_INLINE
benchmark: CXXFLAGS += -O3
benchmark: LDLIBS   += -lbenchmark -pthread
profile:   CPPFLAGS += -DHAVE_INLINE
profile:   CXXFLAGS += -O3 -fno-omit-frame-pointer
profile:   LDLIBS   += -lbenchmark -pthread

LDLIBS  += -lconfig++ -lcuba -lgsl -lgslcblas

.PHONY: all clean run debug benchmark profile force

all: $(EXE)
cluster: $(EXE)
debug: $(DEBUG_EXE)
benchmark: $(BENCHMARK_EXE)
profile: $(BENCHMARK_EXE)

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
