# Machine-specific settings (CXX, include/library paths, ...), see config.mk.example
-include config.mk

SRC_DIR   = src
INC_DIR   = include
BUILD_DIR = build

# Build configuration: release, debug or profile.
# Separate obj dir for each, so switching between them does not trigger rebuilds.
CONFIG ?= release

ifeq ($(CONFIG),release)
  CXXFLAGS += -O3 -DDEBUG=0 -DHAVE_INLINE
  SUFFIX    =
else ifeq ($(CONFIG),debug)
  CXXFLAGS += -g -O0 -DDEBUG=2 -D_GLIBCXX_DEBUG \
	-D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_DEBUG
  SUFFIX    = -debug
else ifeq ($(CONFIG),profile)
  CXXFLAGS += -O3 -g -fno-omit-frame-pointer -DDEBUG=0 -DHAVE_INLINE
  SUFFIX    = -profile
else
  $(error Unknown CONFIG '$(CONFIG)', expected release, debug or profile)
endif

OBJ_DIR       = obj/$(CONFIG)
EXE           = $(BUILD_DIR)/cosPT$(SUFFIX)
BENCHMARK_EXE = $(BUILD_DIR)/bench$(SUFFIX)

SRC     = $(wildcard $(SRC_DIR)/*.cpp)
OBJ     = $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
HEADERS = $(wildcard $(INC_DIR)/*.hpp)

CXXFLAGS += -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion \
	-Wcast-align -Wunused -Wnull-dereference -std=c++17 \
  -Wlogical-op

# Generate header dependency files (.d) alongside the objects
CPPFLAGS += -MMD -MP
LDLIBS   += -lconfig++ -lcuba -lgsl -lgslcblas

# Git hash written to output files, with "-dirty" appended if there are
# uncommitted changes
GIT_SHA := $(shell git describe --always --dirty --abbrev=40 --exclude='*' 2>/dev/null || echo unknown)

.PHONY: all debug benchmark profile clean force

all: $(EXE)

benchmark: $(BENCHMARK_EXE)

debug:
	@$(MAKE) --no-print-directory CONFIG=debug

profile:
	@$(MAKE) --no-print-directory CONFIG=profile benchmark

$(BENCHMARK_EXE): LDLIBS += -lbenchmark -pthread

$(EXE): $(OBJ_DIR)/main.o $(OBJ) | $(BUILD_DIR)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(BENCHMARK_EXE): $(OBJ_DIR)/benchmark.o $(OBJ) | $(BUILD_DIR)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(OBJ_DIR)/compiler_flags | $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: %.cpp $(OBJ_DIR)/compiler_flags | $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# The build time is that of the last compilation of version.o, which is
# recompiled whenever any source file or the git hash changes (but not on a
# no-op make)
$(OBJ_DIR)/version.o: $(SRC_DIR)/version.cpp $(OBJ_DIR)/compiler_flags \
		$(OBJ_DIR)/git_sha $(SRC) main.cpp benchmark.cpp $(HEADERS) | $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DCOSPT_GIT_SHA='"$(GIT_SHA)"' \
		-DCOSPT_BUILD_TIME="\"$$(date '+%d.%m.%Y %R %Z')\"" -c $< -o $@

# These files are only rewritten when their content changes, so that objects
# are rebuilt exactly when the compiler flags or git hash change
$(OBJ_DIR)/compiler_flags: force | $(OBJ_DIR)
	@echo '$(CXX) $(CPPFLAGS) $(CXXFLAGS)' | cmp -s - $@ || echo '$(CXX) $(CPPFLAGS) $(CXXFLAGS)' > $@

$(OBJ_DIR)/git_sha: force | $(OBJ_DIR)
	@echo '$(GIT_SHA)' | cmp -s - $@ || echo '$(GIT_SHA)' > $@

$(BUILD_DIR) $(OBJ_DIR):
	mkdir -p $@

clean:
	$(RM) -r obj
	$(RM) $(addprefix $(BUILD_DIR)/, cosPT cosPT-debug cosPT-profile \
		bench bench-debug bench-profile)

-include $(OBJ:.o=.d) $(OBJ_DIR)/main.d $(OBJ_DIR)/benchmark.d
