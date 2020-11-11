EXE=build/cosPT
DEBUG_EXE=build/debug

SRC_DIR = src
OBJ_DIR = obj
INC_DIR = include

SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ = $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)
HEADERS = $(wildcard $(INC_DIR)/*.h)

GIT = git

CFLAGS += -Wall -Wextra -Wpedantic #-Wconversion

all: CPPFLAGS += -DDEBUG=0 -I/space/ge52sir/local/include/
all: CFLAGS   += -O3
all: LDFLAGS  += -L/space/ge52sir/local/lib/

cluster: CPPFLAGS += -DDEBUG=0 -I./local/include/
cluster: CFLAGS   += -O3
cluster: LDFLAGS  += -L./local/lib/

debug: CPPFLAGS += -DDEBUG=2 -I/space/ge52sir/local/include/
debug: CFLAGS   += -g -O0
debug: LDFLAGS  += -L/space/ge52sir/local/lib/

CPPFLAGS += $(OPTIONS)

LDLIBS  += -lcuba -lgsl -lgslcblas -lm

.PHONY: all clean run debug force

all: $(EXE)
cluster: $(EXE)
debug: $(DEBUG_EXE)

$(EXE): main.o $(OBJ)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(DEBUG_EXE): main.o $(OBJ)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

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
	$(RM) $(OBJ) main.o $(EXE) $(DEBUG_EXE)
