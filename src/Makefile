PROGRAM := a255885
SOURCE := a255885.c

ifeq ($(origin CC), default)
CC := gcc
endif
PYTHON ?= python3

SOURCE_SHA256 := $(shell sha256sum $(SOURCE) | cut -d' ' -f1)
MAKEFILE_SHA256 := $(shell sha256sum Makefile | cut -d' ' -f1)
IDENTITY_CPPFLAGS := -DSOURCE_SHA256=\"$(SOURCE_SHA256)\" -DMAKEFILE_SHA256=\"$(MAKEFILE_SHA256)\"

WARNINGS := -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wdouble-promotion -Wformat=2 -Werror
COMMON := -std=c17 $(WARNINGS) -fopenmp $(IDENTITY_CPPFLAGS)
RELEASE_FLAGS := -O3 -flto -funroll-loops -fomit-frame-pointer -march=native -mtune=native
DEBUG_FLAGS := -O1 -g3 -fno-omit-frame-pointer
LDLIBS := -lm

.PHONY: all release debug asan ubsan tsan clang-check test clean identity

all: release

release: $(PROGRAM)

$(PROGRAM): $(SOURCE) Makefile
	$(CC) $(COMMON) $(RELEASE_FLAGS) $(SOURCE) -o $@ $(LDLIBS)

debug:
	$(CC) $(COMMON) $(DEBUG_FLAGS) $(SOURCE) -o $(PROGRAM)-debug $(LDLIBS)

asan:
	$(CC) $(COMMON) $(DEBUG_FLAGS) -fsanitize=address,undefined $(SOURCE) -o $(PROGRAM)-asan $(LDLIBS)

ubsan:
	$(CC) $(COMMON) $(DEBUG_FLAGS) -fsanitize=undefined $(SOURCE) -o $(PROGRAM)-ubsan $(LDLIBS)

tsan:
	$(CC) $(COMMON) $(DEBUG_FLAGS) -fsanitize=thread $(SOURCE) -o $(PROGRAM)-tsan $(LDLIBS)

clang-check:
	clang $(COMMON) -O2 $(SOURCE) -o $(PROGRAM)-clang $(LDLIBS)

test: release
	$(PYTHON) ../validation/validate_a255885.py --binary ./$(PROGRAM)

identity:
	@sha256sum $(SOURCE) Makefile $(PROGRAM)

clean:
	rm -f $(PROGRAM) $(PROGRAM)-debug $(PROGRAM)-asan $(PROGRAM)-ubsan $(PROGRAM)-tsan $(PROGRAM)-clang
