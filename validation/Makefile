PROGRAM := a255885-independent
SOURCE := a255885_independent.c
CC := gcc
WARNINGS := -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wdouble-promotion -Wformat=2 -Werror
COMMON := -std=c17 $(WARNINGS) -fopenmp
RELEASE_FLAGS := -O3 -flto -funroll-loops -fomit-frame-pointer -march=native -mtune=native
DEBUG_FLAGS := -O1 -g3 -fno-omit-frame-pointer

.PHONY: all release debug asan ubsan test clean

all: release

release: $(PROGRAM)

$(PROGRAM): $(SOURCE) Makefile
	$(CC) $(COMMON) $(RELEASE_FLAGS) $(SOURCE) -o $@

debug:
	$(CC) $(COMMON) $(DEBUG_FLAGS) $(SOURCE) -o $(PROGRAM)-debug

asan:
	$(CC) $(COMMON) $(DEBUG_FLAGS) -fsanitize=address,undefined $(SOURCE) -o $(PROGRAM)-asan

ubsan:
	$(CC) $(COMMON) $(DEBUG_FLAGS) -fsanitize=undefined $(SOURCE) -o $(PROGRAM)-ubsan

test: release
	./$(PROGRAM) self-test 2024

clean:
	rm -f $(PROGRAM) $(PROGRAM)-debug $(PROGRAM)-asan $(PROGRAM)-ubsan
