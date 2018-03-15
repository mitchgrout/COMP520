CC=gcc
CFLAGS=-g -std=c99

all: hash diffs

hash: 
	$(CC) $(CFLAGS) -o hash `find src/hash/ -name "*.c"` `find src/hash/ -name "*.h"`

diffs:
	$(CC) $(CFLAGS) -o maw_diffs `find src/diffs/ -name "*.c"` `find src/diffs/ -name "*.h"`

.PHONY: clean
clean: 
	rm hash maw_diffs
