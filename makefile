CC=gcc
CFLAGS=-g -std=c99

all: hash

hash: 
	$(CC) $(CFLAGS) -o hash `find src/ -name "*.c"` `find src/ -name "*.h"`

.PHONY: clean
clean: 
	rm hash
