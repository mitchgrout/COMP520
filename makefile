CFLAGS=-g -std=c99 -march=native -Ofast
CPPFLAGS=-g -std=c++11 -march=native -Ofast

all: hash diffs trail

hash: 
	gcc $(CFLAGS) -o hash `find src/hash/ -name "*.c"` `find src/hash/ -name "*.h"`

diffs:
	gcc $(CFLAGS) -o maw_diffs `find src/diffs/ -name "*.c"` `find src/diffs/ -name "*.h"`

trail:
	g++ $(CPPFLAGS) -o maw_trail `find src/trail/ -name "*.cpp"` `find src/trail/ -name "*.hpp"` -lm

trail_gen:
	g++ $(CPPFLAGS) -o maw_trail_gen `find src/trail_gen/ -name "*.cpp"` -lm

.PHONY: clean
clean: 
	rm hash maw_diffs maw_trail maw_trail_gen
