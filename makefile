all: hash diffs trail

hash: 
	gcc -g -std=c99 -Ofast -o hash `find src/hash/ -name "*.c"` `find src/hash/ -name "*.h"`

diffs:
	gcc -g -std=c99 -Ofast -o maw_diffs `find src/diffs/ -name "*.c"` `find src/diffs/ -name "*.h"`

trail:
	g++ -g -std=c++11 -march=native -Ofast -o maw_trail `find src/trail/ -name "*.cpp"` `find src/trail/ -name "*.h"` -lm

.PHONY: clean
clean: 
	rm hash maw_diffs maw_trail
