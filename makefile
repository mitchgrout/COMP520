all: hash diffs trail

hash: 
	gcc -g -std=c99 -o hash `find src/hash/ -name "*.c"` `find src/hash/ -name "*.h"`

diffs:
	gcc -g -std=c99 -o maw_diffs `find src/diffs/ -name "*.c"` `find src/diffs/ -name "*.h"`

trail:
	g++ -g -std=c++11 -o maw_trail `find src/trail/ -name "*.cpp"` `find src/trail/ -name "*.h"` -lm

.PHONY: clean
clean: 
	rm hash maw_diffs maw_trail
