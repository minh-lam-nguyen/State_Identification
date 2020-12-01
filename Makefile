CC=gcc -Wall -Wextra

all: SyncTree

SyncTree: obj/SyncTree.o obj/data_structures/set.o
	$(CC) -o bin/$@ $^

test: obj/test.o obj/data_structures/set.o
	$(CC) -o bin/$@ $^

obj/%.o: src/%.c
	$(CC) -o $@ -c $^
