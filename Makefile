CC=gcc -Wall -Wextra

all: SyncTree

SyncTree: src/SyncTree.o src/set.o
	$(CC) -o $@ $^

test: src/test.o src/set.o
	$(CC) -o $@ $^
