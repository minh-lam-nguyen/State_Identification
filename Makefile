CC=g++ -g -Wall -Wextra

all: SyncTree_hashset

SyncTree_hashset: obj/SyncTree.o obj/data_structures/set.o obj/data_structures/fsm.o obj/data_structures/mergseq.o
	$(CC) -o bin/$@ $^

SyncTree_dhashset: obj/SyncTree.o obj/data_structures/set_dynamic.o obj/data_structures/fsm.o obj/data_structures/mergseq.o
	$(CC) -o bin/$@ $^

SyncTree_vector: obj/SyncTree.o obj/data_structures/set_vector.o
	$(CC) -o bin/$@ $^

SyncTree_stlset: obj/SyncTree.o obj/data_structures/set_stl.o
	$(CC) -o bin/$@ $^

SyncTree_unorderedset: obj/SyncTree.o obj/data_structures/unordered_set.o
	$(CC) -o bin/$@ $^

test_hashset: obj/test.o obj/data_structures/set.o
	$(CC) -o bin/$@ $^

test_vector: obj/test.o obj/data_structures/set_vector.o
	$(CC) -o bin/$@ $^

test_stlset: obj/test.o obj/data_structures/set_stl.o
	$(CC) -o bin/$@ $^

obj/%.o: src/%.c
	$(CC) -o $@ -c $^

clean:
	rm bin/*
	rm obj/data_structures/*.o
	rm obj/*.o
