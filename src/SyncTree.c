#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <queue>
#include "data_structures/set.h"
#include "data_structures/fsm.h"
#include "data_structures/mergseq.h"

#define SIZE 150 // res
#define SIZE_Q 100000000

#define INT_SIZE (int)(sizeof(int)*8-1)


// free queue and visited
void clear_Q (int ** q, int taille_q) {
    int i;
    for (i=0; i<taille_q; i++) {
            free(q[i]);
    }
    free(q);
}


int* initNode(FSM fsm){
    // Each node is a list of fsm->s bits
    // If the i-th lowest bit is 1, then the i-th state is in the node.
    // Each of the first bits is hidden in a INT_SIZE byte integer. INT_SIZE is defined as the number of bits needed to encode an integer minus 1 (for the sign bit)
    // For instance, a node containing the states 0, 1, 3 and 5 of a 7 nodes automata would be encoded with the byte 0b00101011.
    // As an integer is encoded with INT_SIZE bits, that number woud be completed with INT_SIZE - 8 zeros on the left
    // If the automata has more than INT_SIZE states, we use two or more integers. The last bit of the first integer is the state 0, the second last is 1, ... the first is the INT_SIZE - 1 -th state, ...
    //
    // Finally, a state also contains a sequence of 1 and 0.
    // That list of states is then followed by a number n and n bits.
    //
    // This function returns a node corresponding to all the states : so it returns a list of integers 0b1111...1111 0b1111...1111 ... 0b1111...1111 0b000...001111...111 The number of ones in the last integer
    // depends on the number of states in the automata : fsm->s % INT_SIZE.
    // The following sequence is of size 0 and is initialized with 0.
    //

    int size = get_s(fsm) / INT_SIZE + 1; // Number of integers needed to encode all the bits
    int* state = (int*)malloc( (SIZE + 1 + size) * sizeof(int) );

    // The initial states contains
    for(int j = 0; j < size - 1; j++){
        state[j] = ((unsigned int)(~0)) >> 1; // Only ones except the first
    }
    int last_size = get_s(fsm) % INT_SIZE;
    state[size - 1] = ((unsigned int)(~0)) >> (INT_SIZE - last_size + 1); // 2^(INT_SIZE - lastsize) - 1

    // Init the size of the sequence following the states and the sequence itself with 0
    for(int j=size; j<SIZE + 1 + size; j++)
        state[j] = 0;
    return state;
}


int* successor(FSM fsm, int* node, int input){
    int size =get_s(fsm) / INT_SIZE + 1;
    int* neigh = (int*)malloc( (SIZE + 1 + size) * sizeof(int) );
    
    int j;
    for (j=0; j<SIZE + 1 + size; j++) {
        neigh[j] = 0;
    }

    int index = 0;
    for(j = 0; j < size; j++){
        int x = node[j];
        while(x != 0){
            if(x % 2 == 1){ // Last bit of x
                int succ = get_succ(fsm, index, input);
                neigh[succ / INT_SIZE] |= (1 << (succ % INT_SIZE));
            }
            index += 1;
            x /= 2;
        }
    }

    return neigh;
}


// Synchronizing tree
int syncTree(FSM fsm, int* res) {

	// queue
	int **queue = (int**)malloc( SIZE_Q * sizeof(int*) );
	int front_q = 0;
	int back_q = 1;
    int encoding_size = get_s(fsm) / INT_SIZE + 1;

	// visited set of states
    Set visited = initSet(encoding_size);

	if (!(queue))
		return -1;

	// add set of all state (init) on queue and visited

	queue[0] = initNode(fsm);

	// bfs
	while (front_q != back_q) {

		int current = front_q;
		front_q++; 


        int* states = queue[current]; 

        int added = add(visited, states);
        if(added == 0)
            continue;
        if(added == -1){
            printf("ERREUR CAPACITE ATTEINTE\n");
            break;
        }

        // size of current seq
            
        int k, tmp_cpt = 0;
        for(int i = 0; i < encoding_size; i++){
            int x = states[i];
            while(x != 0){
                if(x % 2 == 1){
                    tmp_cpt++;
                    if(tmp_cpt == 2)
                        goto end;
                }
                x /= 2;
            }
        }
        end:

        int seq_size = states[encoding_size];

        // if singleton (1 and only one state) : SS
        if (tmp_cpt==1) {
            for(k=0; k<seq_size ;k++)
                res[k] = states[encoding_size+1+k];
            clear_Q(queue, back_q);

            freeSet(visited);
            return seq_size;
        }


		//printf("current: %d / ttl_q: %d \n", current, current + count_q+1);

		// c0 and c1 successors of current
		int *c0 = successor(fsm, states, 0);
		int *c1 = successor(fsm, states, 1);
        int* succs[2] = {c0, c1};


        for(int symbol = 0; symbol < 2; symbol++){
            int* successor = succs[symbol];
            int inVisited = find(visited, successor);

            if(inVisited == 1){
                free(successor);
                continue;
            }
	
            queue[back_q] = successor;
            queue[back_q][encoding_size] = seq_size+1;

            for(k=0; k<seq_size; k++) {
                queue[back_q][k+ encoding_size +1] = states[k+ encoding_size +1];
            }
            queue[back_q][k + encoding_size + 1] = symbol;
            back_q++;

        }
	}

	return -1;
}

struct CompareQ
{
    bool operator ()(int* e, int* f)
    {
        return e[0] > f[0];
    }
};

void clear_priorityQ(std::priority_queue<int*, std::vector<int*>, CompareQ> queue) {
	while (!queue.empty()) {
		free(queue.top());
		queue.pop();
	}
}

// best Sync tree
int bestSearch(FSM fsm, int *res, MergSeq mergs){

	std::priority_queue<int*, std::vector<int*>, CompareQ> queue;
	// visited set of states
   	Set visited = initSet(get_s(fsm));

	// add set of all state (init) on queue
	int* initS = (int*)malloc( (SIZE + 2 + get_s(fsm)) * sizeof(int) );

	initS[0] = 0;
	int j;

	for(j=1; j<get_s(fsm)+1; j++)
            initS[j] = 1;
            
	for(j=get_s(fsm)+1; j<SIZE + 2 + get_s(fsm); j++)
            initS[j] = 0;
            
	queue.push(initS);
		
	while (!queue.empty()) {

		//printf("size: %ld\n", queue.size());
		int* tmp = queue.top();
		//free(queue.top());
		queue.pop();

		/*
		printf("\ncurrent: %d ", tmp[0]);

		for (int tes = 0; tes<get_s(fsm); tes++ )
			printf("%d ", tmp[tes+1]);
		printf("\n");
   		*/
                // tmp = (h + lb) (etats : 1 0 1 0 0 ... 0) (taille de seq) (seq : 1 0 0 0 1 ...) 
		
		int* states = (int*)malloc( sizeof(int) * get_s(fsm) );	
		for (int i=0; i<get_s(fsm); i++) {
			states[i] = tmp[i+1];
			//printf("%d ", states[i]);
		}
		//printf("\n");
		
		//int* states = tmp+1;
		
        // size of current seq
		int seq_size = tmp[get_s(fsm)+1];
		
        int k, tmp_cpt = 0;
        for(k=0; k<get_s(fsm); k++)
            if (states[k] == 1)
                tmp_cpt++;

        // if singleton (1 and only one state) : SS
        if (tmp_cpt==1) {
            for(k=0; k<seq_size ;k++)
                res[k] = tmp[k+get_s(fsm)+2];
            clear_priorityQ(queue);
            freeSet(visited);
            return seq_size;
        }
		
		
		int added = add(visited, states);
		//printf("added: %d\n", added);
        if(added == 0) {
            //free(queue.top());
			//queue.pop();
			free(tmp);
			continue;
        }
        if(added == -1){
            printf("ERREUR CAPACITE ATTEINTE\n");
            break;
        }
               

		// c0 and c1 successors of current
		int *c0 = (int*)malloc( (get_s(fsm)) * sizeof(int) );
		int *c1 = (int*)malloc( (get_s(fsm)) * sizeof(int) );
		
		int j;
		for (j=0; j<get_s(fsm); j++) {
            c0[j] = 0;
            c1[j] = 0;
		}

		for (j=0; j<get_s(fsm); j++) {
            if (states[j]) {
                c0[ get_succ(fsm, j, 0) ] = 1;
                c1[ get_succ(fsm, j, 1) ] = 1;
            }
		}

		/*
		printf("succ: \n");
		for(j=0; j<get_s(fsm); j++)
			printf("%d %d\n", c0[j], c1[j]);
		printf("\n");
		*/

        int* succs[2] = {c0, c1};
                
        for(int symbol = 0; symbol < 2; symbol++){
            int* successor = succs[symbol];
            int inVisited = find(visited, successor);
            //printf("%d\n", inVisited);

            if(inVisited == 1)
                continue;
	
            //printf("addQ: ");

            int* addQ = (int*)malloc( (SIZE + 2 + get_s(fsm)) * sizeof(int) );
            addQ[0] = tmp[get_s(fsm) + 1] + maxMS(mergs, get_s(fsm), successor);
            //printf("%d ", addQ[0]);
            for (k=0; k<get_s(fsm); k++){
            	addQ[k+1] = successor[k];
            	//printf("%d ", addQ[k+1]);
            }
            addQ[get_s(fsm)+1] = seq_size+1;
            //printf("%d ", addQ[get_s(fsm)+1]);
            for (k=0; k<seq_size; k++){
            	addQ[k+get_s(fsm)+2] = tmp[k+get_s(fsm)+2];
            	//printf("%d ", addQ[k+get_s(fsm)+2]);
            }
            addQ[k+get_s(fsm)+2] = symbol;
            //printf("%d\n", addQ[k+get_s(fsm)+2]);

            queue.push(addQ);
            //printf("qtop %d\n", queue.top()[0]);
            
            //free(successor);
		}
		
		//free(states);
		free(c0);
		free(c1);
		
		//free(queue.top());
		//queue.pop();
		free(tmp);
		//printf("size: %ld\n",queue.size());	
	}

	return -1;
}


int test_time() {

	int i;
	//float total = 0;
	double ttl = 0;

	clock_t beg = clock();

	for (i=1; i<=50; i++) {
		    	

        char tmp[SIZE] = "./data/SS_50fsm_n30/fsm_n30_";
		char numb[10];
        sprintf(numb, "%d", i);
        strcat(tmp, numb);
		strcat( tmp, ".fsm" );
	    
        //const char* tmp = "./data/fsm_hss.fsm";
        printf("%s\n", tmp );

		FSM fsm = initFSM(tmp);
		//printFSM(fsm);

		int res[SIZE];

		MergSeq mergs = initMS(fsm);

		//time_t begin = time(NULL);
		clock_t b = clock();

		int taille = syncTree(fsm, res);
		//int taille = bestSearch(fsm, res, mergs);

		//time_t end = time(NULL);
		clock_t e = clock();

		//total += difftime(end,begin);
		ttl += (double)(e-b);

		//printf("dfa %d: length SS: %d time: %f s / %f ms\n", i, taille, difftime(end, begin), (double)(e-b));
		printf( "dfa %d: length SS: %d \t time: %f ms\n\n", i, taille, (double)(e-b) * 1000 / CLOCKS_PER_SEC);
		//printf( "dfa %d: length SS: %d\n", i, taille );
		

		freeMS(mergs, get_s(fsm));

		//int j;
		//printf("SS: ");
		//for (j=0; j<taille; j++){
		//	printf("%d", res[j]);
		//}
		//printf("\n\n");
		
		freeFSM(fsm);

	}

	//printf("\ntotal time: %f ms / average (50 dfa): %f ms\n", ttl, ttl/(i-1));
	//printf("total time: %f seconds / average (50 dfa): %f seconds\n", total, total/(i-1));

	clock_t end = clock();
	double total = (double)(end-beg)* 1000 / CLOCKS_PER_SEC;
	ttl = ttl *1000 / CLOCKS_PER_SEC;
	printf("total time: %f ms / average: %f ms (%d fsm)\n", total, total/(i-1), i-1);
	printf("time without MS time:\ntotal time: %f ms / average: %f ms (%d fsm)\n", ttl, ttl/(i-1), i-1);

	return 0;

}

// test bestSearch
int test_best() {

	char tmp[SIZE] = "./data/fsm_n20_1.fsm";
	//char tmp[SIZE] = "./data/fsm_hss.fsm";
    printf( "%s\n", tmp );

	FSM fsm = initFSM(tmp);
	printFSM(fsm);


	int res[SIZE];

	MergSeq mergs = initMS(fsm);

	printf("best: %d\n", bestSearch(fsm, res, mergs));

	freeMS(mergs, get_s(fsm));
	freeFSM(fsm);

	return 0;
}



// test MergSeq and maxMS
int test_allms() {

	//char tmp[SIZE] = "./data/fsm_hss.fsm";
	char tmp[SIZE] = "./data/fsm_n20_1.fsm";

    printf( "%s\n", tmp );

	FSM fsm = initFSM(tmp);
	printFSM(fsm);

	MergSeq mergs = initMS(fsm);

	int tmpS[5] = {1,1,1,1};
	//int tmpS[5] = {0,1,1,1};
	printf( "maxMS: %d\n", maxMS(mergs, get_s(fsm), tmpS) );

	freeMS(mergs, get_s(fsm));
	freeFSM(fsm);

	return 0;
}

int main() {

	test_time();
	
	//test_best();
	//test_allms();

	return 0;
}
