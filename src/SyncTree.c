#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "data_structures/set.h"

#define NB_SUCC 2 // nb of successors
#define TAILLE_MAX 20 // readfile
#define SIZE 100 // filename and res
#define SIZE_Q 1000000

#define INT_SIZE (int)(sizeof(int)*8-1)

typedef struct {
	int type;
	int s;
	int i;
	int p;
	int succ[][NB_SUCC];
} FSM;


// create FSM with s states, i inputs and the transitions trans
FSM * setFSM(int s, int i, int** trans) {

	FSM* fsm = (FSM*)malloc( sizeof(FSM) + 2*sizeof(int) * s );

	fsm->type = 0;
	fsm->s = s;
	fsm->i = i;
	fsm->p = s*i;

	int j;
	for (j=0; j<s; j++)
	{
		fsm->succ[j][0] = trans[j][0];
		fsm->succ[j][1] = trans[j][1];
	}

	return fsm;
}


// print a FSM
void printFSM(FSM *fsm) {

	printf("\tPrint FSM\n");

	printf("F %d\n", fsm->type);
	printf("s %d\n", fsm->s);
	printf("i %d\n", fsm->i);
	printf("p %d\n", fsm->p);

	int j;
	for(j=0; j<fsm->s; j++)
	{
		printf("%d 0 %d 0\n", j,fsm->succ[j][0]);
		printf("%d 1 %d 0\n", j,fsm->succ[j][1]);
	}

	printf("\tEnd Print FSM\n\n");

}

// free queue and visited
void clear_QV (int ** q, int taille_q) {
    int i;
    for (i=0; i<taille_q; i++) {
            free(q[i]);
    }
    free(q);
}

int* initNode(FSM * fsm){
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

    int size = fsm->s / INT_SIZE + 1; // Number of integers needed to encode all the bits
    int* state = (int*)malloc( (SIZE + 1 + size) * sizeof(int) );

    // The initial states contains
    for(int j = 0; j < size - 1; j++){
        state[j] = ((unsigned int)(~0)) >> 1; // Only ones except the first
    }
    int last_size = fsm->s % INT_SIZE;
    state[size - 1] = ((unsigned int)(~0)) >> (INT_SIZE - last_size + 1); // 2^(INT_SIZE - lastsize) - 1

    // Init the size of the sequence following the states and the sequence itself with 0
    for(int j=size; j<SIZE + 1 + size; j++)
        state[j] = 0;
    return state;
}

int* successor(FSM * fsm, int* node, int input){
    int size = fsm->s / INT_SIZE + 1;
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
                int succ = fsm->succ[index][input];
                neigh[succ / INT_SIZE] |= (1 << (succ % INT_SIZE));
            }
            index += 1;
            x /= 2;
        }
    }

    return neigh;
}


// free FSM
void freeFSM(FSM *fsm) {
	free(fsm);
}


// read FSM from filename
FSM * readFSM(const char * filename) {

	FILE * file;
	if ((file = fopen(filename, "r")) == NULL) {
		printf("Error opening file, exit.\n");
		exit(1);
	}

	char chaine[TAILLE_MAX] = "";

	while (fgets(chaine, TAILLE_MAX, file)[0] != 's') {}

	//int s = chaine[2] - '0';
	char *spt = strtok(chaine, " ");
	spt = strtok(NULL, " ");
	int s = atoi(spt);

	int i = fgets(chaine, TAILLE_MAX, file)[2] - '0';
	
	//printf("%d %d\n", s,i);

	while (fgets(chaine, TAILLE_MAX, file)[0] != 'p') {}
	
	int ** tr = (int**)malloc( s * sizeof(int[2]) );

	int j;
	for(j=0; j<s; j++) {
		tr[j] = (int*)malloc(2 * sizeof(int));
	}

	for (j=0; j<s; j++) {
		// inp 0
		fgets(chaine, TAILLE_MAX, file);
		//printf("%s", chaine);
		//tr[j][0] = chaine[4] - '0';
		spt = strtok(chaine, " ");
		spt = strtok(NULL, " ");
		spt = strtok(NULL, " ");
		tr[j][0] = atoi(spt);

		// inp 1	
		fgets(chaine, TAILLE_MAX, file);		
		//printf("%s\n", chaine);
		//tr[j][1] = chaine[4] - '0';
		spt = strtok(chaine, " ");
		spt = strtok(NULL, " ");
		spt = strtok(NULL, " ");
		tr[j][1] = atoi(spt);
	
	}	
	fclose(file);

	FSM * fsm = setFSM(s,i,tr);

	for (j=0; j<s; j++) {
		free(tr[j]);
	}
	free(tr);

	return fsm;
}


// Synchronizing tree
int syncTree(FSM *fsm, int* res) {


	// queue
	int **queue = (int**)malloc( SIZE_Q * sizeof(int*) );
	int front_q = 0;
	int back_q = 1;
        int encoding_size = fsm->s / INT_SIZE + 1;

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
                    clear_QV(queue, back_q);
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
                        queue[back_q][encoding_size + 1 + k] = states[encoding_size + 1 + k];
                    }
                    queue[back_q][encoding_size + 1 + k] = symbol;
                    back_q++;



                }
	}


	return -1;

}



int main() {

	int i;
	//float total = 0, ttl = 0;
	clock_t beg = clock();

	for (i=1; i<=10; i++) {
	    	
            char tmp[SIZE] = "data/SS_50fsm_n50/fsm_n50_";
		char numb[10];
                sprintf(numb, "%d", i);
                strcat(tmp, numb);
		strcat( tmp, ".fsm" );
	    
                printf("%s\n", tmp );

		FSM * fsm = readFSM(tmp);
                //printFSM(fsm);

		//time_t begin = time(NULL);
		clock_t b = clock();


		int res[SIZE];
		//syncTree(fsm, res);
		int taille = syncTree(fsm, res);

		//time_t end = time(NULL);
		clock_t e = clock();

		//total += difftime(end,begin);
		//ttl += (double)(e-b);

		//printf("dfa %d: length SS: %d time: %f s / %f ms\n", i, taille, difftime(end, begin), (double)(e-b));
		printf( "dfa %d: length SS: %d \t time: %f ms\n\n", i, taille, (double)(e-b) * 1000 / CLOCKS_PER_SEC);
		//printf( "dfa %d: length SS: %d\n", i, taille );
		
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
	double ttl = (double)(end-beg)* 1000 / CLOCKS_PER_SEC;
	printf("total time: %f ms / average: %f ms (%d fsm)", ttl, ttl/(i-1), i-1);

	return 0;

}
