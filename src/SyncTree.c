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
#define DELIM_MS 7


// free queue and visited
void clear_Q (int ** q, int taille_q) {
    int i;
    for (i=0; i<taille_q; i++) {
            free(q[i]);
    }
    free(q);
}

// Synchronizing tree
int syncTree(FSM fsm, int* res) {

	int j;

	// queue
	int **queue = (int**)malloc( SIZE_Q * sizeof(int*) );
	int front_q = 0;
	int back_q = 1;

	// visited set of states
    Set visited = initSet(get_s(fsm), 0);

	if (!(queue))
		return -1;

	// add set of all state (init) on queue and visited
	queue[0] = (int*)malloc( (SIZE + 1 + get_s(fsm)) * sizeof(int) );

	for(j=0; j<get_s(fsm); j++)
            queue[0][j] = 1;
	for(j=get_s(fsm); j<SIZE + 1 + get_s(fsm); j++)
            queue[0][j] = 0;

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
		int seq_size = states[get_s(fsm)];
                    
        int k, tmp_cpt = 0;
        for(k=0; k<get_s(fsm); k++)
            if (states[k] == 1)
                tmp_cpt++;

        // if singleton (1 and only one state) : SS
        if (tmp_cpt==1) {
            for(k=0; k<seq_size ;k++)
                res[k] = states[k+get_s(fsm)+1];
            clear_Q(queue, back_q);
            freeSet(visited);
            return seq_size;
        }


		//printf("current: %d / ttl_q: %d \n", current, current + count_q+1);

		// c0 and c1 successors of current
		int *c0 = (int*)malloc((SIZE + 1 + get_s(fsm)) * sizeof(int));
		int *c1 = (int*)malloc((SIZE + 1 + get_s(fsm)) * sizeof(int));
		
		int j;

		for (j=0; j<SIZE + 1 + get_s(fsm); j++) {
            c0[j] = 0;
            c1[j] = 0;
		}

		for (j=0; j<get_s(fsm); j++) {
            if (states[j]) {
                c0[ get_succ(fsm, j, 0) ] = 1;
                c1[ get_succ(fsm, j, 1) ] = 1;
            }

		}
        int* succs[2] = {c0, c1};


        for(int symbol = 0; symbol < 2; symbol++){
            int* successor = succs[symbol];
            int inVisited = find(visited, successor);

            if(inVisited == 1){
                free(successor);
                continue;

            }
	
            queue[back_q] = successor;
            queue[back_q][get_s(fsm)] = seq_size+1;

            for(k=0; k<seq_size; k++) {
                queue[back_q][k+get_s(fsm)+1] = states[k+get_s(fsm)+1];
            }
            queue[back_q][k + get_s(fsm) + 1] = symbol;
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
   	Set visited = initSet(get_s(fsm), 0);

	// add set of all state (init) on queue
	int* initS = (int*)malloc( (SIZE + 2 + get_s(fsm)) * sizeof(int) );

	initS[0] = 0;
	int j;

	for(j=1; j<get_s(fsm)+1; j++)
            initS[j] = 1;
            
	for(j=get_s(fsm)+1; j<SIZE + 2 + get_s(fsm); j++)
            initS[j] = 0;
            
	queue.push(initS);
		
	//int cpt_lin = 0;
    //int cpt_sort = 0;

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
            //printf("linear: %d \nsorted: %d\n", cpt_lin, cpt_sort);
            return seq_size;
        }
		
		
		int added = add(visited, states);
		//printf("added: %d\n", added);
        if(added == 0) {
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

		int size_c0=0;
		int size_c1=0;

		for (j=0; j<get_s(fsm); j++) {
            if (states[j]) {
            	
            	if (c0[ get_succ(fsm, j, 0) ] == 0 ){
            		c0[ get_succ(fsm, j, 0) ] = 1;
            		size_c0++;
            	} 
            	if (c1[ get_succ(fsm, j, 1) ] == 0 ){
            		c1[ get_succ(fsm, j, 1) ] = 1;
            		size_c1++;
            	} 
            	
        		//c0[ get_succ(fsm, j, 0) ] = 1;
        		//c1[ get_succ(fsm, j, 1) ] = 1;

            }
		}

		/*
		printf("succ: \n");
		for(j=0; j<get_s(fsm); j++)
			printf("%d %d\n", c0[j], c1[j]);
		printf("\n");
		*/

        int* succs[2] = {c0, c1};
        int size_succs[2] = {size_c0, size_c1};
        
                
        for(int symbol = 0; symbol < 2; symbol++){
            int* successor = succs[symbol];
            int inVisited = find(visited, successor);
            //printf("%d\n", inVisited);

            if(inVisited == 1)
                continue;
	

            int* addQ = (int*)malloc( (SIZE + 2 + get_s(fsm)) * sizeof(int) );
            
            if (size_succs[symbol] < DELIM_MS){
            	addQ[0] = tmp[get_s(fsm) + 1] + maxMSLinear(mergs, get_s(fsm), successor);
            	//cpt_lin++;
            }
            else{
	            addQ[0] = tmp[get_s(fsm) + 1] + maxMSSorted(mergs, get_s(fsm), successor);
	            //cpt_sort++;
            }
            
            //printf("addQ: ");
        	//addQ[0] = tmp[get_s(fsm) + 1] + maxMSSorted(mergs, get_s(fsm), successor);

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
		
		free(tmp);
		//printf("size: %ld\n",queue.size());	
	}

	return -1;
}


// best Sync tree with heap and storage
int bestSearchUnique(FSM fsm, int *res, MergSeq mergs){

    std::priority_queue<int*, std::vector<int*>, CompareQ> queue;
    // visited set of states
    Set storage = initSet(get_s(fsm), 1);

    // add set of all state (init) on queue
    int* initS = (int*)malloc( (2 + get_s(fsm)) * sizeof(int) );

    initS[0] = 0;
    int j;

    for(j=1; j<get_s(fsm)+1; j++)
        initS[j] = 1;
        
    for(j=get_s(fsm)+1; j<2 + get_s(fsm); j++)
        initS[j] = 0;
        
    queue.push(initS);
            
    while (!queue.empty()) {

        int* states = queue.top();
        queue.pop();

        // size of current seq
        int seq_size = states[get_s(fsm)+1];
		
        int k, tmp_cpt = 0;
        for(k=0; k<get_s(fsm); k++){
            if (states[k + 1] == 1)
                tmp_cpt++;
            if (tmp_cpt == 3)
                break;
        }

        // if singleton (1 and only one state) : SS
        if (tmp_cpt==1) {
            clear_priorityQ(queue);
            freeSet(storage);
            return seq_size;
        }
        
        // if couple (2 states) : SS = current + MS
        if (tmp_cpt==2) {
            clear_priorityQ(queue);
            freeSet(storage);
            return states[0];
        }
     

        // c0 and c1 successors of current
        int *c0 = (int*)malloc((2 + get_s(fsm)) * sizeof(int) );
        int *c1 = (int*)malloc((2 + get_s(fsm)) * sizeof(int) );
        
        int j;
        for (j=0; j < get_s(fsm); j++) {
            c0[j + 1] = 0;
            c1[j + 1] = 0;
        }

        for (j=0; j< get_s(fsm); j++) {
            if (states[j + 1]) {
                c0[ get_succ(fsm, j, 0) + 1] = 1;
                c1[ get_succ(fsm, j, 1) + 1] = 1;
            }
        }

        int* succs[2] = {c0, c1};
                
        for(int symbol = 0; symbol < 2; symbol++){
            int* successor = succs[symbol];
            int* succInVisited = addOrReturn(storage, successor);

            if(succInVisited == successor){
                successor[0] = seq_size + 1 + maxMSSorted(mergs, get_s(fsm), successor + 1);
                successor[get_s(fsm) + 1] = seq_size + 1;
            }
            else{
                int seq_size_n = succInVisited[get_s(fsm) + 1];
                if(seq_size + 1 < seq_size_n){
                    // UPDATE DANS QUEUE !!!
                    succInVisited[0] -= seq_size_n - (seq_size + 1);
                    succInVisited[get_s(fsm) + 1] -= seq_size_n - (seq_size + 1);
                }
                free(successor);
            }
            
        }
    }

    return -1;
}


int test_time() {

	int i;
	//float total = 0;
	double ttl = 0;

	clock_t beg = clock();

	int noMS= 0;

	for (i=1; i<=50; i++) {
		    	

        char tmp[SIZE] = "./data/SS_50fsm_n30/fsm_n30_";
		char numb[10];
        sprintf(numb, "%d", i);
        strcat(tmp, numb);
		strcat( tmp, ".fsm" );
	    
        //const char* tmp = "./data/fsm_hss.fsm";
        ////printf("%s\n", tmp );

		FSM fsm = initFSM(tmp);
		//printFSM(fsm);


		int res[SIZE];

		MergSeq mergs = initMS(fsm);

		if (mergs == NULL){
			noMS++;
			printf("dfa %d: no SS\n\n", i);
			continue;
		}

		//time_t begin = time(NULL);
		clock_t b = clock();

		//int taille = syncTree(fsm, res);

		int taille = bestSearch(fsm, res, mergs);

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

	i -= noMS;

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

	int tmpS[20] = {0,1,0,1, 1};
	//int tmpS[] = {0,1,1,1};
	printf( "maxMS: %d\n", maxMSLinear(mergs, get_s(fsm), tmpS) );
	printf( "maxMS: %d\n", maxMSSorted(mergs, get_s(fsm), tmpS) );

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
