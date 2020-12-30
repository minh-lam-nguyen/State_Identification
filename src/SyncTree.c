#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <queue>
#include "data_structures/set.h"
#include "data_structures/set_inc.h"
#include "data_structures/fsm.h"
#include "data_structures/mergseq.h"
#include "data_structures/tergseq.h"
#include "data_structures/quadseq.h"
#include "data_structures/shp.h"
#include "data_structures/heap.h"


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

            //clear_Q(queue, back_q);
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
int bestSearch(FSM fsm, MergSeq mergs){

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
int bestSearchUnique(FSM fsm, MergSeq mergs){

    Heap queue = initHeap(0, get_s(fsm)+2);

    // visited set of states
    Set storage = initSet(get_s(fsm), 1);

    // add set of all state (init) on queue
    int* initS = (int*)malloc( (3 + get_s(fsm)) * sizeof(int) );

    initS[0] = 0;
    int j;

    for(j=1; j<get_s(fsm)+1; j++)
        initS[j] = 1;

    initS[get_s(fsm)+1] =0;

    insertToHeap(queue, initS);

    int nbExplore = 0;
    while (sizeOfHeap(queue) > 0) {
        nbExplore++;
        int* states = popMinimum(queue);

        /*
        printf("\nstates: %d ", states[0]);
        for (int ttt = 1; ttt<get_s(fsm)+1; ttt++){
        	printf("%d ", states[ttt]);
        }
        printf("\n");
        printf("size: %d\n", states[get_s(fsm) +1]);
		*/

        // size of current seq
        int seq_size = states[get_s(fsm)+1];
		
        int k, tmp_cpt = 0;
        for(k=0; k<get_s(fsm); k++){
            if (states[k + 1] == 1)
                tmp_cpt++;
            if (tmp_cpt == 3)
                break;
        }

        //printf("%d\n", tmp_cpt);

        // if singleton (1 and only one state) : SS
        if (tmp_cpt==1) {
            printf(">> %d |", nbExplore);
            freeHeap(queue);
            freeSet(storage);
            return seq_size;
        }
        
        // if couple (2 states) : SS = current + MS
        if (tmp_cpt==2) {
            printf(">> %d |", nbExplore);
        	int res = states[0];
            freeHeap(queue);
            freeSet(storage);
            //return states[0];
            return res;
        }

        // c0 and c1 successors of current
        int *c0 = (int*)malloc((3 + get_s(fsm)) * sizeof(int) );
        int *c1 = (int*)malloc((3 + get_s(fsm)) * sizeof(int) );
        
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
                //printf("key succ: %d\n", successor[0]);
                successor[get_s(fsm) + 1] = seq_size + 1;
                //printf("size succ: %d\n", successor[get_s(fsm) + 1]);
                insertToHeap(queue, successor);
            }
            else{
                int seq_size_n = succInVisited[get_s(fsm) + 1];
                //printf("old, new: %d %d\n", seq_size, seq_size_n);
                if(seq_size + 1 < seq_size_n){
                    // UPDATE DANS QUEUE !!!
                    decreaseKey(queue, succInVisited, succInVisited[0] - seq_size_n + (seq_size + 1));
                    succInVisited[get_s(fsm) + 1] -= seq_size_n - (seq_size + 1);
                }
                //free(successor);
            }
            
        }

    }

    return -1;
}


// best Sync tree with heap and storage and MS of triplets
int bestSearchUniqueQ(FSM fsm, QuadSeq quads){

    Heap queue = initHeap(0, get_s(fsm)+2);

    // visited set of states
    Set storage = initSet(get_s(fsm), 1);

    // add set of all state (init) on queue
    int* initS = (int*)malloc( (3 + get_s(fsm)) * sizeof(int) );

    initS[0] = 0;
    int j;

    for(j=1; j<get_s(fsm)+1; j++)
        initS[j] = 1;

    initS[get_s(fsm)+1] =0;

    insertToHeap(queue, initS);

    int nbExplore = 0;
    while (sizeOfHeap(queue) > 0) {
        nbExplore++;
        int* states = popMinimum(queue);

        /*
        printf("\nstates: %d ", states[0]);
        for (int ttt = 1; ttt<get_s(fsm)+1; ttt++){
        	printf("%d ", states[ttt]);
        }
        printf("\n");
        printf("size: %d\n", states[get_s(fsm) +1]);
		*/

        // size of current seq
        int seq_size = states[get_s(fsm)+1];
		
        int k, tmp_cpt = 0;
        for(k=0; k<get_s(fsm); k++){
            if (states[k + 1] == 1)
                tmp_cpt++;
            if (tmp_cpt == 4)
                break;
        }

        //printf("%d\n", tmp_cpt);

        // if singleton (1 and only one state) : SS
        if (tmp_cpt==1) {
            printf(">> %d |", nbExplore);
            freeHeap(queue);
            freeSet(storage);
            return seq_size;
        }
        
        // if couple or triple (2 or 3 states) : SS = current + MS
        if (tmp_cpt <= 3) {
            printf(">> %d |", nbExplore);
            int res = states[0];
            freeHeap(queue);
            freeSet(storage);
            //return states[0];
            return res;
        }

        // c0 and c1 successors of current
        int *c0 = (int*)malloc((3 + get_s(fsm)) * sizeof(int) );
        int *c1 = (int*)malloc((3 + get_s(fsm)) * sizeof(int) );
        
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
                successor[0] = seq_size + 1 + maxQSSorted(quads, get_s(fsm), successor + 1);
                //printf("key succ: %d\n", successor[0]);
                successor[get_s(fsm) + 1] = seq_size + 1;
                //printf("size succ: %d\n", successor[get_s(fsm) + 1]);
                insertToHeap(queue, successor);
            }
            else{
                int seq_size_n = succInVisited[get_s(fsm) + 1];
                //printf("old, new: %d %d\n", seq_size, seq_size_n);
                if(seq_size + 1 < seq_size_n){
                    // UPDATE DANS QUEUE !!!
                    decreaseKey(queue, succInVisited, succInVisited[0] - seq_size_n + (seq_size + 1));
                    succInVisited[get_s(fsm) + 1] -= seq_size_n - (seq_size + 1);
                }
                //free(successor);
            }
            
        }

    }

    return -1;
}


// best Sync tree with heap and storage and MS of triplets
int bestSearchUniqueT(FSM fsm, TergSeq tergs){

    Heap queue = initHeap(0, get_s(fsm)+2);

    // visited set of states
    Set storage = initSet(get_s(fsm), 1);

    // add set of all state (init) on queue
    int* initS = (int*)malloc( (3 + get_s(fsm)) * sizeof(int) );

    initS[0] = 0;
    int j;

    for(j=1; j<get_s(fsm)+1; j++)
        initS[j] = 1;

    initS[get_s(fsm)+1] =0;

    insertToHeap(queue, initS);

    int nbExplore = 0;
    while (sizeOfHeap(queue) > 0) {
        nbExplore++;
        int* states = popMinimum(queue);

        /*
        printf("\nstates: %d ", states[0]);
        for (int ttt = 1; ttt<get_s(fsm)+1; ttt++){
        	printf("%d ", states[ttt]);
        }
        printf("\n");
        printf("size: %d\n", states[get_s(fsm) +1]);
		*/

        // size of current seq
        int seq_size = states[get_s(fsm)+1];
		
        int k, tmp_cpt = 0;
        for(k=0; k<get_s(fsm); k++){
            if (states[k + 1] == 1)
                tmp_cpt++;
            if (tmp_cpt == 4)
                break;
        }

        //printf("%d\n", tmp_cpt);

        // if singleton (1 and only one state) : SS
        if (tmp_cpt==1) {
            printf(">> %d |", nbExplore);
            freeHeap(queue);
            freeSet(storage);
            return seq_size;
        }
        
        // if couple or triple (2 or 3 states) : SS = current + MS
        if (tmp_cpt <= 3) {
            printf(">> %d |", nbExplore);
            int res = states[0];
            freeHeap(queue);
            freeSet(storage);
            //return states[0];
            return res;
        }

        // c0 and c1 successors of current
        int *c0 = (int*)malloc((3 + get_s(fsm)) * sizeof(int) );
        int *c1 = (int*)malloc((3 + get_s(fsm)) * sizeof(int) );
        
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
                successor[0] = seq_size + 1 + maxTSSorted(tergs, get_s(fsm), successor + 1);
                //printf("key succ: %d\n", successor[0]);
                successor[get_s(fsm) + 1] = seq_size + 1;
                //printf("size succ: %d\n", successor[get_s(fsm) + 1]);
                insertToHeap(queue, successor);
            }
            else{
                int seq_size_n = succInVisited[get_s(fsm) + 1];
                //printf("old, new: %d %d\n", seq_size, seq_size_n);
                if(seq_size + 1 < seq_size_n){
                    // UPDATE DANS QUEUE !!!
                    decreaseKey(queue, succInVisited, succInVisited[0] - seq_size_n + (seq_size + 1));
                    succInVisited[get_s(fsm) + 1] -= seq_size_n - (seq_size + 1);
                }
                //free(successor);
            }
            
        }

    }

    return -1;
}


// backward best Sync tree with heap and storage
int bwBestSearchUnique(FSM fsm, Shp shps){

    int size = get_s(fsm);
    Heap queue = initHeap(0, size+2);

    // visited set of states
    SetInc storage = initSetInc(size, 1);

    // add set of singletons on queue

    for(int i = 0; i < size; i++){
        int* initS = (int*)malloc( (3 + size) * sizeof(int) );
        for(int j=1; j<size + 1; j++)
            initS[j] = (j - 1 == i) ? 1 : 0;
        int b = maxShpSorted(shps, size, initS + 1);
        initS[size + 2] = -1;
        if(b == -1){
            free(initS);
            continue;
        }

        initS[size + 1] = 0;
        initS[0] = b;
        insertToHeap(queue, initS);
        addInc(storage, initS);
    }

    int nbExplore = 0;
    while (sizeOfHeap(queue) > 0) {
        nbExplore++;
        int* states = popMinimum(queue);

        // size of current seq
        int seq_size = states[get_s(fsm)+1];
		
        int tmp_cpt = 1;
        for(int k=0; k<get_s(fsm); k++){
            if (states[k + 1] == 1)
                continue;
            tmp_cpt = 0;
            break;
        }

        // if tmp_cpt is 1, states contain all states
        if (tmp_cpt==1) {
            //printf(">> %d |", nbExplore);
            freeHeap(queue);
            freeSetInc(storage);
            return seq_size;
        }
        
        // c0 and c1 preimages of current
        int *c0 = (int*)malloc((3 + size) * sizeof(int) );
        int *c1 = (int*)malloc((3 + size) * sizeof(int) );
        
        for (int j=0; j < size; j++) {
                int s0 = get_succ(fsm, j, 0);
                int s1 = get_succ(fsm, j, 1);
                
                c0[j + 1] = states[s0 + 1];
                c1[j + 1] = states[s1 + 1];
        }
        /*c0[0] = 0;
        c1[0] = 0;
        c0[size + 1] = 0;
        c1[size + 1] = 0;
        c0[size + 2] = 0;
        c1[size + 2] = 0;
        */

        int* preds[2] = {c0, c1};

        for(int symbol = 0; symbol < 2; symbol++){
            int* predecessor = preds[symbol];
            int b = maxShpSorted(shps, size, predecessor + 1);
        
            if(b == -1)
                continue;

            int* predInVisited = addIncOrReturn(storage, predecessor);
            if(predInVisited == predecessor){
                predecessor[0] = seq_size + 1 + b;
                predecessor[size + 1] = seq_size + 1;
                predecessor[size + 2] = -1;
                insertToHeap(queue, predecessor);
            }
            else{
                int seq_size_n = predInVisited[size + 1];
                if(seq_size + 1 < seq_size_n){
                    // UPDATE DANS QUEUE !!!
                    decreaseKey(queue, predInVisited, predInVisited[0] - (seq_size_n - (seq_size + 1)));
                    predInVisited[size + 1] -= seq_size_n - (seq_size + 1);
                }
                
                free(predecessor);
            }
            
        }

    }

    return -1;
}


int test_time() {

	int i;

        // Indique quels algos sont testés
        // Si 1, l'algo est testé sinon il ne l'est pas.
        // 0 : bestSearch, avec set visited, borne MS
        // 1 : bestSearchUnique, avec dict storage, borne MS
        // 2 : bestSearchUniqueT, avec dict storage, borne TS
        // 3 : bwBestSearchUnique, avec dict storage, borne TS
        int algs[] = {0, 1, 1, 1, 0};

	clock_t begPrecomputing, endPrecomputing, b, e;
	int noMS= 0;
        double precTimes[4];
        double times[4];

        for(int algIndex = 0; algIndex < 4; algIndex++){
            if(algs[algIndex] == 0)
                continue;

            for (i=3; i<=50; i++) {
                                    
                begPrecomputing = clock();
                char tmp[SIZE] = "./data/SS_50fsm_n70/fsm";
                //char tmp[SIZE] = "./data/fsm_n20_";
                char numb[10];
                sprintf(numb, "%d", i);
                strcat(tmp, numb);
                strcat( tmp, ".fsm" );

                FSM fsm = initFSM(tmp);
                //printFSM(fsm);

                TergSeq tergs;
                MergSeq mergs;
                QuadSeq quads;
                Shp shps;

                switch(algIndex){
                    case 0:
                    case 1:
                        mergs = initMS(fsm);
                        if (mergs == NULL){
                                noMS++;
                                printf("%d X ", i);
                                continue;
                        }
                        break;
                    case 2:
                        tergs = initTS(fsm);

                        if (tergs == NULL){
                                noMS++;
                                printf("%d X ", i);
                                continue;
                        }
                        break;
                    case 3:
                        quads = initQS(fsm);

                        if (quads == NULL){
                                noMS++;
                                printf("%d X ", i);
                                continue;
                        }
                        break;
                    case 4:
                        mergs = initMS(fsm);
                        shps = initSHP(fsm);

                        if (mergs == NULL){
                                noMS++;
                                printf("%d X ", i);
                                continue;
                        }
                        break;
                }

                endPrecomputing = clock();
                precTimes[algIndex] += (double)(endPrecomputing - begPrecomputing);


                int taille;
                b = clock();
                switch(algIndex){
                    case 0:
                        taille = bestSearch(fsm, mergs);
                        break;
                    case 1:
                        taille = bestSearchUnique(fsm, mergs);
                        break;
                    case 2:
                        taille = bestSearchUniqueT(fsm, tergs);
                        break;
                    case 3:
                        taille = bestSearchUniqueQ(fsm, quads);
                        break;
                    case 4:
                        taille = bwBestSearchUnique(fsm, shps);
                        break;
                }

                e = clock();
                times[algIndex] += (double)(e-b);

                printf("%d ", i);

                printf("%d | ", taille);
                //printf( "dfa %d: length SS: %d \t time: %f ms\n\n", i, taille, (double)(e-b) * 1000 / CLOCKS_PER_SEC);
                
                switch(algIndex){
                    case 0:
                    case 1:
                        freeMS(mergs, get_s(fsm));
                        break;
                    case 2:
                        freeTS(tergs, get_s(fsm));
                        break;
                    case 3:
                        freeQS(quads, get_s(fsm));
                        break;
                    case 4:
                        freeMS(mergs, get_s(fsm));
                        freeSHP(shps, get_s(fsm));
                        break;
                }
                
                freeFSM(fsm);

            }
            printf("\n");
            i -= noMS;
            double total = (precTimes[algIndex] + times[algIndex]) * 1000 / CLOCKS_PER_SEC;
            double time = times[algIndex] * 1000 / CLOCKS_PER_SEC;
            printf("total time: %f ms / average: %f ms (%d fsm)\n", total, total/(i-1), i-1);
            printf("time without MS time:\ntotal time: %f ms / average: %f ms (%d fsm)\n", time, time / (i-1), i-1);
        }


	return 0;

}

// test bestSearch
int test_bw() {

	char tmp[SIZE] = "./data/cerny_fsm/cerny_n5.fsm";
        printf( "%s\n", tmp );

	FSM fsm = initFSM(tmp);
	printFSM(fsm);

	Shp shps = initSHP(fsm);

	printf("best: %d\n", bwBestSearchUnique(fsm, shps));

	freeSHP(shps, get_s(fsm));
	freeFSM(fsm);

	return 0;
}

int main() {

        test_time();
	//test_bw();

	return 0;
}
