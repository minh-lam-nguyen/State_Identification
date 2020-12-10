#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <queue>
#include "data_structures/set.h"

#define NB_SUCC 2 // nb of successors
#define TAILLE_MAX 50 // readfile and filename
#define SIZE 150 // res
#define SIZE_Q 10000000

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

// compute min merging seq for (s1,s2), with s1 <= s2
// return length of merging seq for (s1, s2), -1 if no merging (or error)
int mergSeq(FSM *fsm, int s1, int s2){

	// elem of queue: ( state1, state2, length ) 
	int **queue = (int**)malloc( SIZE_Q * sizeof(int*) );
	if (!queue)
		return -1;
	int front_q = 0;
	int back_q = 1;

	queue[0] = (int*)malloc( 3 * sizeof(int) );
	queue[0][0] = s1;
	queue[0][1] = s2;
	queue[0][2] = 0;

	// elem of visited: ( state1, state2 )
	Set visited = initSet(2);

	while (front_q != back_q) {

		int current = front_q;
		front_q++;
		//printf("\ncurrent/f/b: %d %d %d\n", current, front_q, back_q);

		int *tmp = queue[current];
		//printf("tmp: %d %d %d\n", tmp[0], tmp[1], tmp[2]);

		// s1 and s2 merged 
		if (tmp[0] == tmp[1]){
			int res = tmp[2];
			clear_QV(queue, back_q);
            		freeSet(visited);
			return res;
		}

		// if (j,i) and i<j, always add (i,j) to visited
		if (tmp[0] > tmp[1]){
			int temp = tmp[0];
			tmp[0] = tmp[1];
			tmp[1] = temp;
		}

		//int vis[2] = {tmp[0], tmp[1]};
		int *vis = (int*)malloc( 2 * sizeof(int) );
		vis[0] = tmp[0];
		vis[1] = tmp[1];
		//printf("vis: %d %d\n", vis[0], vis[1]);

		//printf("find: %d\n", find(visited, vis));
		int added = add(visited, vis);
		//printf("added: %d\n", added);

		// already in visited
        if(added == 0) 
            continue;
        // visited too big
        if(added == -1){ 
            printf("ERREUR CAPACITE ATTEINTE\n");
            break;
        }

        // add successors to queue
        for (int i=0; i<fsm->i; i++){
        	//printf("add to q: %d %d %d\n", fsm->succ[tmp[0]][i], fsm->succ[tmp[1]][i], tmp[2]+1);
	        queue[back_q] = (int*)malloc( 3 * sizeof(int) );
        	queue[back_q][0] = fsm->succ[tmp[0]][i];
        	queue[back_q][1] = fsm->succ[tmp[1]][i];
        	queue[back_q][2] = tmp[2]+1;
        	back_q++;
        }

	}

	clear_QV(queue, back_q);
    	freeSet(visited);
	return -1;
}

// compute max merging seq (lower bound) for states s
int maxMS(FSM *fsm, int *s, int** ms){

	int res = -1;
	for (int i=0; i<fsm->s; i++) {
                if(!s[i])
                    continue;
		for (int j=i+1; j<fsm->s; j++) {
			if (!s[j])
                            continue;
                        int tmp = ms[i][j] ;
                        if (tmp > res)
                                res = tmp;	
		}
	}

	return res;
}

// Synchronizing tree
int syncTree(FSM *fsm, int* res) {

	int j;

	// queue
	int **queue = (int**)malloc( SIZE_Q * sizeof(int*) );
	int front_q = 0;
	int back_q = 1;

	// visited set of states
    Set visited = initSet(fsm->s);

	if (!(queue))
		return -1;

	// add set of all state (init) on queue and visited
	queue[0] = (int*)malloc( (SIZE + 1 + fsm->s) * sizeof(int) );

	for(j=0; j<fsm->s; j++)
            queue[0][j] = 1;
	for(j=fsm->s; j<SIZE + 1 + fsm->s; j++)
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
		int seq_size = states[fsm->s];
                    
                int k, tmp_cpt = 0;
                for(k=0; k<fsm->s; k++)
                    if (states[k] == 1)
                        tmp_cpt++;

                // if singleton (1 and only one state) : SS
                if (tmp_cpt==1) {
                    for(k=0; k<seq_size ;k++)
                        res[k] = states[k+fsm->s+1];
                    clear_QV(queue, back_q);
                    freeSet(visited);
                    return seq_size;
                }

		//printf("current: %d / ttl_q: %d \n", current, current + count_q+1);

		// c0 and c1 successors of current
		int *c0 = (int*)malloc((SIZE + 1 + fsm->s) * sizeof(int));
		int *c1 = (int*)malloc((SIZE + 1 + fsm->s) * sizeof(int));
		
		int j;
		for (j=0; j<SIZE + 1 + fsm->s; j++) {
                    c0[j] = 0;
                    c1[j] = 0;
		}

		for (j=0; j<fsm->s; j++) {
                    if (states[j]) {
                        c0[ fsm->succ[j][0] ] = 1;
                        c1[ fsm->succ[j][1] ] = 1;
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
                    queue[back_q][fsm->s] = seq_size+1;

                    for(k=0; k<seq_size; k++) {
                        queue[back_q][k+fsm->s+1] = states[k+fsm->s+1];
                    }
                    queue[back_q][k + fsm->s + 1] = symbol;
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
int bestSearch(FSM *fsm, int *res){

	std::priority_queue<int*, std::vector<int*>, CompareQ> queue;

/*
	// test priority queue
	int test1[5] = {4,0,1,0,1};
	int test2[5] = {3,0,1,0,1};
	printf("empty: %d\n", queue.empty());
	queue.push(test1);
	queue.push(test2);
	printf("empty: %d\n", queue.empty());
	printf("%d\n", queue.top()[0]);
	queue.pop();
	printf("%d\n", queue.top()[0]);
	queue.pop();
	printf("empty: %d\n", queue.empty());
*/

	// visited set of states
   	Set visited = initSet(fsm->s);

	// add set of all state (init) on queue
	int* initS = (int*)malloc( (SIZE + 2 + fsm->s) * sizeof(int) );

	initS[0] = 0;
	int j;
	for(j=1; j<fsm->s+1; j++)
            initS[j] = 1;
            
	for(j=fsm->s+1; j<SIZE + 2 + fsm->s; j++)
            initS[j] = 0;
            

	queue.push(initS);
	
	//printf("%d\n", queue.top()[0]);
	
	
	while (!queue.empty()) {

		int* tmp = queue.top();
		//free(queue.top());
		queue.pop();
		//printf("\ncurrent: %d\n", tmp[0]);
   
                // tmp = (h + lb) (etats : 1 0 1 0 0 ... 0) (taille de seq) (seq : 1 0 0 0 1 ...) 

		/*int* states = (int*)malloc( sizeof(int) * fsm->s );	
		for (int i=0; i<fsm->s; i++) {
			states[i] = tmp[i+1];
			//printf("%d ", states[i]);
		}
		//printf("\n");
		*/

		int* states = tmp+1;
		
                // size of current seq
		int seq_size = tmp[fsm->s+1];
		
                int k, tmp_cpt = 0;
                for(k=0; k<fsm->s; k++)
                    if (states[k] == 1)
                        tmp_cpt++;

                // if singleton (1 and only one state) : SS
                if (tmp_cpt==1) {
                    for(k=0; k<seq_size ;k++)
                        res[k] = tmp[k+fsm->s+2];
//                    clear_QV(queue, back_q);
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
		int *c0 = (int*)malloc( (fsm->s) * sizeof(int) );
		int *c1 = (int*)malloc( (fsm->s) * sizeof(int) );
		
		int j;
		for (j=0; j<fsm->s; j++) {
                    c0[j] = 0;
                    c1[j] = 0;
		}

		for (j=0; j<fsm->s; j++) {
                    if (states[j]) {
                        c0[ fsm->succ[j][0] ] = 1;
                        c1[ fsm->succ[j][1] ] = 1;
                    }
		}
                int* succs[2] = {c0, c1};
                
                for(int symbol = 0; symbol < 2; symbol++){
                    int* successor = succs[symbol];
                    int inVisited = find(visited, successor);

                    if(inVisited == 1){
                        //free(successor);
                        continue;

                    }
			
                    //printf("add: ");
                    int* addQ = (int*)malloc( (SIZE + 2 + fsm->s) * sizeof(int) );
                    addQ[0] = tmp[fsm->s + 1] + maxMS(fsm, successor);
                    //printf("%d ", addQ[0]);
                    for (k=0; k<fsm->s; k++){
                    	addQ[k+1] = successor[k];
                    	//printf("%d ", addQ[k+1]);
                    	}
                    addQ[fsm->s+1] = seq_size+1;
                    //printf("%d ", addQ[fsm->s+1]);
                    for (k=0; k<seq_size; k++){
                    	addQ[k+fsm->s+2] = tmp[k+fsm->s+2];
                    	//printf("%d ", addQ[k+fsm->s+2]);
                    	}
                    addQ[k+fsm->s+2] = symbol;
                    //printf("%d\n", addQ[k+fsm->s+2]);
                    queue.push(addQ);
                    //printf("qtop %d\n", queue.top()[0]);
                    
                    //free(successor);
		}
		
		free(states);
		free(c0);
		free(c1);
		
		//free(queue.top());
		//queue.pop();
		free(tmp);
		//printf("size: %ld\n",queue.size());	
	}

	return -1;
}


int test_ST() {

	int i;
	//float total = 0, ttl = 0;
	clock_t beg = clock();

	for (i=1; i<=50; i++) {
		    	
        char tmp[TAILLE_MAX] = "./data/SS_50fsm_n20/fsm_n20_";
		char numb[10];
        sprintf(numb, "%d", i);
        strcat(tmp, numb);
		strcat( tmp, ".fsm" );

	    
        //const char* tmp = "../data/fsm_hss.fsm";
        ////printf("%s\n", tmp );

		FSM * fsm = readFSM(tmp);
		//printFSM(fsm);

		//time_t begin = time(NULL);
		clock_t b = clock();


		int res[SIZE];
		//syncTree(fsm, res);
		//int taille = syncTree(fsm, res);
		int taille = bestSearch(fsm, res);

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
	printf("total time: %f ms / average: %f ms (%d fsm)\n", ttl, ttl/(i-1), i-1);

	return 0;

}


int test_ms() {

	char tmp[TAILLE_MAX] = "./data/fsm_hss.fsm";
    printf( "%s\n", tmp );

	FSM * fsm = readFSM(tmp);
	printFSM(fsm);

	printf( "MS(0,1): %d\n", mergSeq(fsm, 0, 1) );	
	printf( "MS(0,2): %d\n", mergSeq(fsm, 0, 2) );	
	printf( "MS(0,3): %d\n", mergSeq(fsm, 0, 3) );
	int tmpS[5] = {1,1,1,1,1};
	//int tmpS[5] = {1,1,1,0,0};	
	printf( "maxMS: %d\n", maxMS(fsm, tmpS) );

	return 0;
}

int test_best() {

	char tmp[TAILLE_MAX] = "./data/fsm_n20_1.fsm";
    	printf( "%s\n", tmp );

	FSM * fsm = readFSM(tmp);
	//printFSM(fsm);

	int res[TAILLE_MAX];

	printf("best: %d\n", bestSearch(fsm, res));

	return 0;
}

int main() {

	test_ST();
	//test_ms();
	//test_best();

	return 0;
}
