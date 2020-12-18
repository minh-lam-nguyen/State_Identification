#include "mergseq.h"
#include <stdio.h>
#include <stdlib.h>

struct s_triple{
    int a;
    int b;
    int c;
};
typedef struct s_triple triple;

struct mergseq_ {
    int **ms;
    triple * triples;   
};

// compute min merging seq for (s1,s2), with s1 <= s2
// return length of merging seq for (s1, s2), -1 if no merging (or error)
int mergSeq(FSM fsm, int s1, int s2){

	// elem of queue: ( state1, state2, length ) 
	int **queue = (int**)malloc( get_s(fsm)*get_s(fsm) * sizeof(int*) );

	if (!queue)
		return -1;
	int front_q = 0;
	int back_q = 1;

	queue[0] = (int*)malloc( 3 * sizeof(int) );
	queue[0][0] = s1;
	queue[0][1] = s2;
	queue[0][2] = 0;

	// elem of visited: ( state1, state2 )
	Set visited = initSet(2, 0);

	while (front_q != back_q) {

		int current = front_q;
		front_q++;
		//printf("\ncurrent/f/b: %d %d %d\n", current, front_q, back_q);

		int *tmp = queue[current];
		//printf("tmp: %d %d %d\n", tmp[0], tmp[1], tmp[2]);

		// s1 and s2 merged 
		if (tmp[0] == tmp[1]){
			int res = tmp[2];


			for (int iq=0; iq<back_q; iq++)
				free(queue[iq]);
			free(queue);
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
        for (int i=0; i<get_i(fsm); i++){
        	//printf("add to q: %d %d %d\n", get_succ(fsm, tmp[0], i), get_succ(fsm, tmp[1], i), tmp[2]+1);
	        queue[back_q] = (int*)malloc( 3 * sizeof(int) );
        	queue[back_q][0] = get_succ(fsm, tmp[0], i);
        	queue[back_q][1] = get_succ(fsm, tmp[1], i);
        	queue[back_q][2] = tmp[2]+1;
        	back_q++;
        }

	}

	for (int iq=0; iq<back_q; iq++)
		free(queue[iq]);
	free(queue);
	freeSet(visited);

	return -1;
}

int cmpTriples(const void * a, const void * b){
    triple ta = *(triple *)a;
    triple tb = *(triple *)b;
    return tb.c - ta.c;
}


// compute merging seq for all couple of states and store it in ms
void getMergSeq(FSM fsm, MergSeq mergs){
    int** ms = mergs->ms;
    int s = get_s(fsm);

    int* queue = (int*) malloc(sizeof(int*) * s * s * 2);
    int qIndex = 0;
    int qSize = 0;

    // Init merging sequences
    for(int i = 0; i < get_s(fsm) - 1; i++){
        ms[i][i] = 0;
        queue[2 * qSize] = i;
        queue[2 * qSize + 1] = i;
        qSize++;
        for(int j = i + 1; j < get_s(fsm); j++)
            ms[i][j] = -1;
    }

    while(qIndex != qSize){
        int u = queue[2 * qIndex];
        int v = queue[2 * qIndex + 1];
        qIndex++;

        for(int input = 0; input < get_i(fsm); input++){
            int nbUPred = get_nb_preds(fsm, u, input);
            int nbVPred = get_nb_preds(fsm, v, input);
            for(int pui = 0; pui < nbUPred; pui++){
                int pu = get_pred(fsm, u, input, pui);
                for(int pvi = 0; pvi < nbVPred; pvi++){
                    int pv = get_pred(fsm, v, input, pvi);
                    int ppu = pu;
                    int ppv = pv;
                    if(ppu > ppv){ // Swap
                        ppu = ppu + ppv
                        ppv = ppu - ppv;
                        ppu = ppu - ppv;
                    }
                    if(ms[ppu][ppv] != -1)
                        continue;
                    ms[ppu][ppv] = ms[u][v] + 1;
                    queue[2 * qSize] = ppu;
                    queue[2 * qSize + 1] = ppv;
                }
            }
        }
    }

    free(queue);

    int count = 0;
    for (int i=0; i<get_s(fsm)-1; i++){
        for (int j=i+1; j<get_s(fsm); j++){
            int ms = mergSeq(fsm, i, j);
            if (ms == -1){
                mergs->ms = NULL;
                return;
            }
            mergs->ms[i][j] = ms;
            mergs->triples[count].a = i;
            mergs->triples[count].b = j;
            mergs->triples[count].c = ms[i][j];
            count++;
        }
    }

    qsort(mergs->triples, count, sizeof(triple), cmpTriples);
}


MergSeq initMS(FSM fsm){
    MergSeq mergs = (MergSeq)malloc( sizeof(struct mergseq_) );
    int size = get_s(fsm);

    mergs->ms = (int**)malloc( size * sizeof(int*) );
    mergs->triples = (triple *) malloc(sizeof(triple) * size * (size - 1) / 2);

    for (int i=0; i<size; i++)
            mergs->ms[i] = (int*)malloc( size * sizeof(int) );

    getMergSeq(fsm, mergs);

    if (mergs->ms == NULL)
        return NULL;
    else
        return mergs;
}


void freeMS(MergSeq mergs, int size){
	for(int i=0; i<size; i++)
		free(mergs->ms[i]);
	free(mergs->ms);
    free(mergs->triples);
}



// compute max merging seq (lower bound) for states s
int maxMSLinear(MergSeq mergs, int size, int *s){

    int res = -1;
    for (int i=0; i<size; i++) {
        if(!s[i])
            continue;

        for (int j=i+1; j<size; j++) {
            if (!s[j])
                continue;
            int tmp = mergs->ms[i][j] ;
            if (tmp > res)
                res = tmp;	
        }
    }

    return res;
}

int maxMSSorted(MergSeq mergs, int size, int* s){
    int n = size * (size - 1) / 2;
    for(int i = 0; i < n; i++){
        triple t = mergs->triples[i];
        if(s[t.a] != 1 || s[t.b] != 1)
            continue;
        return t.c;
    }
    return 0;
}
