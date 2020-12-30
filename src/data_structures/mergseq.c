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
    int* ms;
    triple * triples;   
};

int cmpTriples(const void * a, const void * b){
    triple ta = *(triple *)a;
    triple tb = *(triple *)b;
    return tb.c - ta.c;
}


// compute merging seq for all couple of states and store it in ms
void getMergSeq(FSM fsm, MergSeq mergs){
    int s = get_s(fsm);

    int* queue = (int*) malloc(sizeof(int*) * s * s * 2);
    int qIndex = 0;
    int qSize = 0;

    // Init merging sequences
    for(int i = 0; i < get_s(fsm) - 1; i++){
        mergs->ms[i * s + i] = 0;
        queue[2 * qSize] = i;
        queue[2 * qSize + 1] = i;
        qSize++;
        for(int j = i + 1; j < get_s(fsm); j++)
            mergs->ms[i * s + j] = -1;
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
                        ppu = ppu + ppv;
                        ppv = ppu - ppv;
                        ppu = ppu - ppv;
                    }
                    if(mergs->ms[ppu * s + ppv] != -1)
                        continue;
                    mergs->ms[ppu * s + ppv] = mergs->ms[u * s + v] + 1;
                    queue[2 * qSize] = ppu;
                    queue[2 * qSize + 1] = ppv;
                    qSize++;
                }
            }
        }
    }

    free(queue);



    int count = 0;
    for (int i=0; i<get_s(fsm)-1; i++){
        for (int j=i+1; j<get_s(fsm); j++){
            if (mergs->ms[i * s + j] == -1){
                mergs->ms = NULL;
                return;
            }
            mergs->triples[count].a = i;
            mergs->triples[count].b = j;
            mergs->triples[count].c = mergs->ms[i * s + j];
            count++;
        }
    }

    qsort(mergs->triples, count, sizeof(triple), cmpTriples);
}


MergSeq initMS(FSM fsm){
    MergSeq mergs = (MergSeq)malloc( sizeof(struct mergseq_) );
    int size = get_s(fsm);

    mergs->ms = (int*)malloc( size * size * sizeof(int*) );
    mergs->triples = (triple *) malloc(sizeof(triple) * size * (size - 1) / 2);

    getMergSeq(fsm, mergs);

    if (mergs->ms == NULL)
        return NULL;
    else
        return mergs;
}


void freeMS(MergSeq mergs, int size){
    free(mergs->ms);
    free(mergs->triples);
    free(mergs);
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
            int tmp = mergs->ms[i * size + j] ;
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
