#include "tupseq.h"
#include <stdio.h>
#include <stdlib.h>

struct s_couple{
    int index;
    int key;
};

typedef struct s_couple couple;

struct tupseq_ {
    int* ts;
    int* tuples;
    couple* couples;
    int k;
};

int cmpTuples(const void * a, const void * b){
    couple ta = *(couple *)a;
    couple tb = *(couple *)b;
    return tb.key - ta.key;
}


// compute merging seq for all couple of states and store it in ms
void getTupSeq(FSM fsm, TupSeq mergs){
    
    int s = get_s(fsm);
    int count = 1;
    for(int k = 0; k < tups->k; k++)
        count *= size;
    int* queue = (int*) malloc(sizeof(int) * count * k);
    int qIndex = 0;
    int qSize = 0;

    for(int i = 0; i < count; i++){
        tups->ts[i] = -1;
    }

    // Init merging sequences
    for(int i = 0; i < get_s(fsm) - 1; i++){
        int index = 0;
        int counter = 1;
        for(int k = 0; k < tups->k; k++){
            index += i * counter;
            counter *= size;
            queue[i * tups->k + k] =  i;
        }
        tups->ts[index] = 0;
        qSize++;
    }


    while(qIndex != qSize){
        int* tuple = queue + tups->k * qIndex;
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
    */
}


TupSeq initTuS(FSM fsm, int k){
    TupSeq tups = (TupSeq)malloc( sizeof(struct tupseq_) );
    tups->k = k;
    int size = get_s(fsm);

    int count = 1;
    for(int i = 0; i < k; i++)
        count *= size;

    tups->ts = (int*) malloc( count * sizeof(int*) );
    tups->tuples = (int *) malloc(k * count);
    tups->couples = (couple *) malloc(sizeof(couple) * count);
    getTupSeq(fsm, tups);

    if (tups->ts == NULL)
        return NULL;
    else
        return tups;
}


void freeMS(TupSeq tups, int size){
    free(tups->ts);
    free(tups->couples);
    free(tups->tuples);
    free(tups);
}


// compute max merging seq (lower bound) for states s
int maxTuSLinear(TupSeq tups, int size, int *s){

    int res = -1;
    int* t = (int*) malloc(sizeof(int) * tups->k);
    for(int k = 0; k < tups->k; k++)
        t[k] = -1;

    int current = 0;
    while(1){
        
        while(t[current] == -1 || (t[current] != size && s[t[current]] == 0))
            t[current]++;

        if(t[current] == size){
            if(current == 0)
                break;
            t[current] = -1;
            current--;
            continue;
        }

        current++;
        if(current == tups->k){
            int counter = 1;
            int index = 0;
            for(int k = 0; k < tups->k; k++){
                index += t[tups->k - 1 - k] * counter;
                counter *= size;
                if(res < tups->ts[index])
                    res = tups->ts[index];
            }
            current--;
        }

    }
    free(t);

    return res;
}

int maxMSSorted(TupSeq tups, int size, int* s){
    int count = 1;
    for(int k = 0; k < tups->k ; k++)
        count *= size;
    for(int i = 0; i < count; i++){
        couple c = tups->couples[i];
        int* tuple = tups->tuples + c.index * tups->k;
        for(int k = 0; k < tups->k; k++){
            if(s[tuple[k]] == 0)
                goto end;
        }
        return c.key;
end: 
        continue;
    }
    return 0;
}
