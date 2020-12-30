#include "tergseq.h"
#include <stdio.h>
#include <stdlib.h>

struct s_quadruple{
    int a;
    int b;
    int c;
    int d;
};
typedef struct s_quadruple quadruple;

struct tergseq_ {
    int* ts;
    quadruple * quadruples;   
};

int cmpQuadruples(const void * a, const void * b){
    quadruple ta = *(quadruple *)a;
    quadruple tb = *(quadruple *)b;
    return tb.d - ta.d;
}


// compute merging seq for all couple of states and store it in ms
void getTergSeq(FSM fsm, TergSeq tergs){
    int s = get_s(fsm);

    int* queue = (int*) malloc(sizeof(int*) * s * s * s * 2);
    int qIndex = 0;
    int qSize = 0;

    // Init terging sequences
    for(int i = 0; i < get_s(fsm); i++){
        for(int j = i; j < get_s(fsm); j++)
            for(int k = j; k < get_s(fsm); k++)
                tergs->ts[i * s * s + j * s + k] = -1;
        tergs->ts[i * s * s + i * s + i] = 0;
        queue[3 * qSize] = i;
        queue[3 * qSize + 1] = i;
        queue[3 * qSize + 2] = i;
        qSize++;
    }


    while(qIndex != qSize){
        int u = queue[3 * qIndex];
        int v = queue[3 * qIndex + 1];
        int w = queue[3 * qIndex + 2];
        qIndex++;

        for(int input = 0; input < get_i(fsm); input++){
            int nbUPred = get_nb_preds(fsm, u, input);
            int nbVPred = get_nb_preds(fsm, v, input);
            int nbWPred = get_nb_preds(fsm, w, input);
            for(int pui = 0; pui < nbUPred; pui++){
                int pu = get_pred(fsm, u, input, pui);
                int fpvi = 0;
                if(u == v)
                    fpvi = pui;
                for(int pvi = fpvi; pvi < nbVPred; pvi++){
                    int pv = get_pred(fsm, v, input, pvi);
                    for(int pwi = 0; pwi < nbWPred; pwi++){   
                        int pw = get_pred(fsm, w, input, pwi);
                        int ppu = pu;
                        int ppv = pv;
                        int ppw = pw;
                        if(ppu > ppv){ ppu = ppu + ppv; ppv = ppu - ppv; ppu = ppu - ppv;}
                        if(ppu > ppw){ ppu = ppu + ppw; ppw = ppu - ppw; ppu = ppu - ppw;}
                        if(ppv > ppw){ ppv = ppv + ppw; ppw = ppv - ppw; ppv = ppv - ppw;}
                        if(tergs->ts[ppu * s * s + ppv * s + ppw] != -1)
                            continue;
                        tergs->ts[ppu * s * s + ppv * s + ppw] = tergs->ts[u * s * s + v * s + w] + 1;
                        queue[3 * qSize] = ppu;
                        queue[3 * qSize + 1] = ppv;
                        queue[3 * qSize + 2] = ppw;
                        qSize++;
                    }
                }
            }
        }
    }

    free(queue);

    int count = 0;
    for (int i=0; i<get_s(fsm); i++){
        for (int j=i; j<get_s(fsm); j++){
            for (int k=j; k<get_s(fsm); k++){
                if (tergs->ts[i * s * s + j * s + k] == -1){
                    tergs->ts = NULL;
                    return;
                }
                tergs->quadruples[count].a = i;
                tergs->quadruples[count].b = j;
                tergs->quadruples[count].c = k;
                tergs->quadruples[count].d = tergs->ts[i * s * s + j * s + k];
                count++;
            }
        }
    }

    qsort(tergs->quadruples, count, sizeof(quadruple), cmpQuadruples);
}


TergSeq initTS(FSM fsm){
    TergSeq tergs = (TergSeq)malloc( sizeof(struct tergseq_) );
    int size = get_s(fsm);

    tergs->ts = (int*)malloc( size * size * size * sizeof(int*) );
    tergs->quadruples = (quadruple *) malloc(sizeof(quadruple) * size * size * size);

    getTergSeq(fsm, tergs);

    if (tergs->ts == NULL)
        return NULL;
    else
        return tergs;
}


void freeTS(TergSeq tergs, int size){
    free(tergs->ts);
    free(tergs->quadruples);
    free(tergs);
}



int maxTSSorted(TergSeq tergs, int size, int* s){
    int n = size * size * size;
    for(int i = 0; i < n; i++){
        quadruple q = tergs->quadruples[i];
        if(s[q.a] != 1 || s[q.b] != 1 || s[q.c] != 1)
            continue;
        return q.d;
    }
    return 0;
}
