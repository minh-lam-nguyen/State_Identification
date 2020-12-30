#include "quadseq.h"
#include <stdio.h>
#include <stdlib.h>

struct s_quintuple{
    int a;
    int b;
    int c;
    int d;
    int e;
};
typedef struct s_quintuple quintuple;

struct quadseq_ {
    int* qs;
    quintuple * quintuples;   
};

int cmpQuintuples(const void * a, const void * b){
    quintuple ta = *(quintuple *)a;
    quintuple tb = *(quintuple *)b;
    return tb.e - ta.e;
}


// compute merging seq for all couple of states and store it in ms
void getQuadSeq(FSM fsm, QuadSeq quads){
    int s = get_s(fsm);

    int* queue = (int*) malloc(sizeof(int*) * s * s * s * s * 4);
    int qIndex = 0;
    int qSize = 0;

    // Init quading sequences
    for(int i = 0; i < s; i++){
        for(int j = i; j < s; j++)
            for(int k = j; k < s; k++)
                for(int l = k; l < s; l++)
                    quads->qs[i * s * s * s + j * s * s + k * s + l] = -1;
        quads->qs[i * s * s * s + i * s * s + i * s + i] = 0;
        queue[4 * qSize] = i;
        queue[4 * qSize + 1] = i;
        queue[4 * qSize + 2] = i;
        queue[4 * qSize + 3] = i;
        qSize++;
    }


    while(qIndex != qSize){
        int u = queue[4 * qIndex];
        int v = queue[4 * qIndex + 1];
        int w = queue[4 * qIndex + 2];
        int x = queue[4 * qIndex + 3];
        qIndex++;

        for(int input = 0; input < get_i(fsm); input++){
            int nbUPred = get_nb_preds(fsm, u, input);
            int nbVPred = get_nb_preds(fsm, v, input);
            int nbWPred = get_nb_preds(fsm, w, input);
            int nbXPred = get_nb_preds(fsm, x, input);
            int dist = quads->qs[u * s * s * s + v * s * s + w * s + x];
            for(int pui = 0; pui < nbUPred; pui++){
                int pu = get_pred(fsm, u, input, pui);
                int fpvi = 0;
                if(u == v)
                    fpvi = pui;
                for(int pvi = fpvi; pvi < nbVPred; pvi++){
                    int pv = get_pred(fsm, v, input, pvi);
                    for(int pwi = 0; pwi < nbWPred; pwi++){   
                        int pw = get_pred(fsm, w, input, pwi);
                        for(int pxi = 0; pxi < nbXPred; pxi++){   
                            int px = get_pred(fsm, x, input, pxi);
                            int ppu = pu;
                            int ppv = pv;
                            int ppw = pw;
                            int ppx = px;
                            if(ppu > ppv){ ppu = ppu + ppv; ppv = ppu - ppv; ppu = ppu - ppv;}
                            if(ppu > ppw){ ppu = ppu + ppw; ppw = ppu - ppw; ppu = ppu - ppw;}
                            if(ppu > ppx){ ppu = ppu + ppx; ppx = ppu - ppx; ppu = ppu - ppx;}
                            if(ppv > ppw){ ppv = ppv + ppw; ppw = ppv - ppw; ppv = ppv - ppw;}
                            if(ppv > ppx){ ppv = ppv + ppx; ppx = ppv - ppx; ppv = ppv - ppx;}
                            if(ppw > ppx){ ppw = ppw + ppx; ppx = ppw - ppx; ppw = ppw - ppx;}
                            int index = ppu * s * s * s + ppv * s * s + ppw * s + ppx;
                            if(quads->qs[index] != -1)
                                continue;
                            quads->qs[index] = dist + 1;
                            queue[4 * qSize] = ppu;
                            queue[4 * qSize + 1] = ppv;
                            queue[4 * qSize + 2] = ppw;
                            queue[4 * qSize + 3] = ppx;
                            qSize++;
                        }
                    }
                }
            }
        }
    }

    free(queue);

    int count = 0;
    for (int i=0; i<s; i++){
        for (int j=i; j<s; j++){
            for (int k=j; k<s; k++){
                for (int l=k; l<s; l++){
                    if (quads->qs[i * s * s * s + j * s * s + k * s + l] == -1){
                        quads->qs = NULL;
                        return;
                    }
                    quads->quintuples[count].a = i;
                    quads->quintuples[count].b = j;
                    quads->quintuples[count].c = k;
                    quads->quintuples[count].d = l;
                    quads->quintuples[count].e = quads->qs[i * s * s * s + j * s * s + k * s + l];
                    count++;
                }
            }
        }
    }

    qsort(quads->quintuples, count, sizeof(quintuple), cmpQuintuples);
}


QuadSeq initQS(FSM fsm){
    QuadSeq quads = (QuadSeq)malloc( sizeof(struct quadseq_) );
    int size = get_s(fsm);

    quads->qs = (int*)malloc( size * size * size * size * sizeof(int*) );
    quads->quintuples = (quintuple *) malloc(sizeof(quintuple) * size * size * size * size);
    getQuadSeq(fsm, quads);

    if (quads->qs == NULL)
        return NULL;
    else
        return quads;
}


void freeQS(QuadSeq quads, int size){
    free(quads->qs);
    free(quads->quintuples);
    free(quads);
}



int maxQSSorted(QuadSeq quads, int size, int* s){
    int n = size * size * size;
    for(int i = 0; i < n; i++){
        quintuple q = quads->quintuples[i];
        if(s[q.a] != 1 || s[q.b] != 1 || s[q.c] != 1 || s[q.d] != 1)
            continue;
        return q.e;
    }
    return 0;
}
