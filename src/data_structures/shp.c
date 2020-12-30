#include "shp.h"
#include <stdio.h>
#include <stdlib.h>

struct s_couple{
    int a;
    int b;
};

typedef struct s_couple couple;

struct shp_ {
    couple * couples;   
};

int cmpCouples(const void * a, const void * b){
    couple ca = *(couple *)a;
    couple cb = *(couple *)b;
    return ca.b - cb.b;
}


// compute shp for all couple of states and store it in shps
void getShp(FSM fsm, Shp shps){
    int s = get_s(fsm);

    int* queue = (int*) malloc(sizeof(int*) * s);
    for(int i = 0; i < s; i++){
        for(int j = 0; j < s; j++){
            shps->couples[i * s + j].a = j;
            shps->couples[i * s + j].b = s + 1;
        }
        shps->couples[i * s + i].b = 0;

        int qIndex = 0;
        int qSize = 0;
        
        queue[0] = i;
        qSize = 1;

        while(qIndex != qSize){
            int u = queue[qIndex];
            qIndex++;
            int d = shps->couples[i * s + u].b;

            for(int input = 0; input < get_i(fsm); input++){
                int v = get_succ(fsm, u, input);
                if(shps->couples[i * s + v].b != s + 1)
                    continue;
                shps->couples[i * s + v].b = d + 1;
                queue[qSize] = v;
                qSize++;
            }
        }
    }


    free(queue);

    for(int i = 0; i < s; i++)
        qsort(shps->couples + (i * s), s, sizeof(couple), cmpCouples);
    
}


Shp initSHP(FSM fsm){
    Shp shps = (Shp)malloc( sizeof(struct shp_) );
    int size = get_s(fsm);

    shps->couples = (couple *) malloc(sizeof(couple) * size * size);

    getShp(fsm, shps);

    if (shps->couples == NULL)
        return NULL;
    else
        return shps;
}


void freeSHP(Shp shps, int size){
    free(shps->couples);
    free(shps);
}



int maxShpSorted(Shp shps, int size, int* s){
    int max = -1;
    for(int i = 0; i < size; i++){
        int access = 0; // 0 si s n'est pas accessible et 1 sinon
        for(int j = 0; j < size; j++){
            couple c = shps->couples[i * size + j];
            if(c.b == size + 1)
                break;
            if(s[c.a] != 1)
                continue;
            access = 1;
            if(c.b > max)
                max = c.b;
            break;
        }
        if(access == 0)
            return -1;
    }
    return max;
}
