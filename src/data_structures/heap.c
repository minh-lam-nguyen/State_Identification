#include <stdio.h>
#include <stdlib.h>
#include "heap.h"

#define HEAP_SIZE 1048576

struct s_heap{
    int* h[HEAP_SIZE]; // Talbeau contenant HEAP_SIZE éléments
    int kI; // Indice des clefs des éléments (contenus dans les éléments eux-mêmes)
    int aI; // Indice des ancres des éléments (contenus dans les éléments eux-mêmes). L'ancre est l'indice de l'élément dans le tableay h.
    int size; // Taille actuelle du tas
};

Heap initHeap(int keyIndex, int anchorIndex){
    Heap h = (Heap) malloc(sizeof(struct s_heap));
    h->kI = keyIndex;
    h->aI = anchorIndex;
    h->size = 0;
    return h;
}

void freeHeap(Heap h){
    free(h);
}

/**
 * Echange dans h les éléments d'indice eI et fI
 */
void swap(Heap h, int eI, int fI){
    int* e = h->h[eI];
    int* f = h->h[fI];
    h->h[eI] = f;
    h->h[fI] = e;
    e[h->aI] = fI;
    f[h->aI] = eI;
}   

/**
 * Fait remonter l'élément e dans le tas.
 */
void bubbleUp(Heap h, int* e){
    int index = e[h->aI];
    while(1){
        if(index == 0)
            return;
        int fIndex = index / 2;
        int* f = h->h[fIndex];
        if(f[h->kI] > e[h->kI]){
            swap(h, index, fIndex);
            index = fIndex;
        }
        else
            return;
    }
}

/**
 * Fait redescendre l'élément e dans le tas.
 */
void bubbleDown(Heap h, int* e){
    int index = e[h->aI];
    
    while(1){
        int lIndex = index * 2;
        int rIndex = lIndex + 1;
        int fIndex = -1;
        int fKey = -1;
        if(lIndex < h->size){
            int* f = h->h[lIndex];
            if(f[h->kI] < e[h->kI]){
                fIndex = lIndex;
                fKey = f[h->kI];
            }
        }
        if(rIndex < h->size){
            int* f = h->h[rIndex];
            if(f[h->kI] < e[h->kI] && (fKey == -1 || f[h->kI] < fKey)){
                fIndex = rIndex;
            }
        }
        if(fIndex == -1)
            return;
        swap(h, index, fIndex);
        index = fIndex;
    }
}

int insertToHeap(Heap h, int* e){
    if(h->size == HEAP_SIZE)
        return 0;
    h->h[h->size] = e;
    e[h->aI] = h->size;
    h->size++;
    bubbleUp(h, e);
    return 1;
}

void decreaseKey(Heap h, int* e, int k){
    e[h->kI] = k;
    bubbleUp(h, e);
}

int sizeOfHeap(Heap h){
    return h->size;
}

int* popMinimum(Heap h){
    int* min = h->h[0];
    h->size--;
    h->h[0] = h->h[h->size];
    h->h[0][h->aI] = 0;
    bubbleDown(h, h->h[0]);
    return min;
}
