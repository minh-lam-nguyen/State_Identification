#include <stdlib.h>
#include <stdio.h>
#include "data_structures/heap.h"

int main(void){
   Heap s = initHeap(0, 2);

   int a1[4] = {4,1,1,6};
   int a2[4] = {3,2,1,7};
   int a3[4] = {2,3,0,8};
   int a4[4] = {7,4,1,9};
   int a5[4] = {2,5,0,10};

   int d1 = insertToHeap(s, a1); 
   int d2 = insertToHeap(s, a2);
   int d3 = insertToHeap(s, a3);
   int d4 = insertToHeap(s, a4);
   int d5 = insertToHeap(s, a5);

   printf("%d %d %d %d %d\n", d1, d2, d3, d4, d5);

    decreaseKey(s, a1, 1);

    for(int i = 0; i < 5; i++){
        int* a = popMinimum(s);
        printf("%d %d %d %d\n", a[0], a[1], a[2], a[3]);
    }

   freeHeap(s);
}
