#include <stdlib.h>
#include <stdio.h>
#include "set.h"

int main(void){
   Set s = initSet();

   int a1[4] = {1,1,1,1};
   int a2[4] = {1,0,1,1};
   int a3[4] = {1,1,0,1};
   int a4[4] = {1,0,1,1};
   int a5[4] = {1,1,0,1};

   int d1 = add(s, a1, 4);
   int d2 = add(s, a2, 4);
   int d3 = add(s, a3, 4);
   int d4 = add(s, a4, 4);
   int d5 = add(s, a5, 4);

   printf("%d %d %d %d %d\n", d1, d2, d3, d4, d5);

   int e1 = find(s, a1, 4);
   int e2 = find(s, a2, 4);
   int e3 = find(s, a3, 4);
   int e4 = find(s, a4, 4);
   int e5 = find(s, a5, 4);
   printf("%d %d %d %d %d\n", e1, e2, e3, e4, e5);

   freeSet(s);
}
