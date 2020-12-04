#include "set.h"
#include <stdlib.h>
#include <stdio.h>
#include <set>
#include <algorithm>


int GLOBAL_SIZE = 0;

struct Compare
{
    bool operator ()(int * e, int * f)
    {
        for(int i = 0; i < GLOBAL_SIZE; i++){
            if(e[i] < f[i])
                return true;
            if(e[i] > f[i])
                return false;
        }
        return false;
    }
};

struct s_set{
    std::set<int *, Compare> v;
};

Set initSet(int size){
    Set s = new s_set();
    GLOBAL_SIZE = size;
    return s;
}

void freeSet(Set s){
    delete(s);
}

int add(Set s, int* e){
   if(find(s, e)){
       return 0;
   }

   s->v.insert(e);
   return 1;
}

int find(Set s, int* e){
    int found =  s->v.find(e) != s->v.end();
    return found;
}
