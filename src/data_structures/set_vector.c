#include "set.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <algorithm>

struct s_set{
    std::vector<int*> v;
    int size;
};

Set initSet(int size){
    Set s = new s_set();
    s->size = size;
    return s;
}

int equals(int* e, int* f, int size){
    for(int i = 0; i < size; i++){
        if(e[i] != f[i])
            return 0;
    }
    return 1;
}

void freeSet(Set s){
    delete(s);
}

int add(Set s, int* e){
   if(find(s, e))
       return 0;
   s->v.push_back(e);
   return 1;
}

int find(Set s, int* e){
    return std::find_if(s->v.begin(), s->v.end(), [&s, &e](int* f){return equals(e, f, s->size);}) != s->v.end();
}
