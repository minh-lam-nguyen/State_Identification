#include "set.h"
#include <stdlib.h>
#include <stdio.h>
#include <set>
#include <algorithm>


struct elem{
    elem(int* p, int s) : e(p), size(s) {}
    int* e;
    int size;
};

struct Compare
{
    bool operator ()(struct elem e, struct elem f)
    {
        for(int i = 0; i < e.size; i++){
            if(e.e[i] < f.e[i])
                return true;
            if(e.e[i] > f.e[i])
                return false;
        }
        return false;
    }
};

struct s_set{
    std::set<struct elem, Compare> v;
    int size;
};

Set initSet(int size){
    Set s = new s_set();
    s->size = size;
    return s;
}

void freeSet(Set s){
    delete(s);
}

int add(Set s, int* e){
   if(find(s, e)){
       return 0;
   }

   struct elem * f = new elem(e, s->size);
   s->v.insert(*f);
   return 1;
}

int find(Set s, int* e){
    struct elem * f = new elem(e, s->size);
    int found =  s->v.find(*f) != s->v.end();
    delete(f);
    return found;
}
