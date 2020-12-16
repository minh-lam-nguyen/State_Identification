#include "set.h"
#include <stdlib.h>
#include <stdio.h>
#include <unordered_set>
#include <algorithm>

#define HASH_SIZE 4096 

struct elem{
    elem(int* p, int s) : e(p), size(s) {}
    int* e;
    int size;


};

struct Equal
{
    bool operator ()(struct elem e, struct elem f) const
    {
        for(int i = 0; i < e.size; i++){
            if(e.e[i] != f.e[i])
                return false;
        }
        return true;
    }
};

struct Hash
{
    size_t operator ()(struct elem e) const
    {
      size_t res = 0;
      for(int i = 0; i < e.size; i++) {
          res = res * 1997 + i * e.e[i];
      }
      return res % HASH_SIZE;
      //return res % std::numeric_limits<size_t>::max();
    }
};
/*
size_t hash (struct elem &e) {
  size_t res = 0;
  for(int i=0; i < e.size; i++){
    res = res * 1997 + i * e.e[i];
  }
  return res % HASH_SIZE;
}
*/

struct s_set{
    std::unordered_set<struct elem, Hash, Equal> v;
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
   //s->v.emplace(e, s->size);
   return 1;
}

int find(Set s, int* e){
    struct elem * f = new elem(e, s->size);
    int found =  s->v.find(*f) != s->v.end();
    delete(f);
    return found;
}
