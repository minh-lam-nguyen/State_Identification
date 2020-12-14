#include "set.h"
#include <stdlib.h>
#include <stdio.h>

#define HASH_SIZE 4
#define SIZE_V 10

struct s_set{
    int current_hash_index;
    int max_hash;
    int** s; // Tableau contenant HASH_SIZE listes d'int*
    int* c; // Nombre déléments avec chaque hash dans le set
    int size; // Nombre minimum d'entiers dans les int* ajoutés dans le set
};

/**
 * Renvoie le Hash de l'élément e, sachant qu'il est de taille s
 */
unsigned int hash(Set s, int* e, int size) {
    unsigned int res = 0;
    for(int i = 0; i < size; i++) {
        res = (res + i) * 1997 + e[i];
    }
    return res % (s->max_hash);
}

int equals(int* e, int* f, int size){
    for(int i = 0; i < size; i++){
        if(e[i] != f[i])
            return 0;
    }
    return 1;
}

/**
 * Initialisation du Set
 */
Set initSet(int size){
    Set s = (Set) malloc(sizeof(struct s_set));
    s->current_hash_index = -1;
    s->s = (int**) malloc(sizeof(int*) * HASH_SIZE * SIZE_V);
    s->c = (int*) malloc(sizeof(int) * HASH_SIZE);
    
    for (int j=0; j<HASH_SIZE; j++) {
        s->c[j] = 0;
    }
    
    s->size = size;
    s->max_hash = HASH_SIZE;
    return s;
}

/**
 * Libération d'un Set
 */
void freeSet(Set s){
    free(s); // On ne libère pas les éléments de s->s car ils ont été créé ailleurs et seront libérés ailleurs.
}

int double_size(Set s){
    s->current_hash_index += 1;
    int new_hash = 0;
    if(s->current_hash_index < 10){
        int hashes[] = {8191, 16381, 33769, 67547, 135101, 270209, 540433, 1080899, 2161813, 4323629}; // Si on a besoin de plus on est mal barré quand même.
        new_hash = hashes[s->current_hash_index];
    }
    else
        new_hash = 2 * s->max_hash;

    int** double_s = (int**) malloc(new_hash * SIZE_V * sizeof(int*));
    int* double_c = (int*) malloc(new_hash * sizeof(int));
    

    if(double_s == NULL || double_c == NULL)
        return 0;

    int prev_hash = s->max_hash;
    s->max_hash = new_hash;
    for(int i = 0; i < prev_hash; i++){
        for(int j = 0; j < s->c[i]; j++){
            int* e = s->s[i * SIZE_V + j];
            int h = hash(s, e, s->size);
            int c = double_c[h];
            double_s[h * SIZE_V  + c] = e;
            double_c[h]++;
        }
    }
    
    free(s->s);
    free(s->c);
    s->s = double_s;
    s->c = double_c;
    return 1;
}

/**
 * Essaie d'ajouter e à s. Si e est déjà dans s, ne fait rien et renvoie 0. Sinon ajoute e et renvoie 1.
 * e doit être un pointeur sur au moins size entiers créé avec malloc
 * Attention, modifier e après coup peut engendrer des erreurs lors de la recherche ou de l'ajout futur d'éléments.
 */
int add(Set s, int* e){
    while(1){
        int h = hash(s, e, s->size);
        int c = s->c[h];
        for(int i = 0; i < c; i++)
            if(equals(e, s->s[h * SIZE_V + i], s->size))
                return 0;

        if(c == SIZE_V){
            int success = double_size(s);
            if(success == 0)
                return -1; // Capacité max atteinte
        }
        else{
            s->s[h * SIZE_V + c] = e;
            s->c[h]++;
            break;
        }
    }
    
    return 1;
}

/**
 * Vérifie si e est dans s. Renvoie 1 si oui et 0 sinon.
 * e doit être un pointeur sur au moins size entiers
 */
int find(Set s, int* e){
    int h = hash(s, e, s->size);
    int c = s->c[h];
    for(int i = 0; i < c; i++){
        if(equals(e, s->s[h * SIZE_V + i], s->size))
            return 1;
    }
    return 0;
}
