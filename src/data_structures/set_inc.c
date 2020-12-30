#include "set_inc.h"
#include <stdlib.h>
#include <stdio.h>

#define HASH_SIZE 4096 
#define SIZE_V 1024

struct s_set_inc{
    int* s[HASH_SIZE][SIZE_V]; // Tableau contenant HASH_SIZE liste de SIZEV éléments
    int c[HASH_SIZE]; // Nombre déléments avec chaque hash dans le SetInc
    int size; // Nombre minimum d'entiers dans les int* ajoutés dans le SetInc
    int offSetInc;
};

/**
 * Renvoie le Hash de l'élément e, sachant qu'il est de taille s
 */
unsigned int hash(SetInc s, int* e) {
    unsigned int res = 0;
    for(int i = s->offSetInc; i < s->offSetInc + s->size; i++) {
        res = res * 1997 + i * e[i];
    }
    return res % HASH_SIZE;
}

int equals(SetInc s, int* e, int* f){
    for(int i = s->offSetInc; i < s->offSetInc + s->size; i++){
        if(e[i] > f[i])
            return 0;
    }
    return 1;
}

/**
 * Initialisation du SetInc
 */
SetInc initSetInc(int size, int offSetInc){
    SetInc s = (SetInc) malloc(sizeof(struct s_set_inc));
    for (int j=0; j<HASH_SIZE; j++) {
        s->c[j] = 0;
    }
    s->size = size;
    s->offSetInc = offSetInc;
    return s;
}

/**
 * Libération d'un SetInc
 */
void freeSetInc(SetInc s){
    // Libérer les pointeurs
    for(int i=0; i<HASH_SIZE; i++)
        for(int j=0; j<s->c[i]; j++)
            free(s->s[i][j]);
    free(s); // On ne libère pas les éléments de s->s car ils ont été créé ailleurs et seront libérés ailleurs.
}

/**
 * Essaie d'ajouter e à s. Si e est déjà dans s, ne fait rien et renvoie 0. Sinon ajoute e et renvoie 1.
 * e doit être un pointeur sur au moins size entiers créé avec malloc
 * Attention, modifier e après coup peut engendrer des erreurs lors de la recherche ou de l'ajout futur d'éléments.
 */
int addInc(SetInc s, int* e){
    int h = hash(s, e);
    int c = s->c[h];
    for(int i = 0; i < c; i++)
        if(equals(s, e, s->s[h][i]))
            return 0;

    if(c == SIZE_V)
        return -1; // Capacité max atteinte
    
    s->s[h][c] = e;
    s->c[h]++;
    return 1;
}

/**
 * Vérifie si e est dans s. Renvoie 1 si oui et 0 sinon.
 * e doit être un pointeur sur au moins size entiers
 */
int findInc(SetInc s, int* e){
    int h = hash(s, e);
    int c = s->c[h];
    for(int i = 0; i < c; i++){
        if(equals(s, e, s->s[h][i]))
            return 1;
    }
    return 0;
}


int* addIncOrReturn(SetInc s, int* e){
    int h = hash(s, e);
    int c = s->c[h];
    for(int i = 0; i < c; i++)
        if(equals(s, e, s->s[h][i]))
            return s->s[h][i];

    if(c == SIZE_V)
        return NULL; // Capacité max atteinte
    
    s->s[h][c] = e;
    s->c[h]++;
    return e;
}
