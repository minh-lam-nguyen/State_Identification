#ifndef HEAP_H
#define HEAP_H

/**
 * Type représentant un Heap, ce tas contient des int*
 * Chaque int* doit contenir une clef. Un noeud du tas
 * contient une clef plus petite que celle de ses fils.
 * De plus le int* doit contenir un espace vide utilisé par le Heap
 **/
typedef struct s_heap* Heap;

/**
 * Initialisation du Heap. 
 * Chaque int* inséré dans ce heap doit avoir 2 cases spéciales d'indices keyIndex et anchorIndex.
 * La première contient la clef permettant de placer l'élément dans le tas et la seconde contient
 * l'indice du int* dans le tas.
 */
Heap initHeap(int keyIndex, int achorIndex);

/**
 * Libération d'un Heap
 */
void freeHeap(Heap s);

/**
 * Insert l'élément e dans le tas h.
 */
int insertToHeap(Heap h, int* e);

/**
 * Décroit la clef de e à k et replace l'élément dans le tas
 */
void decreaseKey(Heap h, int* e, int k);

/**
 * Renvoie la taille de h.
 */
int sizeOfHeap(Heap h);

/**
 * Supprime la racine du tas et la renvoie. Ce noeud contient la clef la plus petite parmi tous les éléments du tas.
 * h doit être non vide.
 */
int* popMinimum(Heap h);

#endif
