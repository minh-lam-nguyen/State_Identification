#ifndef SET_H
#define SET_H

/**
 * Type représentant un Set, cet ensemble contient des int*
 */
typedef struct s_set* Set;

/**
 * Initialisation du Set
 */
Set initSet();

/**
 * Libération d'un Set
 */
void freeSet(Set s);

/**
 * Essaie d'ajouter e à s. Si e est déjà dans s, ne fait rien et renvoie 0. Sinon ajoute e et renvoie 1.
 * e doit être un pointeur sur au moins size entiers créé avec malloc
 * Attention, modifier e après coup peut engendrer des erreurs lors de la recherche ou de l'ajout futur d'éléments.
 */
int add(Set s, int* e, int size);

/**
 * Vérifie si e est dans s. Renvoie 1 si oui et 0 sinon.
 * e doit être un pointeur sur au moins size entiers
 */
int find(Set s, int* e, int size);

#endif
