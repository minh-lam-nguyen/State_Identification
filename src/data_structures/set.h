#ifndef SET_H
#define SET_H

/**
 * Type représentant un Set, cet ensemble contient des int*
 */
typedef struct s_set* Set;

/**
 * Initialisation du Set, tous les éléments sont des int* de taille au moins size + offset, seuls les size entiers entre les indices offset et offset + size - 1 de chaque éléments sont utilisés pour vérifier l'unicité.
 */
Set initSet(int size, int offset);

/**
 * Libération d'un Set
 */
void freeSet(Set s);

/**
 * Essaie d'ajouter e à s. Si e est déjà dans s, ne fait rien et renvoie 0. Sinon ajoute e et renvoie 1.
 * e doit être un pointeur sur au moins size entiers, où size est la valeur donnée avec la fonction initSet
 * Attention, modifier e après coup peut engendrer des erreurs lors de la recherche ou de l'ajout futur d'éléments.
 */
int add(Set s, int* e);

/**
 * Vérifie si e est dans s. Renvoie 1 si oui et 0 sinon.
 * e doit être un pointeur sur au moins size entiers où size est la valeur donnée à l'initisalisation du set
 */
int find(Set s, int* e);

int* addOrReturn(Set s, int* e);

#endif
