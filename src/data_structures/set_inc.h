#ifndef SETINC_H
#define SETINC_H

/**
 * Type représentant un Set, cet ensemble contient des int*
 */
typedef struct s_set_inc* SetInc;

/**
 * Initialisation du Set, tous les éléments sont des int* de taille au moins size + offset, seuls les size entiers entre les indices offset et offset + size - 1 de chaque éléments sont utilisés pour vérifier l'unicité.
 */
SetInc initSetInc(int size, int offset);

/**
 * Libération d'un Set
 */
void freeSetInc(SetInc s);

/**
 * Essaie d'ajouter e à s. Si e est déjà dans s, ne fait rien et renvoie 0. Sinon ajoute e et renvoie 1.
 * e doit être un pointeur sur au moins size entiers, où size est la valeur donnée avec la fonction initSet
 * Attention, modifier e après coup peut engendrer des erreurs lors de la recherche ou de l'ajout futur d'éléments.
 */
int addInc(SetInc s, int* e);

/**
 * Vérifie si e est dans s. Renvoie 1 si oui et 0 sinon.
 * e doit être un pointeur sur au moins size entiers où size est la valeur donnée à l'initisalisation du set
 */
int findInc(SetInc s, int* e);

int* addIncOrReturn(SetInc s, int* e);

#endif
