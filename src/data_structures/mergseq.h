#ifndef MERGSEQ_H
#define MERGSEQ_H

#include "set.h"
#include "fsm.h"

/*
 * stock les merging sequence
 */
typedef struct mergseq_* MergSeq;

/*
 * Initialisation des MS, calcule les ms du fsm.
 * Initialise aussi un tableau de triplets (i, j, ms(ij)) pour tout couple (i, j) d'états du fsm (où ms(ij) est la taille de la MS entre i et j. Ce tableau est ensuite trié selon la 3e coordonnée. 
 */
MergSeq initMS(FSM fsm);

/*
 * Libération
 */
void freeMS(MergSeq mergs, int size);

/*
 * Retourne le maximum des ms pour le sous-ensemble d'états s en regardant l'ensemble des couples d'états de s un par un.
 */
int maxMSLinear(MergSeq mergs, int size, int *s);

/**
 * Retourne le maximum des ms pour le sous-ensemble d'états s et regardant la liste de triples (i, j, ms(ij)) triée crée lors de l'appel à initMS. Parcours cette liste puis s'arrête au premier triplet où i et j sont dans s.  
 */
int maxMSSorted(MergSeq mergs, int size, int *s);

#endif
