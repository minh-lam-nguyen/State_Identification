#ifndef MERGSEQ_H
#define MERGSEQ_H

#include "set.h"
#include "fsm.h"

/*
 * stock les merging sequence
 */
typedef struct mergseq_* MergSeq;

/*
 * Initialisation des MS, calcule les ms du fsm 
 */
MergSeq initMS(FSM fsm);

/*
 * Libération
 */
void freeMS(MergSeq mergs, int size);

/*
 * Retourne le maximum des ms pour le sous-ensemble d'état s
 */
int maxMS(MergSeq mergs, int size, int *s);

#endif
