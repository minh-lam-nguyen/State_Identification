#ifndef SHP_H
#define SHP_H

#include "fsm.h"

/*
 * stock les merging sequence
 */
typedef struct shp_* Shp;

/*
 * Initialisation des SHP, calcule les shp du fsm.
 */
Shp initSHP(FSM fsm);

/*
 * Libération
 */
void freeSHP(Shp shps, int size);

/**
 * Retourne le maximum des ms pour le sous-ensemble d'états s
 */
int maxShpSorted(Shp shps, int size, int *s);

#endif
