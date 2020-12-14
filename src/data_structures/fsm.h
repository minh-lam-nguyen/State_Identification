#ifndef FSM_H
#define FSM_H

/**
 * Type FSM
 */
typedef struct fsm_* FSM;

/**
 * Initialisation du FSM à partir d'un fichier
 */
FSM initFSM(const char * filename);


/**
 * Libération du FSM
 */
void freeFSM(FSM f);


/**
 * Affiche le FSM 
 */
void printFSM(FSM f);

/**
 * get nb states s
 */
int get_s(FSM f);

/**
 * get nb input i
 */
int get_i(FSM f);

/**
 * get successor of state s after applying input i
 */
int get_succ(FSM f, int s, int i);


#endif
