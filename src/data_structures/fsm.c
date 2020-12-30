#include "fsm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_NB_STATES 500
#define NB_SUCC 2 // nb of successors
#define TAILLE_MAX 50 // readfile and filename


struct fsm_{
	int type;
	int s;
	int i;
	int p;
	int succ[MAX_NB_STATES][NB_SUCC];
        int predsSizes[MAX_NB_STATES][NB_SUCC];
        int preds[MAX_NB_STATES][NB_SUCC][MAX_NB_STATES];
};


// create FSM with s states, i inputs and the transitions trans
FSM setFSM(int s, int i, int** trans) {

	FSM f = (FSM)malloc( sizeof(struct fsm_) + 2*sizeof(int) * s );

	f->type = 0;
	f->s = s;
	f->i = i;
	f->p = s*i;

	int j;
	for (j=0; j<s; j++)
            for(int k = 0; k < i; k++)
                f->predsSizes[j][k] = 0;
	
        for (j=0; j<s; j++){
            for(int k = 0; k < i; k++){
	        int succ = f->succ[j][k] = trans[j][k];
                f->preds[succ][k][f->predsSizes[succ][k]] = j;
                f->predsSizes[succ][k]++;
            }
	}   
	
	return f;
}


// free FSM
void freeFSM(FSM f) {
	free(f);
}


// read FSM from filename
FSM initFSM(const char * filename) {

	FILE * file;
	if ((file = fopen(filename, "r")) == NULL) {
		printf("Error opening file, exit.\n");
		exit(1);
	}

	char chaine[TAILLE_MAX] = "";

	while (fgets(chaine, TAILLE_MAX, file)[0] != 's') {}

	//int s = chaine[2] - '0';
	char *spt = strtok(chaine, " ");
	spt = strtok(NULL, " ");
	int s = atoi(spt);

	int i = fgets(chaine, TAILLE_MAX, file)[2] - '0';
	
	//printf("%d %d\n", s,i);

	while (fgets(chaine, TAILLE_MAX, file)[0] != 'p') {}
	
	int ** tr = (int**)malloc( s * sizeof(int[2]) );

	int j;
	for(j=0; j<s; j++) {
		tr[j] = (int*)malloc(2 * sizeof(int));
	}

	for (j=0; j<s; j++) {
		// inp 0
		fgets(chaine, TAILLE_MAX, file);
		//printf("%s", chaine);
		//tr[j][0] = chaine[4] - '0';
		spt = strtok(chaine, " ");
		spt = strtok(NULL, " ");
		spt = strtok(NULL, " ");
		tr[j][0] = atoi(spt);

		// inp 1	
		fgets(chaine, TAILLE_MAX, file);		
		//printf("%s\n", chaine);
		//tr[j][1] = chaine[4] - '0';
		spt = strtok(chaine, " ");
		spt = strtok(NULL, " ");
		spt = strtok(NULL, " ");
		tr[j][1] = atoi(spt);
	
	}	
	fclose(file);

	FSM f = setFSM(s,i,tr);

	for (j=0; j<s; j++) {
		free(tr[j]);
	}
	free(tr);

	return f;
}


// print a FSM
void printFSM(FSM f) {

	printf("\tPrint FSM\n");

	printf("F %d\n", f->type);
	printf("s %d\n", f->s);
	printf("i %d\n", f->i);
	printf("p %d\n", f->p);

	int j;
	for(j=0; j<f->s; j++)
	{
		printf("%d 0 %d 0\n", j,f->succ[j][0]);
		printf("%d 1 %d 0\n", j,f->succ[j][1]);
	}

	printf("\tEnd Print FSM\n\n");
}


// return nb of states s
int get_s(FSM f){
	return f->s;
}

// return nb of input i
int get_i(FSM f){
	return f->i;
}

// return successor of s after applying input i
int get_succ(FSM f, int s, int i) {
	return f->succ[s][i];
}

int get_nb_preds(FSM f, int s, int i){
    return f->predsSizes[s][i];
}

int get_pred(FSM f, int s, int i, int pred_index){
    return f->preds[s][i][pred_index];
}
