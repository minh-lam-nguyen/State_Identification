#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NB_SUCC 2 // nb of successors
#define SIZE 100 //
#define SIZE_Q 20000
#define SIZE_V 20000
#define TAILLE_MAX 20 // readfile

typedef struct {
	int type;
	int s;
	int i;
	int p;
	int succ[][NB_SUCC];
} FSM;


// create FSM with s states, i inputs and the transitions trans
FSM * setFSM(FSM *fsm, int s, int i, int** trans) {

	fsm = malloc( sizeof(FSM) + 2*sizeof(int) * s );

	fsm->type = 0;
	fsm->s = s;
	fsm->i = i;
	fsm->p = s*i;

	int j;
	for (j=0; j<s; j++)
	{
		fsm->succ[j][0] = trans[j][0];
		fsm->succ[j][1] = trans[j][1];
	}
	
	return fsm;

}


// free FSM
void freeFSM(FSM *fsm) {
	int i;
	for (i=0; i<fsm->s; i++) {
		free(fsm->succ[i]);
	}
	free(fsm);
}


// read FSM from filename
FSM * readFSM(const char * filename) {

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
	
	int ** tr = malloc( s * sizeof(int[2]) );

	int j;
	for(j=0; j<s; j++) {
		tr[j] = malloc(2 * sizeof(int));
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

	FSM * fsm = setFSM(fsm,s,i,tr);

	for (j=0; j<s; j++) {
		free(tr[j]);
	}
	free(tr);

	return fsm;
}


// print a FSM
void printFSM(FSM *fsm) {

	printf("\tPrint FSM\n");

	printf("F %d\n", fsm->type);
	printf("s %d\n", fsm->s);
	printf("i %d\n", fsm->i);
	printf("p %d\n", fsm->p);

	int j;
	for(j=0; j<fsm->s; j++)
	{
		printf("%d 0 %d 0\n", j,fsm->succ[j][0]);
		printf("%d 1 %d 0\n", j,fsm->succ[j][1]);
	}

	printf("\tEnd Print FSM\n\n");

}

void clearQV(int ** q, int ** v, int taille) {
	printf("d");
	int i;
	printf("m");
	for (i=0; i<taille-1; i++) {
		free(v[i]);
		free(q[i]);
		printf("%d %d \n",i, taille);
	}

	free(q);
	free(v);
	printf("f\n");
}


// Synchronizing tree
int syncTree(FSM *fsm, int* res) {

	// queue
	int **queue = malloc( SIZE_Q * sizeof(int*) );
	int front_q = 0;
	int count_q = 1;

	//int ttl = 1;

	// visited set of states
	int **visited = malloc( SIZE_V * sizeof(int*) );
	int count_v = 1;

	//if (queue && visited) {
	//	int i;
	//	for(i=0; i<SIZE_Q; i++) {
	//		queue[i] = malloc( (SIZE + 1 + fsm->s) * sizeof(int) );
	//		visited[i] = malloc( fsm->s * sizeof(int) );
	//	}
	//}

	if (!(queue && visited)) {
		return -1;
	}

	// add set of all state (init) on queue and visited
	queue[0] = malloc( (SIZE + 1 + fsm->s) * sizeof(int) );
	visited[0] = malloc( fsm->s * sizeof(int) );
	int j;
	for(j=0; j<fsm->s; j++) {
		queue[0][j] = 1;
		visited[0][j] = 1;
	}
	queue[0][j] = 0;

	// bfs
	while (count_q) {


		int current = front_q;
		front_q++; 
		count_q--;
		//printf("nb of visited: %d\n", count_v);
		//printf("current: %d\n", current);

		// c0 and c1 successors of current
		int *c0 = malloc(fsm->s * sizeof(int));
		int *c1 = malloc(fsm->s * sizeof(int));
		
		int j;
		for (j=0; j<fsm->s; j++) {
			c0[j] = 0;
			c1[j] = 0;
		}

		for (j=0; j<fsm->s; j++) {
			if (queue[current][j]) {
				c0[ fsm->succ[j][0] ] = 1;
				c1[ fsm->succ[j][1] ] = 1;
			}
		}


		// check if c0/c1 already visited
		int isVisited0 = 0;
		int isVisited1 = 0;
		
//		j = 0; 
//		while (!isVisited0 && !isVisited1 && j<count_v){
//			int k;
//			isVisited0 =1;
//			isVisited1 =1;
//
//			for(k=0; k< fsm->s;k++){
//				if (visited[j][k] != c0[k]){
//					isVisited0 =0;
//				}
//				if (visited[j][k] != c1[k]){
//					isVisited1 =0;
//				}
//			}
//
//			j++;
//		}

		j = 0; 
		while (!isVisited0 && j<count_v){
			int k;
			isVisited0 =1;

			for(k=0; k< fsm->s;k++){
				if (visited[j][k] != c0[k]){
					isVisited0 =0;
				}
			}

			j++;
		}

		j = 0; 
		while (!isVisited1 && j<count_v){
			int k;
			isVisited1 =1;

			for(k=0; k< fsm->s;k++){
				if (visited[j][k] != c1[k]){
					isVisited1 =0;
				}
			}

			j++;
		}


	

		// size of current seq
		int seq_size = queue[current][fsm->s];
		//printf("%d\n", seq_size);

		// c0 not in visited
		if (!isVisited0) {
			int k, tmp_cpt = 0;
			for(k=0; k<fsm->s; k++) {
				if (c0[k] == 1) {
					tmp_cpt++;
				}
			}

			// if singleton (1 and only one state) : SS
			if (tmp_cpt==1) {
				for(k=0; k<seq_size ;k++) {
					res[k] = queue[current][k+fsm->s+1];
				}
				res[k] = 0;
				//clearQV(queue, visited, count_v);
				return seq_size+1;
			}else{
			// else : add to visited and queue
				count_q++;
				queue[current+count_q] = malloc( (SIZE + 1 + fsm->s) * sizeof(int) );
				visited[count_v] = malloc( fsm->s * sizeof(int) );

				for(k=0; k<fsm->s; k++) {
					visited[count_v][k] = c0[k];
					queue[current+count_q][k] = c0[k];
				}
				count_v++;
				queue[current+count_q][k] = seq_size+1;


				for(k=0; k<seq_size; k++) {
					queue[current+count_q][k+fsm->s+1] = queue[current][k+fsm->s+1];
				}
				queue[current+count_q][k + fsm->s + 1] = 0;
				//count_q++;
				//ttl++;

			}
		}

		// c1 not in visited
		if (!isVisited1) {
			// if singleton (1 and only one state)
			int k, tmp_cpt = 0;
			for(k=0; k<fsm->s; k++) {
				if (c1[k] == 1) {
					tmp_cpt++;
				}
			}
			
			// if singleton (1 and only one state) : SS
			if (tmp_cpt==1) {
				for(k=0; k<seq_size ;k++) {
					res[k] = queue[current][k+fsm->s+1];
				}
				res[k] = 1;
				//clearQV(queue, visited, count_v);
				return seq_size+1;
			}else{
			// else : add to visited and queue
				count_q++;
				queue[current+count_q] = malloc( (SIZE + 1 + fsm->s) * sizeof(int) );
				visited[count_v] = malloc( fsm->s * sizeof(int) );

				for(k=0; k<fsm->s; k++) {
					visited[count_v][k] = c1[k];
					queue[current+count_q][k] = c1[k];
				}
				count_v++;
				queue[current+count_q][k] = seq_size+1;

				for(k=0; k<seq_size; k++) {
					queue[current+count_q][k+fsm->s+1] = queue[current][k+fsm->s+1];
				}
				queue[current+count_q][k + fsm->s + 1] = 1;
				//count_q++;
				//ttl++;
			}
		}



		//free(queue[current]);

		free(c0);
		free(c1);
	}

	//clearQV(queue, visited, current, count_v);

	return 0;

}



int main() {

/*
	int s = 4;
	int i = 2;
	int tr[4][2] = { {1,3}, {0,1}, {3,0}, {3,2} };
	//int tr[4][2] = { {3,3}, {3,1}, {3,0}, {3,2} };
	//int tr[4][2] = { {1,1}, {0,1}, {3,1}, {3,1} };
	
	
	FSM* fsm = setFSM(fsm,s,i,tr);
	printFSM(fsm);


	time_t begin = time(NULL);

	int res[SIZE];
	int taille = syncTree(fsm, res);

	//sleep (2);
	time_t end = time (NULL);

	printf("time: %f seconds \n", difftime(end,begin));

	printf("length SS: %d\n", taille);
	int j;
	printf("SS: ");
	for (j=0; j<taille; j++){
		printf("%d", res[j]);
	}
	printf("\n");
*/

/*
	FSM * fsm = readFSM("fsm_hss.fsm");

	printFSM(fsm);

	time_t begin = time(NULL);

	int res[SIZE];
	int taille = syncTree(fsm, res);

	time_t end = time (NULL);

	printf("time: %f seconds \n", difftime(end,begin));

	printf("length SS: %d\n", taille);
	int j;
	printf("SS: ");
	for (j=0; j<taille; j++){
		printf("%d", res[j]);
	}
	printf("\n\n");

	freeFSM(fsm);

*/

	int i;
	//float total = 0, ttl = 0;
	clock_t beg = clock();

	for (i=1; i<25; i++) {
		char tmp[SIZE] = "SS_50fsm_n30/fsm_n30_";
		char numb[10];
		strcat( tmp, itoa(i, numb, 10) );
		//strcat( tmp, itoa(12, numb, 10) ); // longest SS
		//strcat( tmp, itoa(34, numb, 10) ); //no SS
		strcat( tmp, ".fsm" );

		FSM * fsm = readFSM(tmp);
		//printFSM(fsm);

		//time_t begin = time(NULL);
		//clock_t b = clock();


		int res[SIZE];
		//syncTree(fsm, res);
		int taille = syncTree(fsm, res);

		//time_t end = time(NULL);
		//clock_t e = clock();

		//total += difftime(end,begin);
		//ttl += (double)(e-b);

		//printf("dfa %d: length SS: %d time: %f s / %f ms\n", i, taille, difftime(end, begin), (double)(e-b));
		printf( "dfa %d: length SS: %d\n", i, taille );
		
		//int j;
		//printf("SS: ");
		//for (j=0; j<taille; j++){
		//	printf("%d", res[j]);
		//}
		//printf("\n\n");

		freeFSM(fsm);

	}

	//printf("\ntotal time: %f ms / average (50 dfa): %f ms\n", ttl, ttl/(i-1));
	//printf("total time: %f seconds / average (50 dfa): %f seconds\n", total, total/(i-1));

	clock_t end = clock();
	double ttl = (double)(end-beg);
	printf("total time: %f ms / average: %f ms (%d fsm)", ttl, ttl/(i-1), i-1);

	return 0;

}