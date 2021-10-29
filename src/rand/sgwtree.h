/*
 * generate degree sequence from degree profile
 */
INT *gendegsequence(INT *N, INT size, gsl_rng *rgen) {
	INT i, j, p;	
	INT *out;

	// result will be stored in an array of integers
	out = (INT *) calloc(size, sizeof(INT));
	if(out == NULL) {
		// memory allocation error
		fprintf(stderr, "Memory allocation error in function gendegsequence\n");
		exit(-1);
	}

	// create initial array with sequential list of outdegrees
	for(i=0, p=0; i<size; i++)
		for(j=0; j<N[i]; j++, p++)
			out[p]=i;

	// shuffle 
	gsl_ran_shuffle (rgen, out, size, sizeof (INT));

	return out;
}

/*
 * Cyclically shifts a vector to make it a degree sequence of a tree
 */

int cycshift(INT *D, INT size) {
	INT i, j;
	INT indmin;
	long long sum, min;
	INT *tmp;

	for(i=0, sum=0, min=0, indmin=0; i<size; i++) {	
		sum += (long long) D[i] - 1;
		if(sum < min) {
			min = sum;
			indmin = i;
		}
	}

	if(indmin < size-1) {
		tmp = (INT *) calloc(size, sizeof(INT));
		if(tmp == NULL) {
			// memory allocation error
			fprintf(stderr, "Memory allocation error in function cycshift\n");
			exit(-1);
		}	
		
		// save reordered list in temporary array
		for(i=indmin+1,j=0; i<size; i++,j++)
			tmp[j] = D[i];
		for(i=0; i<indmin+1; i++,j++)
			tmp[j] = D[i];  // j has now the correct value, despite compiler warnings

		// copy back
		for(i=0; i<size; i++)
			D[i] = tmp[i];

		// free memory
		free(tmp);
	}

	return 0;
}
				
/*
 * Compute tree from outdegree sequence - bfs order
 */
struct graph *deg2bfstree(INT *D, INT len) {
	struct graph *G;
	INT pos, i, j;

	// initialize graph
   	G = newgraph(len); 
	
	// add edges such that the i th entry of the outdegree list corresponds
	// to the outdegree of the i th vertex in bfs order and such that
	// the first vertex in the neighbourhood list of any non-root vertex is 
	// its parent

	if(len>0) {
		G->arr[0]->height = 0;	
		G->root = G->arr[0];
	}
	for(i=0, pos=1; i<len; i++) {
		// save vertex degree for later use
		G->arr[i]->deg = D[i];
		for(j=0; j<D[i]; j++, pos++) {
			addEdge(G->arr[i], G->arr[pos]);
			// save vertex height for later use
			G->arr[pos]->height = G->arr[i]->height + 1;
		}
	}	
	

	return G;
}



/*
 * Compute tree from outdegree sequence - dfs order
 */
struct graph *deg2dfstree(INT *D, INT len) {
	struct graph *G;
	INT i;

	G = deg2bfstree(D, len);	// construct tree in BFS-order
								// sets deg and height attributes
	G->bfs = G->arr;			// set bfs-order list
	G->arr = dfsorder(G, G->root); // calculate DFS-order
	
	// change vertex ids to dfs order
	for(i=0; i<G->num; i++)
		G->arr[i]->id = i;

	return G;
}




