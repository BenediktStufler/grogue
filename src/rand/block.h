/*
 * Simulate a random graph from a block-stable class
 */
int grog(struct cmdarg *comarg, gsl_rng **rgens) {
	INT *degprofile;		// outdeg profile of tree
	INT *D;				// degree sequence of tree
	INT *Gdegprofile;		// degree profile of graph
	INT counter;
	char *cname;
	struct graph *G;
	INT i,j;

	INT ***cdeco;
	INT **cbeta;
	INT max;


	//debug
	G=NULL;

	//initialize precomputed values
	if(comarg->method == 2) {
		cactiinit();
	}


	for(counter=1; counter <= comarg->num; counter++) {	
		if(comarg->method == 1) {	// class of trees
			/* generate degree profile */
			degprofile = binbpoisson(comarg->size, comarg->size-1, rgens);

			/* generate degree sequence */
			D = gendegsequence(degprofile, comarg->size, rgens[0]);

			/* cyclically shift sequence */
			cycshift(D, comarg->size);

			// already sets the height and deg attribute of each vertex
			// list of vertices is now in dfs order, bfs order is set
			G = deg2dfstree(D, comarg->size);

			// deg attribute is set to outdegree; need to set it to degree
			// root corresponds to i=0 and is skipped on purpose
			for(i=1; i<G->num; i++)
				G->arr[i]->deg = G->arr[i]->deg + 1;

		} else if(comarg->method == 2) {	// class of cacti graphs
			/* generate degree profile */
			degprofile = tbinb_prec(comarg->size, comarg->size-1, cactiq, cactimax, comarg->threads, rgens);

			/*debug*/
			//printf("Done getting degree list\n");

			/* generate decorations */
			// set maximum
			for(i=0, max=0; i<comarg->size; i++) 
				if(degprofile[i]>0) max = i;

			//debug
			//printf("%ld\n", max);

			// get deco
			cdeco = cactideco(degprofile, max, comarg->threads, rgens);
	
			// debug	
			/*
			for(i=0;i<=max;i++)
				printf("N[%ld]= %ld\n", i, degprofile[i]);
			for(i=1; i<=max; i++) {
				for(j=0; j<degprofile[i]; j++) {
					printf("%ld - %ld,", i, cdeco[i][j][0]);
					for(int k = 1; k<= cdeco[i][j][0]; k++) {
						printf("%ld,",cdeco[i][j][k]);
					}
					printf("\n");
				}
			}
			*/

			/* generate degree sequence */
			D = gendegsequence(degprofile, comarg->size, rgens[0]);

			/* cyclically shift sequence */
			cycshift(D, comarg->size);

			/* create list of decorations */
			cbeta = cactibeta(D, cdeco, comarg->size, max);

			/* create graph */
			G = cactibij(D, cbeta, comarg->size);

			/* free decoration */
			for(i=1; i<=max; i++) {
				for(j=0; j<degprofile[i]; j++) {
					free(cdeco[i][j]);
				}
				if( cdeco[i] != NULL) free(cdeco[i]);
			}
			free(cdeco);
			free(cbeta);
		} else {
			// method not supported yet
			fprintf(stderr, "Chosen graph class is not supported yet.\n");
			exit(-1);
		}

		// shuffle ids?

		// output stuff
		if( comarg->Tdegfile || comarg->Toutfile  || comarg->Theightfile || comarg->Tcentfile ) {
			/* set height, vertex degrees, bfs order, disconnected warning flag */
			if(G->bfs == NULL) G->bfs = bfsorder(G, G->root, 1, 1);


			/* output graph if requested */
			if( comarg->Toutfile ) {
				cname = convname(comarg->outfile, counter, comarg->num, comarg->Tnum);
				outgraph(G,cname);
				free(cname);
			}

			/* output vertex degree profile if requested */
			if( comarg->Tprofile ) {
				cname = convname(comarg->profile, counter, comarg->num, comarg->Tnum);
				// assumes that the deg parameters
				// have already been set
				Gdegprofile = makedegprofile(G);

				outdegprofile(Gdegprofile, comarg->size, cname);
				free(cname);
				free(Gdegprofile);
			}

			/* output degree sequence if requested */
			if( comarg->Tdegfile ) {
				cname = convname(comarg->degfile, counter, comarg->num, comarg->Tnum);
				outdegseq(G, cname);
				free(cname);
			}
			
			/* output height sequence if requested */
			if( comarg->Theightfile ) {
				cname = convname(comarg->heightfile, counter, comarg->num, comarg->Tnum);
				outheightseq(G, cname);
				free(cname);
			}

			/* Calculate closeness centrality if requested */
			if( comarg->Tcentfile ) {
				cname = convname(comarg->centfile, counter, comarg->num, comarg->Tnum);
				threadedcentrality(G, 0, G->num, comarg->threads);
				outcent(G, cname);
				free(cname);
			}
		}

		// clean up
		if(G != NULL) free_graph(G);
		free(D);
		free(degprofile);
	}


	// free precomputed values... (add function)

	return 0;
}


