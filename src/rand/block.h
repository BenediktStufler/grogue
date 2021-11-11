/*
 * Simulate a random graph from a block-stable class
 */
int grog(struct cmdarg *comarg, gsl_rng **rgens) {
	INT *degprofile = NULL;		// outdeg profile of tree
	INT *D = NULL;				// degree sequence of tree
	INT *Gdegprofile = NULL;		// degree profile of graph
	INT counter;
	char *cname;
	struct graph *G = NULL;
	INT i,j;

	INT ***cdeco = NULL;
	INT **cbeta = NULL;

	struct outerdec **odeco = NULL;

	INT max;


	//initialize precomputed values
	if(comarg->method == 2) {
		cactiinit();
	} else if(comarg->method == 3) {
		outerinit();
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
		} else if(comarg->method == 3) {	// class of outerplanar graphs
			/* generate degree profile */
			degprofile = tbinb_prec(comarg->size, comarg->size-1, outerq, outermax, comarg->threads, rgens);

			// set maximum
			for(i=0, max=0; i<comarg->size; i++) 
				if(degprofile[i]>0) max = i;
			
			// get deco
			odeco = outerdeco(degprofile, max, comarg->threads, rgens);

			// debug
			/*
			printf("max = %ld\n", max);
			printf("N[%d]= %ld\n", 0, degprofile[0]);
			printf("N[%d]= %ld\n", 1, degprofile[1]);
			for(i=2;i<=max;i++)  {
				if(degprofile[i]>0) {
					printf("N[%ld]= %ld\n", i, degprofile[i]);
					for(j=0;j<degprofile[i];j++) {
						printf("components = %ld\n", odeco[i][j].numcomp);
						for(int k=0; k<2*i; k++) 
							printf("%ld,", odeco[i][j].stor[k]);
						printf("\n");
					}
				}
			}
			*/



			/* generate degree sequence */
			D = gendegsequence(degprofile, comarg->size, rgens[0]);

			/* cyclically shift sequence */
			cycshift(D, comarg->size);

			/* debug
			printf("D = ");
			for(i=0;i<comarg->size;i++)
				printf("%ld, ", D[i]);
			printf("\n");
			*/

			/* create graph */
			G = outerbij(D, odeco, comarg->size, max);

			/* free decoration */
			for(i=2; i<=outerprec; i++) {
				if( odeco[i] != NULL) free(odeco[i]);
			}
			for(i=outerprec+1; i<=max; i++) {
				for(j=0; j<degprofile[i]; j++) {
					free(odeco[i][j].stor);
				}
				if( odeco[i] != NULL) free(odeco[i]);
			}
			free(odeco);


		} else {
			// method not supported yet
			fprintf(stderr, "Chosen graph class is not supported yet.\n");
			exit(-1);
		}

		// shuffle ids?

		// output stuff
		if( comarg->Tdegfile || comarg->Toutfile  || comarg->Theightfile || comarg->Tcentfile || comarg->Tprofile || comarg->Tmaxheightfile || comarg->Tmaxdegfile ) {
			/* set height, vertex degrees, bfs order, disconnected warning flag */
			if(G->bfs == NULL) G->bfs = bfsorder(G, G->root, 1, 1);


			/* output graph if requested */
			if( comarg->Toutfile ) {
				cname = convname(comarg->outfile, counter, comarg->num, comarg->Aoutfile);
				outgraph(G,cname,comarg->Woutfile);
				free(cname);
			}

			/* output vertex degree profile if requested */
			if( comarg->Tprofile ) {
				cname = convname(comarg->profile, counter, comarg->num, comarg->Aprofile);
				// assumes that the deg parameters
				// have already been set
				Gdegprofile = makedegprofile(G);

				outdegprofile(Gdegprofile, comarg->size, cname, comarg->Wprofile);
				free(cname);
				free(Gdegprofile);
			}

			/* output degree sequence if requested */
			if( comarg->Tdegfile ) {
				cname = convname(comarg->degfile, counter, comarg->num, comarg->Adegfile);
				outdegseq(G, cname, comarg->Wdegfile);
				free(cname);
			}

			/* output maximum degree if requested */
			if( comarg->Tmaxdegfile ) {
				cname = convname(comarg->maxdegfile, counter, comarg->num, comarg->Amaxdegfile);
				outmaxdeg(G, cname, comarg->Wmaxdegfile);
				free(cname);
			}

			/* output height sequence if requested */
			if( comarg->Theightfile ) {
				cname = convname(comarg->heightfile, counter, comarg->num, comarg->Aheightfile);
				outheightseq(G, cname, comarg->Wheightfile);
				free(cname);
			}

			/* output maximal height if requested */
			if( comarg->Tmaxheightfile ) {
				cname = convname(comarg->maxheightfile, counter, comarg->num, comarg->Amaxheightfile);
				outmaxheight(G, cname, comarg->Wmaxheightfile);
				free(cname);
			}



			/* Calculate closeness centrality if requested */
			if( comarg->Tcentfile ) {
				cname = convname(comarg->centfile, counter, comarg->num, comarg->Acentfile);
				threadedcentrality(G, 0, G->num, comarg->threads);
				outcent(G, cname, comarg->Wcentfile);
				free(cname);
			}
		}

		// clean up
		if(G != NULL) free_graph(G);
		if(D != NULL) free(D);
		if(degprofile != NULL) free(degprofile);
	}


	// free precomputed values... (add function)

	return 0;
}


