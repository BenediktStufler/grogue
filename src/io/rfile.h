
/*
 * Read connected graph from file instead of generating the graph at random
 */
int rfile(struct cmdarg *comarg) {
	struct graph *G;
	INT *degprofile;




	/* read graph from file */
	G = parsegraphml(comarg->infile, comarg->vid);

	if(G->num == 0) return 0;	 // nothing to do if there are no vertices

	if(G->root == NULL) G->root = G->arr[0]; // set root if none was specified

	/* set height, vertex degrees, bfs order, disconnected warning flag */
	G->bfs = bfsorder(G, G->root, 1, 1);
	if(G->disconnected) {
		fprintf(stderr, "Error: graph from input file is disconnected. Disconnected graphs are not supported at the moment.\n"); 
		return -1;
	}

	
	/* output degree profile if requested */
	if( comarg->Tprofile ) {
		degprofile = makedegprofile(G);	// assumes that the deg parameters
										// have already been set
		outdegprofile(degprofile, G->num, comarg->profile, 'w');
		free(degprofile);
	}
	
	/* output degree sequence if requested */
	if( comarg->Tdegfile ) {
		outdegseq(G, comarg->degfile, 'w');
	}

	/* output maximum degree if requested */
	if( comarg->Tmaxdegfile ) {
		outmaxdeg(G, comarg->maxdegfile, 'w');
	}

	/* output height sequence if requested */
	if( comarg->Theightfile ) {
		outheightseq(G, comarg->heightfile, 'w');
	}

	/* output maximal height if requested */
	if( comarg->Tmaxheightfile ) {
		outmaxheight(G, comarg->maxheightfile, 'w');
	}


	/* Calculate closeness centrality if requested */
	if( comarg->Tcentfile ) {
		threadedcentrality(G, 0, G->num, comarg->threads);
		outcent(G, comarg->centfile, 'w');
	}



	/* clean up */
	free_graph(G);

	return 0;
}
