/*
 * data that gets passed to a thread
 */
struct cactiarg {
	INT ***deco;
        gsl_rng *rgen;
        INT *counter;
        pthread_mutex_t *mut;
        INT max;
	INT *degprofile;
	INT *dcounter;
};			


/*
 * Generate cacti graph from decorated tree
 */
struct graph *cactibij(INT *D, INT **cbeta, INT len) {
	struct graph *G;	
	INT pos, i, j, k;

	if(len <= 0) {
		fprintf(stderr, "Error in function cactibij. Length value <= 0.\n");
		exit(-1);
	}

	// initialize graph
   	G = newgraph(len); 
	G->root = G->arr[0];

	// debug
	/*
	for(i=0;i<len;i++) {
		if( D[i] == 0 ) {
			printf("i=%ld, D[i] =%ld, cbeta[i]=NULL\n", i, D[i]);
		} else {
			printf("i=%ld, D[i] =%ld, cbeta[i]=", i, D[i]);
			for(j=0;j<=D[i];j++)
				printf("%ld,",cbeta[i][j]);
			printf("\n");
		}
	}
	*/


	for(i=0, pos=1; i<len; i++) {
		if( D[i] > 0 ) {
			for(j=1; j<=cbeta[i][0]; j++) {
				if( cbeta[i][j] == 1 ) {
					addEdge(G->arr[i], G->arr[pos]);
					pos++;
				} else {
					addEdge(G->arr[i], G->arr[pos]);
					for(k=1; k<cbeta[i][j]; k++)  {
						addEdge(G->arr[pos], G->arr[pos+1]);
						pos++;
					}
					addEdge(G->arr[i], G->arr[pos++]);
				}
			}
		}
	}

/*
	for(i=0, pos=1; i<len; i++) {
		// save vertex degree for later use
		for(j=0; j<D[i]; j++, pos++) {
			addEdge(G->arr[i], G->arr[pos]);
			// save vertex height for later use
			G->arr[pos]->height = G->arr[i]->height + 1;
		}
	}	
*/
	return G;

}



/* 
 * Create list of decorations 
 */
INT **cactibeta(INT *D, INT ***cdeco, INT len, INT max) {
	INT i;
	INT **cactibeta;
	INT *counter;
	
       	cactibeta = (INT **) calloc(len, sizeof(INT *));
       	counter = (INT *) calloc(max+1, sizeof(INT));	// initializes to zero
	if(cactibeta == NULL || counter == NULL) {
		// memory allocation error
		fprintf(stderr, "Memory allocation error in function cactibeta\n");
		exit(-1);
	}

	for(i=0;i<len;i++) {
		if(D[i]==0) {
			cactibeta[i] = NULL;
		} else {
			//debug
			//printf("cactibeta[%ld] = cdeco[%ld][%ld]\n", i, D[i], counter[D[i]]);

			cactibeta[i] = cdeco[D[i]][counter[D[i]]];
			counter[D[i]]++;
		}
	}

	free(counter);

	return cactibeta;
}



/*
 * Generate blocks for cacti graph decoration
 */
void *cactiblock(void *dim) {
	// retrieve data that got passed to thread
	INT max=((struct cactiarg *)dim)->max;
	INT *degprofile=((struct cactiarg *)dim)->degprofile;
	INT *counter=((struct cactiarg *)dim)->counter;
	INT ***deco=((struct cactiarg *)dim)->deco;

	gsl_rng *rgen = ((struct cactiarg *)dim)->rgen;
	pthread_mutex_t *mut = ((struct cactiarg *)dim)->mut;

	INT i;
	INT size = 0;

	INT N;
	INT *K;
	INT *dcounter = ((struct cactiarg *)dim)->dcounter;


	/* calls to calloc need to be protected by mutex */
	pthread_mutex_lock(mut);
	K = (INT *) calloc(max, sizeof(INT));
	if(K == NULL) {
		// memory allocation error
		fprintf(stderr, "Memory allocation error in function cactiblock\n");
		exit(-1);
	}
	pthread_mutex_unlock(mut);

	/*
	 * B'(z) = z + z^2/(2*(1-z))
	 *
	 * t an arbitrary constant satisfying C^bullet(rho_C) < t < rho_B
	 * 	with C^bullet(rho_C) = 0.4563... and rho_B = 1
	 * N = Poisson( B'(t))
	 * then for each 0 <= i < N:
	 * if( coin_flip( t/B'(t) ) ) then K_i = 1
	 * else K_i = 2 + Geom(t)
	 *
	 */

	/* Mathematica code:
	Bp = z + z^2/(2*(1-z))
	t = 0.95
	lambda = N[ Bp /.{z->t}, 200]
	p = N[ (z/Bp) /.{z->t}, 200]
	nu = N[1-t, 200]
	Bpp = D[Bp, z]
	exsize = N[ z*Bpp  /.{z->t}, 200]
	
		                2
         	z (2 - 2 z + z )
	exsize= ----------------
		          2
           	2 (-1 + z)
	*/

	/* choose value for parameter t that optimizes target size */


	/*double t = 0.8;
	 * --> target size about 19.3139*/
	double lambda = 2.4;
	double p = (double) 1.0 / (double) 3.0;
	double nu = 0.2;

	/*double t = 0.9;
	 * --> target size about 45.45
	double lambda = 4.95;
	double p = 0.181818;
	double nu = 0.1;
	*/

	/*double t = 0.92;
	 * --> target size about 72.335
	double lambda = 0.4968;
	double p = 0.148148;
	double nu = 0.4968;
	*/


	/* double t = 0.95;
	 * --> target size about 190.475
	double lambda = 9.975;	// B'(t)
	double p = 0.0952381;	// t / B'(t)
	double nu = 0.05;	// 1-t
	*/

	/* double t = 0.99;
	 * --> target size about 4950.49
	double lambda = 49.995;
	double p = 0.019802;
	double nu = 0.01;
	*/

	

	while(*counter > 0) {
		// take sample with size smaller than max
		N = gsl_ran_poisson(rgen, lambda);
		size = N;

		//debug
		//printf("max: %ld, poisson value: %ld\n",max, N);

		if(size>max) continue;
		for(i=0;i<N; i++) {
			if( gsl_ran_bernoulli(rgen, p) == 1 ) {
				K[i] = 1;
			} else {
				// beware that gsl_ran_geometric(nu) has density
				// p_k = nu (1-nu)^{k-1}, k>=1
				K[i] = 1 + gsl_ran_geometric(rgen, nu);
				size += K[i]-1;
			}
		}
		// debug
		// printf("poisson=%ld, size=%ld\n", N, size);

		if(size>max || size<=0) continue;	// we ignore structures of size 0

		// check counter
		if(dcounter[size] < degprofile[size]) {
			// yay, we found a valid decoration

			/* begin of part that is partially locked by mutex */
			pthread_mutex_lock(mut);
			// check again if we still need to do this
			if(dcounter[size] < degprofile[size]) {		

				//debug
				//printf("counter: %ld, size: %ld, dcounter[size] = %ld, degprofile[size] = %ld\n", *counter, size, dcounter[size], degprofile[size]);

				// save solution 
				deco[size][dcounter[size]][0] = N;
				for(i=1;i<=N;i++)
					deco[size][dcounter[size]][i] = K[i-1];

				// update counters
				dcounter[size] += 1;	
				*counter -= 1;

				//debug
				//if(dcounter[size]==degprofile[size]) printf("%ld size complete, all of %ld samples\n", size, degprofile[size]);
			}
			// unlock mutex - there's still work to do
			pthread_mutex_unlock(mut);
		}
	}

	// we found enough configurations 
				
	// time to clean up...
	free(K);

	// ... and go home
	return (void *) 0;
}




/*
 * Generate decorations for cacti graphs
 */
INT ***cactideco(INT *degprofile, INT max, int numThreads, gsl_rng **rgens) {
	INT i,j;
	struct cactiarg *argList;   	// arguments for the separate threads
	pthread_t *th;			// array of threads
	void *ret;
	INT *dcounter;

	// mutex for thread synchronization
	pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;


	// sanity check
	if(max <= 0) {
		// nothing to decorate here
		fprintf(stderr, "Error in function cactideco. Nothing to decorate here.\n");
		exit(-1);
	}


	// allocate space for pointers to graphs
	INT ***deco = calloc(sizeof(INT **), max+1);
	if(deco == NULL) {
		// memory allocation error
		fprintf(stderr, "Memory allocation error in function cactideco\n");
		exit(-1);
	}
	deco[0] = NULL;	// no graphs of size zero
	for(i=1; i<=max; i++) {
		if(degprofile[i] >0) {
			deco[i] = calloc(sizeof(INT *), degprofile[i]);
			if(deco[i] == NULL) {
				// memory allocation error
				fprintf(stderr, "Memory allocation error in function cactideco\n");
				exit(-1);
			}
			for(j=0; j<degprofile[i]; j++) {
				deco[i][j] = calloc(sizeof(INT), i+1);
				if(deco[i][j] == NULL) {
					// memory allocation error
					fprintf(stderr, "Memory allocation error in function cactideco\n");
					exit(-1);
				}
			}
		} else {
			deco[i] = NULL;
		}
	}




	/* optimization for small values */
	// trivial for size = 1
	for(j=0; j<degprofile[1]; j++) {
		deco[1][j][0] = 1;
		deco[1][j][1] = 1;
	}

	INT choice;
	for(i=2; i<=cactiprec && i<=max; i++) {
		for(j=0; j<degprofile[i]; j++) {
			//debug
			//printf("%ld %ld\n",i,j);

			choice = gsl_ran_discrete(rgens[0],cactigen[i]);
			memcpy(deco[i][j],cactipart[i][choice],(cactipart[i][choice][0]+1)*sizeof(INT));
		}
	}


	// the variable counter will only be read / modified in an area 
	// protected by the mutex mut
	INT counter=0;
	for(i=cactiprec+1; i<=max; i++)
		counter += degprofile[i];


	if(counter > 0) {
		INT *truncdegprofile = calloc(sizeof(INT), max+1); // initializes to zero
		if(truncdegprofile  == NULL) {
			// memory allocation error
			fprintf(stderr, "Memory allocation error in function cactideco\n");
			exit(-1);
		}
		for(i=cactiprec+1; i<=max; i++)
			truncdegprofile[i] = degprofile[i];

		dcounter = (INT *) calloc(max+1, sizeof(INT));	// initializes to zero
		if(dcounter == NULL) {
			// memory allocation error
			fprintf(stderr, "Memory allocation error in function balls in boxes\n");
			exit(-1);
		}



		// pack list of arguments
		argList = calloc(sizeof(struct targ), numThreads);
		for(i=0; i<numThreads; i++) {
			argList[i].deco = deco;
			argList[i].max = max;
			argList[i].degprofile = truncdegprofile;
			argList[i].rgen = rgens[i];
			argList[i].counter = &counter;
			argList[i].dcounter = dcounter;
			argList[i].mut = &mut;
		}



		/* launch threads */
		th = calloc(sizeof(pthread_t), numThreads);
		for(i=0; i<numThreads; i++) {
			if(pthread_create(&th[i], NULL, &cactiblock, &argList[i] )) {
				fprintf(stderr, "Error launching thread number %"STR(FINT)"\n", i);
				exit(-1);
			}
		}


		/* wait for threads to finish */
		for(i=0; i<numThreads; i++) {
			pthread_join(th[i], &ret);
			if(ret) {
				fprintf(stderr, "Error executing thread number %"STR(FINT)"\n", i);
				exit(-1);
			}
		}


		/* clean up */
		free(argList);
		free(th);
		free(truncdegprofile);
		free(dcounter);
	}

	/* clean up */
	pthread_mutex_destroy(&mut);

	return deco;	
}


