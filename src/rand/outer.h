/*
 * data that gets passed to a thread
 */

struct outerdec {
	INT numcomp;
	INT *stor;
};

struct outeredge {
	struct vertex *start;
	struct vertex *end;
};

struct outerarg {
	struct outerdec **deco;
        gsl_rng *rgen;
        INT *counter;
        pthread_mutex_t *mut;
        INT max;
	INT *degprofile;
	INT *dcounter;
	INT tid;		// thread id
};			


/*
 * Generate outerplanar graph from decorated tree
 */
struct graph *outerbij(INT *D, struct outerdec **odec, INT len, INT max) {
	struct graph *G;	
	struct vertex **pos;
	struct vertex *start;
	struct vertex *end;
	INT i, j, l;
	long int e; 	// may be assigned negative values
			// We use index of array instead of pointer because pointing
			// one below start of array is undefined behaviour.
			// (Side note: pointing one above end of array is allowed.)
	INT *cyc;
	struct outerdec *dec;

	if(len <= 0) {
		fprintf(stderr, "Error in function outerbij. Length value <= 0.\n");
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

	INT *v = calloc(max+1, sizeof(INT));	// initializes to zero
	struct outeredge *estack = malloc(sizeof(struct outeredge)*2*max);
	if(v == NULL || estack == NULL) {
		// memory allocation error
		fprintf(stderr, "Memory allocation error in function outerbij\n");
		exit(-1);
	}


	for(i=0, pos=G->arr+1; i<len; i++) {
		if( D[i] > 0 ) {
			if(D[i]==1) {
				addEdge(G->arr[i], *pos);
				pos++;
			} else {
				// pop decoration of correct size
				dec = &odec[D[i]][v[D[i]]];
				cyc = dec->stor;
				v[D[i]]++;
				
				for(j=0;j<dec->numcomp;j++) {
					// init stack
					// push initial edge
					estack[0].start = G->arr[i];
					estack[0].end = *pos;
					e=0;

					// add initial edge to graph
					addEdge(G->arr[i], *pos);
					pos++;

					// build dissection
					while(e >= 0) {
						if(*cyc==0) {
							// pop edge (abridged)
							e--;
						} else {
							// add path of length *cyc>=2
							// from start to end
						
							// pop edge and push new edge in one go
							start = estack[e].start;
							end = estack[e].end;
							estack[e].end = *pos;
							addEdge(start, *pos);

							for(l=1; l<*cyc-1;l++) {
								// push new edge
								e++;
								estack[e].start = *pos;
								estack[e].end = *(pos+1);

								addEdge(*pos, *(pos+1));
								pos++;
							}
							e++;
							estack[e].start = *pos;
							estack[e].end = end;

							addEdge(*pos, end);
							pos++;
						}
						cyc++;
					}
				}
			}
		}
	}

	free(v);
	free(estack);
	return G;

}




/*
 * Generate blocks for outer graph decoration
 */
void *outerblock(void *dim) {
	// retrieve data that got passed to thread
	INT max=((struct outerarg *)dim)->max;
	INT *degprofile=((struct outerarg *)dim)->degprofile;
	INT *counter=((struct outerarg *)dim)->counter;
	struct outerdec **deco=((struct outerarg *)dim)->deco;

	gsl_rng *rgen = ((struct outerarg *)dim)->rgen;
	pthread_mutex_t *mut = ((struct outerarg *)dim)->mut;
	//INT tid = ((struct outerarg *)dim)->tid;

	long int i;
	INT size = 0;

	INT N;
	INT *K;
	INT *dcounter = ((struct outerarg *)dim)->dcounter;
	INT v;
	long int trec;	// may be assigned negative values

	/* calls to calloc need to be protected by mutex */
	pthread_mutex_lock(mut);
	K = (INT *) calloc(2*max, sizeof(INT));
	if(K == NULL) {
		// memory allocation error
		fprintf(stderr, "Memory allocation error in function outerblock\n");
		exit(-1);
	}
	pthread_mutex_unlock(mut);

	/*
	 * 2*B'(z) = z + D
	 * D = z + D^2/(1-D)
	 *
	 *
	 * t an arbitrary constant satisfying C^bullet(rho_C) < t < rho_B
	 * 	with C^bullet(rho_C) = 0.170765... and rho_B = 3 - 2 Sqrt[2] \approx 0.17157333...
	 * N = Poisson( B'(t))
	 * then for each 0 <= i < N:
	 * if( coin_flip( t/(2*B'(t)) ) then 
	 * 	W_i = [0]
	 * 	K_i = 1
	 * else
	 * 	W_i = outdegree sequencs of GWT
	 * 	K_i = number of leaves of this GWT
	 *
	 * 	offspring distribution:
	 * 		if coin_flip( t/D(t) ) then 0
	 * 		else 2 + Geom( D(t) )
	 *


	solO = Solve[w == z + w^2/(1-w), w]
	Df = w /. solO[[1]]
	Bp = (z + Df)/2 // Simplify
	phi = Exp[Bp]

	(* Mathematica can't solve psi==1. Use different approach from AoP paper *)
	(* psi = z*D[phi,z]/phi *)
	sol = Solve[3*z^4 -28*z^3 +70*z^2-58*z+8==0, z];
	y = z /. sol[[3]];
	norm = Simplify[phi /. {z->y}];

	phiD = 1/(1- z/(1-z)) // Simplify
	psiD = z*D[phiD,z]/phiD //Simplify
	solD = Solve[psiD==1,z]
	yD = z /. solD[[1]] // Simplify
	rhoD = (z / phiD) /.{z -> yD} // Simplify

	(* determine constants *)
	t = 0.171
	lambda = N[ Bp /.{z->t}, 200]
	p = N[ (z/(2*Bp)) /.{z->t}, 200]
	pp = N[ (z/Df) /.{z->t}, 200]
	nu = N[1-Df/.{z->t}, 200]
	Bpp = D[Bp, z]
	exsize = N[ z*Bpp  /.{z->t}, 200]
	
	                       3 - z
	         z (5 + --------------------)
	                Sqrt[1 + (-6 + z) z]
	exsize = ----------------------------
	                      8

	*/

	/* choose value for parameter t that optimizes target size 
	 * double t = 0.171572;
	 * --> target size about 27.3686
	//double lambda = 0.231954;
	//double p = 0.36984;
	//double pp = 0.586899;
	//double nu = 0.707663;
	*/


	/* optimization
	//INT min = 2;
	//INT fail = 0;
	//INT tfail = 0;
	//INT watch = 10;
	*/	

	INT target = 0;		// 0 corresponds to rhoD

	while(*counter > 0) {
		//adjust target size
		/*
		if(dcounter[watch]==0) {
			// we are probably using a non-optimal target size
			for(i=min;i<=max;i++) {
				if(dcounter[i]>0) {
					min = i;
					watch = i;
					target = i+tid;
					break;
				}
			}
		}
		*/



		// take sample with size smaller than max
		N = gsl_ran_poisson(rgen, oboltz[target].lambda);
		//debug
		//printf("max: %ld, poisson value: %ld\n",max, N);
		if(N>max) {
			// optimization
			// fail++;	
			continue;
		}

		size = N;
		v=0;
		for(i=0; i<N; i++) {
			if( gsl_ran_bernoulli(rgen, oboltz[target].p) == 1 ) {
				K[v] = 0;
				v++;
			} else {
				// beware that gsl_ran_geometric(nu) has density
				// p_k = nu (1-nu)^{k-1}, k>=1
				trec=0;
				while(trec >=  0) {
					if( gsl_ran_bernoulli(rgen, oboltz[target].pp) == 1 )  {
						K[v] = 0;
						v++;
						trec--;
					} else {
						K[v] = 1 + gsl_ran_geometric(rgen, oboltz[target].nu);
						size += K[v]-1;
						trec += K[v]-1;
						if(size > max) {
							break;
						}
						v++;
					}
				}
			}
		}
		// debug
		// printf("poisson=%ld, size=%ld\n", N, size);

		if(size>max || size<=1) {
			// optimization
			//fail++;	
			continue;	// we ignore structures of size 0 or 1
		}

		// check counter
		if(dcounter[size] != 0) {
			// yay, we found a valid decoration
		
			// debug	
			//printf("target %ld, config %ld, previous fails %ld\n", target, size, fail);
			// optimization
			// tfail += fail;		
			// fail = 0;

			/* begin of part that is partially locked by mutex */
			pthread_mutex_lock(mut);
			// check again if we still need to do this
			if(dcounter[size] != 0) {		

				/*debug
				printf("counter: %ld, size: %ld, dcounter[size] = %ld, degprofile[size] = %ld\n", *counter, size, dcounter[size], degprofile[size]);
				printf("components = %ld, v = %ld\n", N, v);
				printf("K = ");
				for(i=0;i<v; i++)
					printf("%ld, ", K[i]);
				printf("\n");
				*/

				// save solution 
				deco[size][degprofile[size]-dcounter[size]].numcomp = N;
				memcpy(deco[size][degprofile[size]-dcounter[size]].stor, K, v*sizeof(INT) );

				// update counters
				dcounter[size] -= 1;	
				*counter -= 1;

				//debug
				//if(dcounter[size]==degprofile[size]) printf("%ld size complete, all of %ld samples\n", size, degprofile[size]);
			}
			// unlock mutex - there's still work to do
			pthread_mutex_unlock(mut);
		} else {
			// optimization
			//fail++;
		}
	}

	// debug
	//printf("total fails: %ld\n", tfail);

	// we found enough configurations 
				
	// time to clean up...
	free(K);

	// ... and go home
	return (void *) 0;
}




/*
 * Generate decorations for outerplanar graphs
 */
struct outerdec **outerdeco(INT *degprofile, INT max, int numThreads, gsl_rng **rgens) {
	INT i,j;
	struct outerarg *argList;   	// arguments for the separate threads
	pthread_t *th;			// array of threads
	void *ret;
	INT *dcounter;

	// mutex for thread synchronization
	pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;


	// sanity check
	if(max <= 0) {
		// nothing to decorate here
		fprintf(stderr, "Error in function outerdeco. Nothing to decorate here.\n");
		exit(-1);
	}


	// allocate space for pointers to graphs
	struct outerdec **deco = malloc(sizeof(struct outerdec *)*(max+1));
	if(deco == NULL) {
		// memory allocation error
		fprintf(stderr, "Memory allocation error in function outerdeco\n");
		exit(-1);
	}
	deco[0] = NULL;	// no graphs of size zero
	deco[1] = NULL; // will be hard-coded
	for(i=2; i<=max; i++) {
		if(degprofile[i] >0) {
			deco[i] = malloc(sizeof(struct outerdec) * degprofile[i]);
			if(deco[i] == NULL) {
				// memory allocation error
				fprintf(stderr, "Memory allocation error in function outerdeco\n");
				exit(-1);
			}
			for(j=0; j<degprofile[i]; j++) {
				deco[i][j].stor = malloc(sizeof(INT)*2*i);
				if(deco[i][j].stor == NULL) {
					// memory allocation error
					fprintf(stderr, "Memory allocation error in function outerdeco\n");
					exit(-1);
				}
			}
		} else {
			deco[i] = NULL;
		}
	}




	/* optimization for small values */
	/*
	INT choice;
	for(i=2; i<=outerprec && i<=max; i++) {
		for(j=0; j<degprofile[i]; j++) {
			//debug
			//printf("%ld %ld\n",i,j);

			choice = gsl_ran_discrete(rgens[0],outergen[i]);
			memcpy(deco[i][j],outerpart[i][choice],(outerpart[i][choice][0]+1)*sizeof(INT));
		}
	}
	*/


	// the variables counter and dcounter will only be read / modified
	// in an area protected by the mutex mut
	dcounter = malloc(sizeof(INT) * (max+1));
	if(dcounter == NULL) {
		// memory allocation error
		fprintf(stderr, "Memory allocation error in function outerdeco\n");
		exit(-1);
	}
	INT counter=0;
	for(i=outerprec+1; i<=max; i++) {
		dcounter[i] = degprofile[i];
		counter += degprofile[i];
	}


	if(counter > 0) {
		INT *truncdegprofile = calloc(max+1, sizeof(INT)); // initializes to zero
		if(truncdegprofile  == NULL) {
			// memory allocation error
			fprintf(stderr, "Memory allocation error in function outerdeco\n");
			exit(-1);
		}
		for(i=outerprec+1; i<=max; i++)
			truncdegprofile[i] = degprofile[i];


		// pack list of arguments
		argList = calloc(numThreads, sizeof(struct targ));
		for(i=0; i<numThreads; i++) {
			argList[i].deco = deco;
			argList[i].max = max;
			argList[i].degprofile = truncdegprofile;
			argList[i].rgen = rgens[i];
			argList[i].counter = &counter;
			argList[i].dcounter = dcounter;
			argList[i].mut = &mut;
			argList[i].tid = i;
		}



		/* launch threads */
		th = calloc(numThreads, sizeof(pthread_t));
		for(i=0; i<numThreads; i++) {
			if(pthread_create(&th[i], NULL, &outerblock, &argList[i] )) {
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


