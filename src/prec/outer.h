/*
 * precomputation of probability weights for outerplanar graphs
 */

INT outermax;
double *outerq;
INT outerprec;
double **outerp;
INT *outernum;
INT ***outerpart;
gsl_ran_discrete_t **outergen;

struct outerboltz{
	double t;
	double lambda;
	double p;
	double pp;
	double nu;
};

struct outerboltz *oboltz;

void outerinit() {


/*
 * Weights for balls in boxes model
 * The first 10k weights were obtained using Mathematica's Series command
 * Computation took about 1 day for this first phase
 * The remaining 90k weights were obtained using a singularity expansion 
 * up to order 50
 * Relative error: (original / approximation) = 1 + O(n^(-50))
 * Calculating xi[10000] / approximation yields
 * 1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000037880732521372734346576927930278`200.
 *
*/



outermax = 100001;
outerq = (double *) malloc(outermax*sizeof(double));
if(!outerq) {
	fprintf(stderr, "Memory allocation error in function outerinit\n");
	exit(-1);
}

#include "outer_q.h"



/*
 * Preprocessing works great for cacti graphs, but effectiveness appears to be
 * limited for outerplanar graphs. We keep the structures here for future
 * optimizations.
 */


// we preprocess uniform generators for decorations from 2 until outerprec
outerprec = 1;

// array that holds the number of decorations of a given size
outernum = (INT *) malloc((outerprec+1)*sizeof(INT));
if(!outernum) {
	fprintf(stderr, "Memory allocation error in function outerinit\n");
	exit(-1);
}

memcpy(outernum, (INT [])
{0,0,0,0,0,0,0},
(outerprec+1)*sizeof(INT));





// array that holds the probabilities for each of the decorations
outerp = (double **) malloc((outerprec+1)*sizeof(double *));
if(!outerp) {
	fprintf(stderr, "Memory allocation error in function outerinit\n");
	exit(-1);
}

for(INT i=0; i<=outerprec; i++) {
	outerp[i] = (double *) malloc( outernum[i] * sizeof(double));
	if(!outerp[i]) {
		fprintf(stderr, "Memory allocation error in function outerinit\n");
		exit(-1);
	}
}


//memcpy(outerp[0], (double []){1.0}, outernum[0]*sizeof(double));
//memcpy(outerp[1], (double []){1.0}, outernum[1]*sizeof(double));
//memcpy(outerp[2], (double []){0.5,0.5}, outernum[2]*sizeof(double));



// array that holds the actual decorations
outerpart = (INT ***) malloc((outerprec+1)*sizeof(INT **));
if(!outerpart) {
	fprintf(stderr, "Memory allocation error in function outerinit\n");
	exit(-1);
}

for(INT i=0; i<=outerprec; i++) {
	outerpart[i] = (INT **) malloc( outernum[i] * sizeof(INT *));
	if(!outerpart[i]) {
		fprintf(stderr, "Memory allocation error in function outerinit\n");
		exit(-1);
	}
}



// the random number generators

outergen = (gsl_ran_discrete_t **) malloc((outerprec+1)*sizeof(gsl_ran_discrete_t *));
if(!outergen) {
	fprintf(stderr, "Memory allocation error in function outerinit\n");
	exit(-1);
}

// set them to NULL for on demand generation later on
for(INT i=0; i<=outerprec; i++) {
	outergen[i] = NULL;
}

for(INT i=2; i<=outerprec; i++) {
	outergen[i] = gsl_ran_discrete_preproc(outernum[i], outerp[i]);
}


oboltz = malloc(100001*sizeof(struct outerboltz));
if(!oboltz) {
	fprintf(stderr, "Memory allocation error in function outerinit\n");
	exit(-1);
}

#include "outer_boltz_part1.h"
#include "outer_boltz_part2.h"

}
