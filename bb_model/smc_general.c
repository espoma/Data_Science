#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "bb_smc.h"


/* ==========================================================
Last update: 26.4.20 -- v1.1

ABC-SMC (parallel version) 
Author: Marco Esposito
=========================================================== */

#define n_param 1
#define size_data 10

#define a0 0.5
#define b0 0.5
#define std_ker 0.4
#define striding_size 2
#define n 10



// ARGUMENT NUMBER 1: time at which the simulation starts
// ARGUMENT NUMBER 2: seed
// ARGUMENT NUMBER 3: number of particles to consider for the simulation (for shorter runs)


int main(int argc, char *argv[]) {

	/* ==================== VARIABLES DECLARATION ======================== */
	FILE *w_theta;
	char weights_theta[500];

	int time = atoi(argv[1]);
	int seed = atoi(argv[2]); 
	int n_particles = atoi(argv[3]);

	double weight;
	int status;
	double distance;

	int k;
	double** theta;
	theta = malloc((n_param) * sizeof(double*));
	for (k=0; k<n_param; k++) {
		theta[k] = calloc(n_particles, sizeof(double));
	}

	double* weights;
	weights = malloc((n_particles) * sizeof(double));

	double* theta_star;
	theta_star = malloc((n_param) * sizeof(double));

	double* std_kernel;
	std_kernel = malloc((n_param) * sizeof(double));

	std_kernel[0] = std_ker;

	double r2, aux;
	int i, l, i_th, l_temp;
	int sample_accepted;
	/* =============================================== */

	/* If the number of arguments is 4, then argv[3] is the number of particles to consider.
	This is useful when a shorter run is needed (necessary to specify as weights have to 
	be normalized accordingly!) */ 
	if (argc==4) {sprintf(weights_theta, "w_theta_time_%d_part_%d.txt", time, n_particles);}
	else {
		printf("argc = %d\n", argc);
		printf("Usage: %s\n", argv[0]);
		printf("      %s [parameters] \n", argv[0]);
		exit(1);
	}

	// int n_param = 1;


	/* ============== IMPORT THE FILE WITH THE DATA =========== */
	int y[size_data];
	status = import_data_BB(y, size_data);
	if (status == 1) {
		exit(1);
	}
	double mean_y = mean_int(y, size_data);
	/* ========================================================= */


	/* ==================== READ THE DISTANCE FROM A FILE ===================== */
	double eps_t;
	status = read_epsilon(eps_t, time);
	if (status == 1) {
		exit(1);
	}
	/* ========================================================================= */


	/* ================== set up GSL RNG =================== */
	gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set (r, seed+11); // set seed
	/* ===================================================== */

	sample_accepted = 1;

	/* =============================================================================== */
	/* =============================== START THE ABC ROUND =========================== */
	/* =============================================================================== */
	while (sample_accepted==0) {
	
		i_th = -1;

		aux = 0.;
		w_theta = fopen(weights_theta, "r");
		if (w_theta == NULL) {printf("File %s not found!!!\n", weights_theta);}
		fscanf(w_theta, "%*[^\n]");
		r2 = RND;
		for (i=0; i<n_particles; i++) {
			fscanf(w_theta, "%lf %lf", &weights[i], &theta[0][i]);
			aux += weights[i];
			if (r2<aux) {
				if (i_th<0) {
					i_th = i;
				}
			}
		}
		fclose(w_theta);


		double* a;
		double* b;

		a = malloc((n_param) * sizeof(double));
		b = malloc((n_param) * sizeof(double));

		a[0] = a0; b[0] = b0;	
		
		
		/* YOU CAN SIMPLY PERTURB THE PARTICLE WITH A KERNEL AT WILL WITHOUT DOING 
		THE LOOP OVER PARAMETERS AS IT'S JUST ONE IN THIS CASE. POSSIBLY IMPLEMENT 
		THIS IN A FUNCTION AND PUT IT IN THE .h FILE */
		l_temp = 0;
		//Perturb the i_th particles via a gaussian kernel with std = std_kernel
		for (l=0; l<n_param; l++) {
			theta_star[l] = gsl_ran_gaussian(r, std_kernel[l]) + theta[l][i_th];
			//If the prior on theta_star[l] is zero when calculated in theta_star[l], 
			//the new particle is rejected and we start over 
			if (prior_BB(l, theta_star[l], a0, b0) == 0.) {
				break;
			}
			l_temp = l;
		}
		if (l_temp < n_param-1) {continue;}


		/* HERE YOU CAN JUST IMPLEMENT THIS WITH A BINARY FUNCTION sample_accepted = verify_sample(distance, ..) */
		// -- If the distance function is too large, start over
		// -- If all the conditions are fulfilled, the sample is accepted
		int g[size_data];
		for (d=0; d<size_data; d++) {g[d] = gsl_ran_binomial(r, theta_star[0], n);}
		
		if (mean_int(g, size_data) - mean_y == 0) {
			sample_accepted = 1; 
		}
		else {continue;} // Start over if distance function not equal to zero




		/* ================ CALCULATE THE WEIGHTS ===================== */
		double tot_prior_pdf = 1.;
		double kernel_value_pdf[n_particles];
		for (i=0; i<n_particles; i++) {kernel_value_pdf[i] = 1.;}

		for (l=0; l<n_param; l++) {
			tot_prior_pdf *= prior_BB(l, theta_star[l], a0, b0);
		}

		double denominator = 0.;
		for (i=0; i<n_particles; i++) {
			for (l=0; l<n_param; l++) {
				kernel_value_pdf[i] *= gsl_ran_gaussian_pdf(theta_star[l] - theta[l][i], std_kernel[l]);
			}
			denominator += weights[i] * kernel_value_pdf[i];
		}	
		/* =============================================================== */

		weight = tot_prior_pdf/(float)(denominator);

		/* ============== CREATE THE FILE WITH WEIGHTS AND PARAMETERS ============= */
		write_particles_BB(theta_star, time, seed, striding_size, weight);
		/* ===================================================================== */


		free(a);
		free(b);

	} //END OF THE WHILE LOOP -- PARTICLE WAS ACCEPTED




	/* ====================== CREATE THE FILE WITH THE DISTANCE ============== */
	write_epsilon(distance, time, seed, striding_size);
	/* ======================================================================= */



	for (k=0; k<n_param; k++) {
		free(theta[k]);
	}
	free(theta);
	free(weights);
	free(theta_star);
	free(std_kernel);


	return 0;
} //End of the program
