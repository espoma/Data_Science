#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "bb_smc.h"


/* ==========================================================
Last update: 1.5.20 -- v2.0

ABC-SMC (parallel version)
Author: Marco Esposito
=========================================================== */

// =========== Define some constants =============
#define size_data 50

#define a0 0.5
#define b0 0.5
#define std_ker 0.3
#define striding_size 1
#define n 20
// =====================================



// ARGUMENT NUMBER 1: time at which the simulation starts
// ARGUMENT NUMBER 2: seed
// ARGUMENT NUMBER 3: number of particles to consider
// for the simulation (lesser for shorter runs)


int main(int argc, char *argv[]) {


	if (argc != 4) {

		printf("argc = %d\n", argc);
		printf("Usage: %s\n", argv[0]);
		printf("      %s [parameters] \n", argv[0]);
		exit(1);

	}


	/* ==================== VARIABLES DECLARATION ======================== */
	int time = atoi(argv[1]);
	int seed = atoi(argv[2]);
	int n_particles = atoi(argv[3]);

	double theta_star;
	double theta;
	double weight;
	double epsilon_t;
	int status;

	double* weights;
	weights = malloc((n_particles) * sizeof(double));

	int sample_accepted;
	/* ==================================================================== */


	/* ============== IMPORT DATA =========== */
	int y[size_data];
	status = import_data_BB(y, size_data);
	if (status == 1) {
		printf("Data not found\n");
		exit(1);
	}
	double mean_y = mean_int(y, size_data);
	// printf("(MAIN) -- mean_y: %f\n", mean_y);
	/* ======================================== */


	/* =========== READ THRESHOLD ============ */
	status = read_threshold(&epsilon_t, time);
	if (status == 1) {
		printf("Threshold not found\n");
		exit(1);
	}
	/* ======================================== */


	/* ================== set up GSL RNG =================== */
	gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set (r, seed+101); // set seed
	/* ===================================================== */

	sample_accepted = 1;

	/* =========================================================== */
	/* ================= START THE ABC ROUND ===================== */
	/* =========================================================== */
	while (sample_accepted == 1) {

		// Select particle
		status = select_particle(weights, &theta, time, n_particles, r);
		if (status == 1) {
			printf("File with particles at time %d, %d particles not present\n", time, n_particles);
			exit(1);
		}

		// Perturb the parameter
		status = ker_perturb(&theta_star, r, theta, a0, b0, std_ker);
		if (status == 1) {
			continue;
		}


		// - If the distance function is too large, start over
		// - If all the conditions are fulfilled, the sample is accepted
		sample_accepted = verify_par_BB(theta_star, mean_y, r, size_data, n, epsilon_t, time, seed);
		if (sample_accepted == 1) {
			continue;
		}


	} //END OF ABC ROUND (WHILE LOOP) -- PARTICLE WAS ACCEPTED

	calculate_weight_BB(&weight, weights, theta, theta_star, n_particles, a0, b0, std_ker);

	write_particles_BB(theta_star, time, seed, striding_size, weight);

	free(weights);

	return 0; //return of main

} //End of the program
