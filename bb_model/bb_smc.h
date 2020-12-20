#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#define RND gsl_rng_uniform(r)

double gsl_ran_gaussian(const gsl_rng *r, double sigma);
double gsl_ran_gaussian_pdf(double x, double sigma);
double gsl_ran_beta_pdf(double x, double a, double b);
double gsl_ran_flat(const gsl_rng *r, double a, double b);
double gsl_ran_flat_pdf(double x, double a, double b);
unsigned int gsl_ran_binomial(const gsl_rng * r, double p, unsigned int n);


double mean_int(int array[], int array_size) {
	double avg = 0;
	int h;
	for (h=0; h<array_size; h++) {
		avg += array[h];
	}
	avg /= array_size;
	return avg;
}


// double ABC_distance_BLR(double y[], double y_s[], int size_data) {
// 	double dist = 0;
// 	int k;
// 	for (k=0; k<size_data; k++) {
// 		dist += (y[k] - y_s[k]) * (y[k] - y_s[k]);
// 	}
// 	return dist;
// }



// double prior(int l, double theta_star, double a, double b) {
// 	double res;
// 	if (l==0) {
// 		res = gsl_ran_flat_pdf(theta_star, a, b);
// 		return res;
// 	}
// 	if (l==1) {
// 		res = gsl_ran_flat_pdf(theta_star, a, b);
// 		return res;
// 	}
// 	else {return -1;}
// }



double prior_BB(double theta_star, double a, double b) {

	double res;
	res = gsl_ran_beta_pdf(theta_star, a, b);
	return res;

}

// double sort(double number[], int n) {

// 	int i, j;
// 	double a;
//     for (i = 0; i < n; i++) {
//         for (j = i + 1; j < n; j++) {
//             if (number[i] > number[j]) {
//                 a =  number[i];
//                 number[i] = number[j];
//                 number[j] = a;
//             }
//         }
//     }
//     return 0;
// }


int select_particle(double *weights, double *theta, int time, int n_particles, const gsl_rng *r) {

	FILE * w_theta;
	char weights_theta[500];
	double r2;
	double aux = 0.;
	int i = 0;
	int indicator = -1;

	double* theta_vec;
	theta_vec = malloc((n_particles) * sizeof(double*));

	sprintf(weights_theta, "w_theta_time_%d_part_%d.txt", time, n_particles);
	w_theta = fopen(weights_theta, "r");
	if (w_theta == NULL) {
		printf("File %s not found!!!\n", weights_theta);
		return 1;
	}
	fscanf(w_theta, "%*[^\n]");

	r2 = RND;
	for (i=0; i<n_particles; i++) {
		fscanf(w_theta, "%lf %lf", &weights[i], &theta_vec[i]);
		aux += weights[i];
			if (r2 > aux && indicator < 0) {
				*theta = theta_vec[i];
				indicator = 1;
		}
	}

	// for (i=0; i<n_particles; i++) {
	// 	fscanf(w_theta, "%lf %lf", &weights[i], &theta_vec[i]);
	// 	aux += weights[i];
	// 	if (r2<aux) {
	// 		if (i_th<0) {
	// 			i_th = i;
	// 		}
	// 	}
	// }

	fclose(w_theta);
	free(theta_vec);

	return 0;

}


void calculate_weight_BB(double *weight, double *weights, double theta, double theta_star, int n_particles, double a0, double b0, double std_ker) {
	// Calculates the new weights for the successive ABC round
	int i;
	double tot_prior_pdf = 1.;
	double kernel_value_pdf[n_particles];
	for (i=0; i<n_particles; i++) {kernel_value_pdf[i] = 1.;}

	tot_prior_pdf *= prior_BB(theta_star, a0, b0);

	double denominator = 0.;
	for (i=0; i<n_particles; i++) {
		kernel_value_pdf[i] *= gsl_ran_gaussian_pdf(theta_star - theta, std_ker);
		denominator += weights[i] * kernel_value_pdf[i];
	}
	*weight = tot_prior_pdf/(float)(denominator);

}


int ker_perturb(double *theta_star, const gsl_rng *r, double theta, double a0, double b0, double std_ker) {
	//Perturb the i_th parameter theta via a gaussian kernel with std = std_ker
	*theta_star = gsl_ran_gaussian(r, std_ker) + theta;
	//If the prior on theta_star[l] is zero when calculated in theta_star[l],
	//the new particle is rejected and we start over
	if (prior_BB(*theta_star, a0, b0) == 0.) {
		return 1;
	}

	return 0;

}


int read_threshold(double* epsilon_t, int time) {

	FILE * epsilon;
	char epsilon_file[500];
	sprintf(epsilon_file, "epsilon_%d.txt", time);
	epsilon = fopen(epsilon_file, "r");
	if (epsilon == NULL) {
		printf("File %s not found!!!\n", epsilon_file);
		return 1;
	}
	fscanf(epsilon, "%lf", epsilon_t);
	fclose(epsilon);

	return 0;

}

// Distance function equal to zero for sample to be accepted
int verify_par_BB(double theta_star, double mean_y, const gsl_rng *r, int size_data, int n, double epsilon_t, int time, int seed) {

	int k;
	int g[size_data];
	double mean_g;
	for (k=0; k<size_data; k++) {g[k] = gsl_ran_binomial(r, theta_star, n);}
	mean_g = mean_int(g, size_data);
	printf("mean g: %f\nmean_y: %f\nepsilon_t: %f\n\n", mean_g, mean_y, epsilon_t);
	if ( fabs(mean_g - mean_y) <= epsilon_t ) {
		FILE* epsilon_file;
		char epsilon_char[500];
		sprintf(epsilon_char, "epsilon_%d_%d.txt", time+1, seed);
		epsilon_file = fopen(epsilon_char, "w");
		fprintf(epsilon_file, "%f", fabs(mean_g - mean_y));
		fclose(epsilon_file);
		return 0;
	}

	return 1;

}



void write_particles_BB(double theta_star, int time, int seed, int striding_size, double weight) {
	// Writes ABC particles on file
	FILE *output_file;
	char output_char[500];

	sprintf(output_char, "w_theta_time_%d_%d.txt", time+1, 1+(seed-1)/striding_size);
	output_file = fopen(output_char, "a");
	fprintf(output_file, "%.10f %.5f\n", weight, theta_star);
	fclose(output_file);

}


int import_data_BB(int *y, int size_data) {
	// Import data for beta-binomial model
	int d;

	FILE *file_data;
	char file_name[500];

	sprintf(file_name, "data_bb.txt");
	file_data = fopen("data_bb.txt", "r");
	if (file_data == NULL) {
		return 1;
	}
	for (d=0; d<size_data; d++) {
		fscanf(file_data, "%d", &y[d]);
	}
	fclose(file_data);

	return 0;

}




void write_epsilon(double distance, int time, int seed, int striding_size) {
	// Write epsilon for successive ABC round to file
	FILE *epsilon_;
	char eps_char[500];
	sprintf(eps_char, "epsilon_%d_%d.txt", time+1, 1+(seed-1)/striding_size);
	epsilon_ = fopen(eps_char, "a");
	fprintf(epsilon_, "%f\n", distance);
	fclose(epsilon_);

}


int read_epsilon_3d(double *eps_t, int time) {
	// Read 3d epsilon
	FILE *epsilon_;
	char eps_char[500];

	sprintf(eps_char, "epsilon_%d.txt", time);
	epsilon_ = fopen(eps_char, "r");
	if (epsilon_ == NULL) {
		return 1;
	}

	fscanf(epsilon_, "%lf %lf %lf", &eps_t[0], &eps_t[1], &eps_t[2]);
	fclose(epsilon_);

	return 0;

}



void write_epsilon_3d(double *distance, int time, int seed, int striding_size) {

	FILE *epsilon_;
	char eps_char[500];
	sprintf(eps_char, "epsilon_%d_%d.txt", time+1, 1+(seed-1)/striding_size);
	epsilon_ = fopen(eps_char, "a");
	fprintf(epsilon_, "%f %f %f\n", distance[0], distance[1], distance[2]);
	fclose(epsilon_);

}
