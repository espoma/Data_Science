module load gsl
gcc -g -Wall -o bb_smc_10 smc_bbmodel.c -I$GSL_HOME/include -L$GSL_HOME/lib -lgsl -lgslcblas -lm


