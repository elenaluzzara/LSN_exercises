#include <iostream> 
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include "random.h"

using namespace std;

double * metropolis(int M_step, double x0, double delta, double mu0, double sigma0);

void metropolis_param(int M, double x0, double delta, double delta_mu, double delta_sigma, double mu0, double sigma0, double temp);

double H_ave(int M, double x[], double mu, double sigma);

double psi_trial(double x, double mu, double sigma);

double prob_integr(double x, double mu, double sigma);

double H_integranda (double x, double mu, double sigma);

void data_blocking(int M, int N, double x[], double mu, double sigma);

double error(double AV[],double AV2[],int i);

