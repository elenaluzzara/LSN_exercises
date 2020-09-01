/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // methods
  void SetRandom(int * , int, int); //per inizializzarlo 
  void SaveSeed(); 
  double Rannyu(void); //generazione uniforme tra 0 e 1
  double Rannyu(double min, double max); //generare numero uniforme tra min e max
  double Gauss(double mean, double sigma); //box muller
  double Exp(double lamba);
  double Lorentz(double mu, double gamma);
};

#endif // __Random__

double * metropolis(int M, double x0, double delta, double mu, double sigma);


double H_ave(int M, double x[], double mu, double sigma);


double psi_trial(double x, double mu, double sigma);

double prob_integr(double x, double mu, double sigma);

double H_integranda (double x, double mu, double sigma);

void data_blocking(int M, int N, double x[], double mu, double sigma);

double error(double AV[],double AV2[],int i);


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
