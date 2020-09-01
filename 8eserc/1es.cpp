#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include "random.h"

using namespace std; 

int main (int argc, char *argv[]){ 

   int M = 5000000;
   int N = 100;

   double delta = 0;
   double sigma0 = 0;
   double mu0 = 0;
   double x0 = 0;

   double * x = new double[M];
   
   ifstream ReadInput;

   ReadInput.open("input.dat");
   ReadInput >> mu0;
   cout << "The parameter mu is: " << mu0 << endl;
   ReadInput>> sigma0;
   cout << "The parameter sigma is: " << sigma0 << endl;
   ReadInput >> delta;
   cout << "The amplitude of the intervall needed for the metropolis algorithm is: " << delta << endl;
   ReadInput.close();


   x = metropolis(M, x0, delta, mu0, sigma0);

   double sum = 0;
   for(int i = 0; i  < M; i ++){
   sum +=  H_integranda(x[i], mu0, sigma0);
   }

   int risposta;
   cout << "Do you want to print this value in a file (1 for yes, 0 for no)? "; 
   cin >> risposta; 
   ofstream minimo; 
   if(risposta == 1){
   minimo.open("griglia_SigmaMu.out", ios::app); 
   minimo << mu0 << setw(12) << sigma0 << setw(12) << sum/M << endl; 
   minimo.close();
   }


   data_blocking(M, N, x, mu0, sigma0); 

 delete[] x;

return 0; 

}

/////+

double * metropolis(int M, double x0, double delta, double mu, double sigma){

  double acc = 0;
  double alpha = 0;
  double r = 0;
  int n = 0;

  Random rnd;

   int seed[4]; //questa prima parte va fatto sempre per implementare i semi di completamento
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in"); 
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


   double *x = new double[M];

   ofstream WriteDensity;
   WriteDensity.open("ProbDensity.xyz");

   for(int i = 0; i < M; i++){

     if(i == 0){
     x[i] = x0;    
     }

     else{

     x[i] = x[i-1] + rnd.Rannyu(-delta, delta); //matrice di trasferimento T

     acc = abs(pow(psi_trial(x[i],mu,sigma),2))/abs(pow(psi_trial(x[i-1],mu,sigma),2));
     if(acc < 1) alpha = acc;
     else alpha = 1;
     r = rnd.Rannyu();

	if(r > alpha){ 
	x[i] = x[i-1]; 
	n++;
        }
    
     }
        WriteDensity << x[i] << endl;
   }
   WriteDensity.close();

  double accept = 100*(1 - n*pow(M,-1));
  cout << "The percentage of the accepted configuration is : " << accept << "%" << endl; 

 rnd.SaveSeed();
return x;
}


double H_ave(int M, double x[], double mu, double sigma){

   double sum = 0;
   for(int i = 0; i  < M; i ++){
   sum +=  H_integranda(x[i], mu, sigma);
   }
  
 return sum/M;
} 


double psi_trial(double x, double mu, double sigma){
     double psi = exp(-pow(x-mu,2)/(2*pow(sigma,2))) + exp(-pow(x+mu,2)/(2*pow(sigma,2)));
 return psi;
}

double prob_integr(double x, double mu, double sigma){
    double psi2 = pow(psi_trial(x, mu, sigma),2); 
 return psi2;
}

double H_integranda (double x, double mu, double sigma){

  double exp_meno = exp(-pow(x-mu,2)/(2*pow(sigma,2)));
  double exp_piu = exp(-pow(x+mu,2)/(2*pow(sigma,2)));
  double potenziale = pow(x,4) - 2.5*pow(x,2); 
  double deriv2_psi = pow(x+mu,2)*exp_piu/pow(sigma,4) + pow(x-mu,2)*exp_meno/pow(sigma,4) - exp_meno/pow(sigma,2) - exp_piu/pow(sigma,2); 
 
  double H_psi = -0.5*deriv2_psi + potenziale*psi_trial(x,mu,sigma); 

 return H_psi/psi_trial(x,mu,sigma); 
}

void data_blocking(int M, int N, double x[], double mu, double sigma){

  int L = M/N; 
  double sum = 0;
  int k = 0;

  double ave[N] = {0};
  double av2[N]= {0};
  double sum_prog[N]= {0};
  double su2_prog[N] = {0}; 
  double error_prog[N]= {0};
  int asse_x[N] = {0};

   for(int i=0; i<N; i++){
	sum = 0; 
	for (int j=0; j<L; j++){
		k =j+i*L;
		sum+= H_integranda(x[k], mu, sigma);
	}
	ave[i] = sum/L; 
	av2[i] = pow(ave[i],2);
  }
  
  for(int i =0; i < N; i++){
	for(int j =0; j < i+1; j++){
		sum_prog[i] += ave[j];
		su2_prog[i] += av2[j];
	}
	sum_prog[i]/=(i+1);
	su2_prog[i]/=(i+1);
	error_prog[i] = error(sum_prog,su2_prog,i);
        asse_x[i] = i*L;
  }

   ofstream out;
   out.open("H_ave.out");
   if (out.is_open()){
	for (int i = 0; i < N; i++){
        out << asse_x[i] << " " << sum_prog[i] << " " << error_prog[i] << endl; 
   	} 
   } else cerr << "PROBLEM: Unable to open file.out" << endl;

  out.close();

return;
}

double error(double AV[],double AV2[],int i) 
{
  double varianza = 0;
  if (i == 0) return 0; 
  else {
  varianza = sqrt((AV2[i] - pow(AV[i],2))/i);
  return varianza;
  }

}











