#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "lib.h"

using namespace std; 

int main (int argc, char *argv[]){ 

   int M = 1000000;
   int N = 100;
   int L = M/N; 

   int tipo_matr=0;//per scegliere tra distr uniforme e gaussiana
   int state=0; //per scegliere se campionare stato fond o 1 stato eccitato
   int equil =0; //per scegliere se equilibrare oppure se è stato già fatto

   int M_equ = 1000; //step per equilibrare
   double x_equ[M_equ]; //coordinate che prendo da config.final 
   double y_equ[M_equ];
   double z_equ[M_equ];

   double acc = 0;
   double alpha = 0;
   int n = 0; 

   Random rnd;

   double * x = new double[M];
   double * y = new double[M];
   double * z = new double[M];
   double * r = new double[M];
     
   double delta=0;
   double sigma=0;

   double x_0 =0;
   double y_0 =0;
   double z_0 =0;

   ifstream ReadInput;
   ifstream in; //per ripartire da configurazione di simulazione precedente
   ofstream out; //per salvare coordinate e vedere quando si è raggiunto equilibrazione
   ofstream WriteDensity; //per scrivere x y e z dopo equilibrazione

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

   ReadInput.open("input.dat");
   ReadInput >> state;
   if(state==0) cout << "The state is the ground state." << endl;
   else cout << "The state is the excited state 2p" << endl;
   ReadInput >> tipo_matr; 
   if(tipo_matr==0) cout << "The transfer matrix is uniform" << endl;
   else cout << "The transfer matrix is gaussian" << endl;
   ReadInput >> delta; 
   if(tipo_matr==0) cout << "The amplitude of the uniform intervall is: " << delta << endl;
   ReadInput >> sigma;
   if(tipo_matr!=0) cout << "The amplitude of the Gaussian is: " << sigma << endl;
   ReadInput >> x_0;
   ReadInput >> y_0;
   ReadInput >> z_0;
   cout << "The initial points are: x=" << x_0 << " y=" << y_0 << " z=" << z_0 << endl;
   ReadInput >> equil; 
   if(equil!=0) cout << "Equilibration.." << endl;
   ReadInput.close();

   BlockingMethod blk(N,L);

   if(equil == 0){ //scrivo le coord x y z trovate con metropolis
     if(tipo_matr == 0) WriteDensity.open("risultati/ProbDensity_unif_" + to_string(state) + ".xyz");
     else WriteDensity.open("risultati/ProbDensity_gauss_" + to_string(state) + ".xyz");
   }

   if(equil!=0) M = M_equ;

   for(int i = 0; i < M; i++){

	if(i==0){
		 if(equil != 0){ //punti iniziali
     			 x[i] = x_0;
     			 y[i] = y_0;
     			 z[i] = z_0;
     		 }

		 else{
  		  cout << "Reading the first coordinates from the last in config.final" << endl << endl;
  	 	  in.open("config.final");
		  for(int i=0; i<M_equ;i++){
  		  	in >> x_equ[i] >> y_equ[i] >> z_equ[i];
		  }
		  in.close();

     		  x[i] = x_equ[M_equ -1];
     		  y[i] = y_equ[M_equ -1];
     		  z[i] = z_equ[M_equ -1];
		 }
	}

     	else{

		if(tipo_matr==0){
      	 	  x[i] = x[i-1] + rnd.Rannyu(-delta, delta); //matrice di trasferimento T
      	 	  y[i] = y[i-1] + rnd.Rannyu(-delta, delta);
      	 	  z[i] = z[i-1] + rnd.Rannyu(-delta, delta);
		}
     		else{
     	 	  x[i] = x[i-1] + rnd.Gauss(0, sigma); //matrice di trasferimento T
     		  y[i] = y[i-1] + rnd.Gauss(0, sigma);
      		  z[i] = z[i-1] + rnd.Gauss(0, sigma);
     		}

     		if(state == 0) acc = p_gs(x[i], y[i], z[i])/p_gs(x[i-1], y[i-1], z[i-1]);
     		else acc = p_1e(x[i], y[i], z[i])/p_1e(x[i-1], y[i-1], z[i-1]);

     		if(acc < 1) alpha = acc;
     		else alpha = 1;

		if(rnd.Rannyu() > alpha){ 
	  	 x[i] = x[i-1]; 
	   	 y[i] = y[i-1]; 
	  	 z[i] = z[i-1]; 
	  	 n++;
        	}
	}

	if(equil == 0) WriteDensity << x[i] << "   " << y[i] << "   " << z[i] << endl;
   }
   if(equil == 0) WriteDensity.close();

  double accept = 100*(1 - n*pow(M,-1));
  cout << "The percentage of the accepted configurations is: " << accept << "%" << endl;   

  double sum =0;
  for(int i=0; i<M;i++){
     r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
     if(i>100) sum+= r[i]; //rimuovo i primi 100 per una media meno influenzata dai primi passi
  }

  if(equil!=0) cout << "Equilibration: the mean value of r after " << M_equ << " steps is: " << sum/(M_equ-100) << endl; 

  if(equil!=0){
     out.open("risultati/equilibration_r" + to_string(state) + ".dat");
     for(int i=0; i<M_equ;i++){
  	out << i << "  " << r[i] << endl; 
     }
     out.close();
     out.open("config.final");
     for(int i=0; i<M_equ;i++){
  	out << x[i] << "  " << y[i] << "  " << z[i] << endl; 
     }
     out.close();
  }

  if(equil == 0){
    blk.BlkMethod(r);
    if(tipo_matr == 0) blk.Stampa("risultati/r_unif_" + to_string(state) + ".dat");
    else blk.Stampa("risultati/r_gauss_" + to_string(state) + ".dat");
  }


 rnd.SaveSeed();

 delete[] x;
 delete[] y;
 delete[] z;
 delete[] r;



return 0; 

}


