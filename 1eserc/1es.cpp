#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "lib.h"

using namespace std; 

int main (int argc, char *argv[]){ 

   int M = 100000;
   int N = 100;
   int L = M/N; 
   double casual[M]; 
   double varianza[M];

   double intervallo = 0;
   int n_int = 100; //numero di suddivisioni dell'intervallo 0,1
   int count[n_int];
   int n_throws = 10000; //numero di throws per ogni n_int intervalli
   int interation = 100; 
   double chi2;
   double chi2_tot[interation];
   int x[interation];

   Random rnd;

   BlockingMethod blk_ave(N,L); 
   BlockingMethod blk_var(N,L); 

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

//parte 1
   for(int i=0; i<M; i++){
      casual[i] = rnd.Rannyu();
      }

   blk_ave.BlkMethod(casual, "integrale_r.dat");

//parte 2 
   for(int i=0; i<M; i++){
      varianza[i] = pow(casual[i]-0.5,2);
      }

   blk_var.BlkMethod(varianza, "varianza_r.dat");

//parte 3
   intervallo = pow(n_int, -1); //suddivido intervallo 0,1

    for (int k=0; k<interation; k++){ //ripeto test 100 volte
        
        for (int i=0;i<n_int;i++) count[i]=0;
        
        for(int i=0; i<n_throws; i++){
            	 int j=rnd.Rannyu()*n_int; //indice del subintervallo in cui finisce c
            	 count[j]++;
        }
        
	chi2 = 0;
       	for (int i=0; i<n_int; i++){
            	  chi2+=pow(count[i]-(n_throws/n_int),2)/(n_throws/n_int);
        }

        chi2_tot[k]=chi2;
	x[k] = k; 
    }

   data_chi2("chi2.dat", x, chi2_tot, interation);

 rnd.SaveSeed();

return 0; 

}
