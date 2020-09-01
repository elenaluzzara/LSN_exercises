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
  int k = 0;
  int nstep = 100;
  double sum = 0;
  double t[N];

  double S_T[M];
  double _Z=0;
  double price[M]; 
  int S0 = 100; 
  int T = 1; 
  int K = 100; 
  double r = 0.1; 
  double sigma = 0.25; 
 

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

//diretto C
   BlockingMethod blk_dirC(N,L);

   for(int i = 0; i < M; i++){
     _Z = rnd.Gauss(0,1); 
     S_T[i] = S0*exp((r-pow(sigma,2)/2)*T + sigma*_Z*sqrt(T)); 
   }

   for(int i=0; i<M;i++){
     price[i] =exp(-r*T)*max(0,S_T[i] - K);
   }

   blk_dirC.BlkMethod(price);
   blk_dirC.Stampa("direct_C.dat");

//diretto P
   BlockingMethod blk_dirP(N,L);

   for(int i=0; i<M;i++){
     price[i] =exp(-r*T)*max(0, K - S_T[i]);
   } 

   blk_dirP.BlkMethod(price);
   blk_dirP.Stampa("direct_P.dat");

//C discreto
   
   BlockingMethod blk_discrC(N,L);

   for(int i = 1; i < nstep; i++){
     t[0] = T*pow(nstep,-1);
     t[i] = t[i-1] + T*pow(100,-1);
   }

   for(int i = 0; i < M; i++){
    double S[nstep];
    for(int j= 0; j < nstep; j++){
    
     if(j == 0){
       S[j] = S0;
     }

     else{   
       _Z = rnd.Gauss(0,1);
       S[j] = S[j-1]*exp((r-pow(sigma,2)/2)*(t[j]-t[j-1]) + sigma*_Z*sqrt(t[j]-t[j-1])); 
     }

    }

    S_T[i] = S[nstep-1]; //final asset
   }

   for(int i=0; i<M;i++){
     price[i] =exp(-r*T)*max(0,S_T[i] - K);
   }

   blk_discrC.BlkMethod(price);
   blk_discrC.Stampa("discrete_C.dat");


//P discreto 
   BlockingMethod blk_discrP(N,L);

   for(int i=0; i<M;i++){
     price[i] =exp(-r*T)*max(0, K - S_T[i]);
   }

   blk_discrP.BlkMethod(price);
   blk_discrP.Stampa("discrete_P.dat");


 rnd.SaveSeed();

return 0; 

}
