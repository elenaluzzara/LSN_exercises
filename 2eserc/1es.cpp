#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "lib.h"

using namespace std; 

int main (int argc, char *argv[]){ 

   int M = 160000;
   int N = 100;
   int L = M/N; 

   double casual[M]; 
   double integranda[M];

   double sampling[M];
   double integranda_impsamp[M];

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

   BlockingMethod blk(N,L); //per uniform distribution
   BlockingMethod blk_impsamp(N,L); //per l'importance sampling

   for(int i=0; i<M; i++){
      casual[i] = rnd.Rannyu();
   }
   
//uniform distribution

   for(int i = 0; i < M; i++){
	integranda[i] = M_PI*0.5*cos(M_PI*casual[i]*0.5);
   }

   blk.BlkMethod(integranda);
   blk.Stampa("integrale_unifdistr.dat");

//importance sampling


   for(int i=0; i<M;i++){
      sampling[i] = 1 - sqrt(1-casual[i]);   
   }

  for(int i = 0; i < M; i++){
	integranda_impsamp[i] = M_PI*0.5*cos(sampling[i]*M_PI*0.5)/(2*(1-sampling[i]));
  } 

   blk_impsamp.BlkMethod(integranda_impsamp);
   blk_impsamp.Stampa("integrale_impsamp.dat");


 rnd.SaveSeed();

return 0; 

}
