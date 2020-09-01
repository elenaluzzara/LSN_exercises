#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "lib.h"

using namespace std; 

int main (int argc, char *argv[]){ 

int M = 10000;
int N[4]={1,2,10,100};
double SN_stard[4][M];
double SN_exp[4][M];
double SN_lor[4][M];
double sum[4];

double casual[M]; 
int x[M];

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

//stardard dice 

 for(int k = 0; k < 4; k++){ //per le 4 diverse somme

  for(int i = 0; i < M; i++){

    sum[k] = 0;

	for(int j=0; j<N[k]; j++){
      	  casual[j] = rnd.Rannyu();
      	  sum[k]+=casual[j];
        }

    SN_stard[k][i]=sum[k]/N[k];
    x[i] = i;
  }

 }

 data_histo("standard_dice.txt", x, SN_stard[0], SN_stard[1], SN_stard[2], SN_stard[3], M);

//exponential dice 

 for(int k = 0; k < 4; k++){

  for(int i = 0; i < M; i++){

    sum[k] = 0;

	for(int j=0; j<N[k]; j++){
           casual[j] = rnd.Exp(1);
      	   sum[k]+=casual[j];
        }

    SN_exp[k][i]=sum[k]/N[k];
    x[i] = i;
  }

 }

 data_histo("exp_dice.txt", x, SN_exp[0], SN_exp[1], SN_exp[2], SN_exp[3], M);

//lorentzian dice 

 for(int k = 0; k < 4; k++){

  for(int i = 0; i < M; i++){

    sum[k] = 0;

	for(int j=0; j<N[k]; j++){
      	  casual[j] = rnd.Lorentz(0,1);
      	  sum[k]+=casual[j];
        }

    SN_lor[k][i]=sum[k]/N[k];
    x[i] = i;
  }

 }

 data_histo("lorentz_dice.txt", x, SN_lor[0], SN_lor[1], SN_lor[2], SN_lor[3], M);

 rnd.SaveSeed();

return 0; 

}
