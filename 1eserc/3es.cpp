#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "lib.h"

using namespace std; 

int main (int argc, char *argv[]){ 

   int M = 200000;
   int N = 100;
   int L = M/N;
   double l = 0.8;
   double d = 1.;
   double punto[M];
   double lato[M];
   double x_circ[M];
   int count = 0;
   int k = 0;

   double pigreco[N]; //N];

   double x=0;
   double y=0; 
//double cerchio_pigreco;

   Random rnd;

   BlockingMethod blk(N,L); 


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

//calcolo pigreco con ago di buffon
   for(int i=0; i<M; i++){ //punto in cui cade una punta dell'ago
      punto[i] = rnd.Rannyu(0, d);
   }

   for(int i=0; i<M;){
      double angle = 0;
      x = rnd.Rannyu(-l, l);
      y = rnd.Rannyu(-l, l);
	if((pow(x,2) + pow(y,2)) < l*l){ //metodo accept-reject per calcolare l'inclinazione dell'ago

	  x_circ[i] = x*l*(1./sqrt(pow(x,2) + pow(y,2))); //coordinata x sulla circonferenza ricavata proiettando lungo il raggio il punto (x,y) estratto
	  i++;
	} 
   }

  for(int i=0; i<N; i++){
	count = 0;
	for (int j=0; j<L; j++){ 

		k =j+i*L;
		if(x_circ[k] + punto[k] <= 0 || x_circ[k] + punto[k] >= d){ //se questo vale l'ago ha intersecato l'asse x = 0 o x = d 
		  count++;
		}
	}
        pigreco[i] = (2*L*l)/(count*d); 
  }
  
  blk.SetAve(pigreco); //pigreco è un array con la media nei vari blocchi
  blk.BlkMethod_pigreco("pigreco_buffon.dat"); //qui calcolo l'errore per blocco

 rnd.SaveSeed();

return 0; 

}
