#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "lib.h"
#include "random.h"

using namespace std;

int main(int argc, char *argv[]){

  int which = 0;

  int nstep = 100; 
  int nblk = 60;
  int nstep_blk = 1000; //per blk
  int indice = 0;

  Random *rnd = new Random();

  ofstream out;

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
            rnd->SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

  cout << "Cities on a circumference (0) or cities inside a square (1)?";
  cin >> which;

  Path pop(rnd,which);
  double temp = 20;

	for(int i=0;i<100;i++){ //100 temperature step

  		for(int istep = 0; istep < nstep; istep++){ 
			pop.Mutazione(temp); //mutazione
			pop.outL(which,temp); //scrivere best L su output
		}
		temp = temp*0.9;
		nstep = nstep + 150;
	}
	
  
  pop.BestPath(which); 


 return 0;

}





