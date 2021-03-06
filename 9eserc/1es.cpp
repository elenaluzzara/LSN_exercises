#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "lib.h"
#include "random.h"

using namespace std;

int main(int argc, char *argv[]){

// così dovrei creare 100 array (dovrebbero bastare 100, ma si possono aumentare se non ci sono risultati soddisfacenti, dovrebbe rimanere costante il numero della popolazione poi magari ridurlo quando si sono trovati ) ognuno con 32 elementi con un certo ordine
  int which = 0;

  int nstep = 2000;

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

  pop.InitialPath(which);

  		for(int istep = 0; istep < nstep; istep++){ 

			pop.Mutazione(); //operatori di mutazione
			pop.NuovaPop(); //operatori di crossover e creazione di una nuova popolazione
			pop.BestL(which); //stampa output di best L
			pop.AverageL(which); //stampa output della media delle prime 50 best L
		}	

  pop.BestPath(which); 


 return 0;

}





