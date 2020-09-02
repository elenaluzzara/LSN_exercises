#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <vector>
#include <climits>
#include "TSP_lib.h"
#include "mpi.h" 
 
using namespace std;

int main(int argc, char *argv[]){

   int nstep = 2000; 
   int n_migr= 80; //ogni quanto i nodi si scambiano info sul best path
   int size, rank;
   int exchange, exchange2, exchange3;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Status stat[4];
   MPI_Request req, req2;

   int *bestpath= new int[32];
   
   int itag=1;

   Random *Rnd=new Random();

   int p1, p2;
   ifstream Primes("Primes");
        if (Primes.is_open()){
            for(int i=0; i<=rank;i++) {Primes >> p1>> p2;} //ogni rank prende una riga di Primes divers
        } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   int seed[4];
   ifstream in; 
   in.open("seed.in");
   in >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   in.close(); 

   Rnd->SetRandom(seed,p1,p2); //inizializzo Rnd in ogni rank in modo diverso


   Path pop(Rnd,rank);


   for(int istep = 0; istep < nstep; istep++){ 
	pop.Mutazione();
	pop.NuovaPop();
	pop.BestL(rank);
	if(istep%n_migr==0){ //ogni 50 step cerco il bestpath dei nodi e faccio scambio
 
		bestpath = pop.ReadBestPath(); //scrivo il best path 

		exchange = int(rand()%3)+1;
		if(exchange == 1) {exchange2=2; exchange3=3;}
		else if(exchange == 2) {exchange2=1; exchange3=3;}
		else {exchange2=1; exchange3=2;}

		if(rank == 0){ //faccio scambio tra rank 0 e rank==exchange
			MPI_Isend(bestpath, 32, MPI_INTEGER, exchange, itag, MPI_COMM_WORLD, &req);
			MPI_Recv(bestpath, 32, MPI_INTEGER, exchange, itag, MPI_COMM_WORLD, &stat[1]);
		        }
		else if(rank == exchange){
				MPI_Send(bestpath, 32, MPI_INTEGER, 0, itag, MPI_COMM_WORLD);
				MPI_Recv(bestpath, 32, MPI_INTEGER, 0, itag, MPI_COMM_WORLD, &stat[0]);
			}
		else if(rank == exchange2){ //faccio scambio tra rank==exchange2 e rank==exchange3
				MPI_Isend(bestpath, 32, MPI_INTEGER, exchange3, itag, MPI_COMM_WORLD, &req2);
				MPI_Recv(bestpath, 32, MPI_INTEGER, exchange3, itag, MPI_COMM_WORLD, &stat[3]);
			}
		else if(rank == exchange3){
				MPI_Send(bestpath, 32, MPI_INTEGER, exchange2, itag, MPI_COMM_WORLD);
				MPI_Recv(bestpath, 32, MPI_INTEGER, exchange2, itag, MPI_COMM_WORLD, &stat[2]);
			}

		pop.SetBestPath(bestpath); //sovrascrivo il bestpath scambiato
	}
	
  }

  pop.BestPath(rank); 

 MPI::Finalize();

 return 0;

}

