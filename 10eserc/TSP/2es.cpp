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
   int n_migr = 50; //ogni quanto i nodi si scambiano info sul best path
   int size, rank, scambio, scambio2, scambio3;
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
            for(int i=0; i<=rank;i++) {Primes >> p1>> p2;}
        } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   int seed[4];
   ifstream in; 
   in.open("seed.in");
   in >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   in.close(); 

  Rnd->SetRandom(seed,p1,p2);


  Path pop(Rnd,rank);


  for(int istep = 0; istep < nstep; istep++){ 
	pop.Mutazione();
	pop.NuovaPop();
	pop.BestL(rank);
	if(istep%n_migr==0){ //ogni 100 step leggo dai file di output qual è il best path e faccio scambio tra nodi
		bestpath = pop.ReadBestPath();

		scambio = int(rand()%2) + 1;
		if(scambio == 1) {scambio2=2; scambio3=3;}
		else if(scambio == 2) {scambio2=1; scambio3=3;}
		else {scambio2=1; scambio3=2;}

		if(rank == 0){
			MPI_Isend(bestpath, 32, MPI_INTEGER, scambio, itag, MPI_COMM_WORLD, &req);
			MPI_Recv(bestpath, 32, MPI_INTEGER, scambio, itag, MPI_COMM_WORLD, &stat[1]);
		        }
		else if(rank == scambio){
				MPI_Send(bestpath, 32, MPI_INTEGER, 0, itag, MPI_COMM_WORLD);
				MPI_Recv(bestpath, 32, MPI_INTEGER, 0, itag, MPI_COMM_WORLD, &stat[0]);
			}
		else if(rank == scambio2){
				MPI_Isend(bestpath, 32, MPI_INTEGER, scambio3, itag, MPI_COMM_WORLD, &req2);
				MPI_Recv(bestpath, 32, MPI_INTEGER, scambio3, itag, MPI_COMM_WORLD, &stat[3]);
			}
		else if(rank == scambio3){
				MPI_Send(bestpath, 32, MPI_INTEGER, scambio2, itag, MPI_COMM_WORLD);
				MPI_Recv(bestpath, 32, MPI_INTEGER, scambio2, itag, MPI_COMM_WORLD, &stat[2]);
			}

		pop.SetBestPath(bestpath);
	}
	
  }

  pop.BestPath(rank); 

 MPI::Finalize();

 return 0;

}

