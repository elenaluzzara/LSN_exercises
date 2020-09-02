#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include "lib.h"

using namespace std; 

int main (int argc, char *argv[]){ 

int nstep = 101; //step fatti dal random walk
int M = 10000; //esprimenti totali

int a = 1; //larghezza dello step

int *x[nstep];
int *y[nstep];
int *z[nstep];
double rN[nstep];
double rN2[nstep];
double var[nstep];
double error[nstep];

double *theta[nstep];
double *fi[nstep];
double *x_cont[nstep];
double *y_cont[nstep];
double *z_cont[nstep];
double rN_cont[nstep];
double rN2_cont[nstep];
double var_cont[nstep];
double error_cont[nstep];

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

//nstep numero di step e N numero di random walk con nstep 


  for(int i=0; i<nstep; i++){
	x[i] = new int[M];
	y[i] = new int[M];
	z[i] = new int[M];
	theta[i] = new double[M];
	fi[i] = new double[M];
	x_cont[i] = new double[M];
	y_cont[i] = new double[M];
	z_cont[i] = new double[M];
  }
  
  for(int k=0; k < M; k++){
   
    for(int i=0; i < nstep; i++){
      int scelta = 0;
      int sceltax = 0; 
      int sceltay = 0; 
      int sceltaz = 0; 

      if(i == 0){ //inizio
        x[i][k] = 0;
        y[i][k] = 0;
        z[i][k] = 0;
      }

      else{

        scelta = int(rnd.Rannyu(0,3)); //scelta tra 0 1 2 per x y z

     	if(scelta == 0){
          sceltax = int(rnd.Rannyu(0,2)); //scelta tra 0 1 per -1 o +1
		if(sceltax == 0){
	      	  x[i][k] = x[i-1][k] -a;
		}
		else{
	      	  x[i][k] = x[i-1][k] +a;	
		}	
	}
	else{
	  x[i][k] = x[i-1][k];
	}

     	if(scelta == 1){
          sceltay = int(rnd.Rannyu(0,2)); //scelta tra 0 1 per -1 o +1
		if(sceltay == 0){
	      	  y[i][k] = y[i-1][k] -a;
		}
		else{
	      	  y[i][k] = y[i-1][k] +a;	
		}
	}
        else{
	  y[i][k] = y[i-1][k];
	}

     	if(scelta == 2){
          sceltaz = int(rnd.Rannyu(0,2)); //scelta tra 0 1 per -1 o +1
		if(sceltaz == 0){
	      	  z[i][k] = z[i-1][k] -a;
		}
		else{
	      	  z[i][k] = z[i-1][k] +a;	
		}
	}
	else{
	  z[i][k] = z[i-1][k];
	}
     }
    }

  } 

 double sum;
 double sum2;

 for(int i=0; i < nstep; i++){ //calcolo rN e errore all'ultimo blocco
    sum = 0;
    sum2 = 0;

   for(int k= 0; k < M; k++){
     sum+=pow(x[i][k],2) + pow(y[i][k],2) + pow(z[i][k],2);
     sum2+=pow(pow(x[i][k],2) + pow(y[i][k],2) + pow(z[i][k],2),2);
   }
   
    rN[i] = sqrt(sum/M);
    rN2[i]= sum2/M;
    var[i] = sqrt((rN2[i] - pow(rN[i],4))/M);
    error[i] = sqrt(var[i]);
  
 }

  ofstream out;
  out.open("rN2_RWdiscrete.dat");
  for(int i=0; i<nstep;i++){
	out << i <<"  "<< rN[i] << "  " << error[i] << endl;
  }
  out.close();


//continuum


  for(int k=0; k < M; k++){
   
    for(int i=0; i < nstep; i++){

	theta[i][k] = acos(rnd.Rannyu(0,2)-1); //secondo distribution sin(x)/2
    	fi[i][k] = acos(rnd.Rannyu(0,2)-1);
	
	if(rnd.Rannyu() > 0.5){ //perchè fi deve essere estratto tra 0 e 2*pigreco
    		fi[i][k] = acos(rnd.Rannyu(0,2)-1) + M_PI;
	}		
    

      if(i == 0){
      	x_cont[i][k] = 0;
      	y_cont[i][k] = 0;
      	z_cont[i][k] = 0;
      }

      else{   //coordinate polari
      	x_cont[i][k] = x_cont[i-1][k] + a*cos(fi[i][k])*sin(theta[i][k]);
      	y_cont[i][k] = y_cont[i-1][k] + a*sin(fi[i][k])*sin(theta[i][k]);
      	z_cont[i][k] = z_cont[i-1][k] + a*cos(theta[i][k]);
      }
    }

  }


  for(int i=0; i < nstep; i++){
  sum = 0; 
  sum2 = 0;

   for(int k= 0; k < M; k++){
   sum+=pow(x_cont[i][k],2) + pow(y_cont[i][k],2) + pow(z_cont[i][k],2);
   sum2+=pow(pow(x_cont[i][k],2) + pow(y_cont[i][k],2) + pow(z_cont[i][k],2),2);
   }
   
  rN_cont[i] = sqrt(sum/M);
  rN2_cont[i]= sum2/M;
  var_cont[i] = sqrt((rN2_cont[i] - pow(rN_cont[i],4))/M);
  error_cont[i] = sqrt(var_cont[i]);
  
 }


  out.open("rN2_RWcontinuum.dat");
  for(int i=0; i<nstep;i++){
	out << i <<"  "<< rN_cont[i] << "  " << error_cont[i] << endl;
  }
  out.close();

return 0; 

}
   






















