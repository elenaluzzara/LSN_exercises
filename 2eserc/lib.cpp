#include <iostream> 
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include "lib.h"

using namespace std;

BlockingMethod :: BlockingMethod(){

  sum=0;
  nblk=0;
  nstep=0;
  sum_prog = new double[nblk];
  su2_prog = new double[nblk];
  error_prog = new double[nblk];
  ave = new double[nblk];
  av2 = new double[nblk];
  k=0;

  for(int i=0; i< nblk; i++){
	sum_prog[i] = 0;
	su2_prog[i] = 0;
	error_prog[i] = 0;
	ave[i] = 0;
	av2[i] = 0;
  }
}

BlockingMethod :: BlockingMethod(int N,int L){

  sum=0;
  nblk=N;
  nstep=L;
  sum_prog = new double[nblk];
  su2_prog = new double[nblk];
  error_prog = new double[nblk];
  ave = new double[nblk];
  av2 = new double[nblk];
  k=0;

  for(int i=0; i< nblk; i++){
	sum_prog[i] = 0;
	su2_prog[i] = 0;
	error_prog[i] = 0;
	ave[i] = 0;
	av2[i] = 0;
  }
}

BlockingMethod :: ~BlockingMethod(){}


void BlockingMethod :: BlkMethod(double *casual){

   for(int i=0; i<nblk; i++){
	sum = 0; 
	for (int j=0; j<nstep; j++){
		k =j+i*nstep;
		sum+=casual[k];
	}
	ave[i] = sum/nstep; 
	av2[i] = pow(ave[i],2);
  }

  for(int i =0; i < nblk; i++){
	for(int j =0; j < i+1; j++){
		sum_prog[i] += ave[j];
		su2_prog[i] += av2[j];
	}
	sum_prog[i]/=(i+1);
	su2_prog[i]/=(i+1);
	error_prog[i] = ErrBlk(i);
  }


 return;
}

void BlockingMethod :: Stampa(const char* nomefile){

  ofstream out;
  out.open(nomefile);
  for(int i=0; i<nblk;i++){
	out << i*nstep <<"  "<< sum_prog[i] << "  " << error_prog[i] << endl;
  }
	
  out.close();
return;
}



void BlockingMethod :: BlkMethod_pigreco(const char* nomefile){ //per il calcolo con il pigreco dove c'è già la media di pigreco calcolata con L throws per N volte

   ofstream out;

   for(int i=0; i<nblk; i++){
	av2[i] = pow(ave[i],2);
  }
  
  for(int i =0; i < nblk; i++){
	for(int j =0; j < i+1; j++){
		sum_prog[i] += ave[j];
		su2_prog[i] += av2[j];
	}
	sum_prog[i]/=(i+1);
	su2_prog[i]/=(i+1);
	error_prog[i] = ErrBlk(i);
  }

  out.open(nomefile);
  for(int i=0; i<nblk;i++){
	out << i*nstep <<"  "<< sum_prog[i] << "  " << error_prog[i] << endl;
  }
	
  out.close();
 return;

}

double BlockingMethod :: ErrBlk(int iblk){

  double varianza = 0;
  if (iblk == 0)
	return 0; 
  else {
  varianza = sqrt((su2_prog[iblk] - pow(sum_prog[iblk],2))/iblk);
  return varianza;
  }

}

//

void data_chi2(const char* nomefile,  int x[], double dati[], int dime)
{
   ofstream out;
   out.open(nomefile);
   if (out.is_open()){
	for (int i = 0; i < dime; i++){
        out << x[i] << " " << dati[i] << endl; 
   	} 
   } else cerr << "PROBLEM: Unable to open file.out" << endl;

  out.close();
  return;
}

void data_histo(const char* nomefile,  int x[], double dati0[], double dati1[], double dati2[], double dati3[], int dime)
{
   ofstream out;
   out.open(nomefile);
   if (out.is_open()){
	for (int i = 0; i < dime; i++){
        out << x[i] << " " << dati0[i] << " " << dati1[i] << " " << dati2[i] << " " << dati3[i] << endl; 
   	} 
   } else cerr << "PROBLEM: Unable to open file.out" << endl;

  out.close();
  return;
}

