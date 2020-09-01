#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include "lib.h"

using namespace std;

int main()
{ 

  int M = 500000; 
  int N[101];   
  int L[101]; 

  double *valore_pres=new double[M];
  double *valore_epot=new double[M]; 
  double *error_blk_press=new double[101];
  double *error_blk_epot=new double[101];

  BlockingMethod blk_press[101];
  BlockingMethod blk_epot[101];

  ifstream in; 
  ofstream out;

  int phase=0;  

  cout << "Of which phase do you want to calculate the uncertainties as a function of the number of block?" << endl; 
  cout << "Type 0 for solid, 1 for liquid and 2 for gas" << endl;
  cin >> phase;

  if(phase==0) in.open("risultati/sol_instant_pres.out"); 
  if(phase==1) in.open("risultati/liqu_instant_pres.out"); 
  if(phase==2) in.open("risultati/gas_instant_pres.out"); 
  for(int i=0; i < M; i++){
	in >> valore_pres[i];
  }
  in.close(); 

  if(phase==0) in.open("risultati/sol_instant_epot.out"); 
  if(phase==1) in.open("risultati/liqu_instant_epot.out"); 
  if(phase==2) in.open("risultati/gas_instant_epot.out"); 
  for(int i=0; i < M; i++){
	in >> valore_epot[i];
  }
  in.close(); 


  for(int i=0; i<101; i++){
	L[i]= 10 + i*49.9; 
	N[i] = M/L[i]; 
  	blk_press[i].Set_nstep(L[i]);
  	blk_press[i].Set_nblk(N[i]);
  	blk_epot[i].Set_nstep(L[i]);
  	blk_epot[i].Set_nblk(N[i]);
  }

  for(int i=0; i<101; i++){
	int nblk = N[i] -1;
	blk_press[i].BlkMethod(valore_pres); //blocking method
	blk_epot[i].BlkMethod(valore_epot);
	error_blk_press[i] = blk_press[i].Get_ErrorProg()[nblk]; 
	error_blk_epot[i] = blk_epot[i].Get_ErrorProg()[nblk];
  }

  if(phase==0) out.open("risultati/sol_blkave_press.out"); 
  if(phase==1) out.open("risultati/liqu_blkave_press.out"); 
  if(phase==2) out.open("risultati/gas_blkave_press.out");  
  for(int i=0; i<101; i++){
	out << L[i] << "  " << error_blk_press[i] << endl; 
  }
  out.close();


  if(phase==0) out.open("risultati/sol_blkave_epot.out"); 
  if(phase==1) out.open("risultati/liqu_blkave_epot.out"); 
  if(phase==2) out.open("risultati/gas_blkave_epot.out");   
  for(int i=0; i<101; i++){
	out << L[i] << "  " << error_blk_epot[i] << endl; 
  }
  out.close();

  delete [] valore_pres;
  delete [] valore_epot; 

return 0; 

}




  

