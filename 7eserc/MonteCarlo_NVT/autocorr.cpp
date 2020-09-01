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
  
  double *valore_pres=new double[M];
  double *valore_epot=new double[M]; 
  double *energy_autocorr=new double[1000];
  double *pression_autocorr=new double[1000];

  ifstream in; 
  ofstream out;

  int phase=0;  

  cout << "Of which phase do you want to calculate the autocorrelation?" << endl; 
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

  for (int t=0; t<1000;t++){ 
        energy_autocorr[t] = autocorrelation(valore_epot,t, M); //t èla distanza in step tra due misure
        pression_autocorr[t] = autocorrelation(valore_pres,t,M);
  }

  if(phase==0) out.open("risultati/sol_press_autocor.out"); 
  if(phase==1) out.open("risultati/liqu_press_autocor.out"); 
  if(phase==2) out.open("risultati/gas_press_autocor.out");  
  for(int i=0; i<1000; i++){
	out << pression_autocorr[i] << endl; 
  }
  out.close();


  if(phase==0) out.open("risultati/sol_epot_autocor.out"); 
  if(phase==1) out.open("risultati/liqu_epot_autocor.out"); 
  if(phase==2) out.open("risultati/gas_epot_autocor.out");   
  for(int i=0; i<1000; i++){
	out << energy_autocorr[i] << endl; 
  }
  out.close();

  delete [] valore_pres;
  delete [] valore_epot; 

return 0; 

}
