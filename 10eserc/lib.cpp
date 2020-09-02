#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "lib.h"
#include "random.h"

using namespace std;

//Path :: Path(){}

Path :: Path(Random *rnd, int which){

  _rand = rnd;
  n_city = 32; 
  
  	path = new int [n_city];
	old_path = new int[n_city];
	x = new double[n_city];
	y = new double[n_city];
	L = 0;


	if(which==0){
		for(int i=0;i<n_city;i++){
			double theta=_rand->Rannyu(0,2*M_PI);
				x[i]=cos(theta);
				y[i]=sin(theta);
		}
	}


	else{
		for(int i=0;i<n_city;i++){
			x[i]=_rand->Rannyu(0.,2.);
			y[i]=_rand->Rannyu(0.,2.);
		}
	}

		
	for(int i=0;i<n_city;i++){
		path[i]=i+1;   
	}


	RandomPath();
	Check();

	L = ComputeL();		

}



	

void Path :: RandomPath() {

 
	int step=(int)_rand->Rannyu(1.,30.);	//genero casualmente il punto del vettore in cui devo scambiare il valore selezionato è quello successivo.
	int r=(int)_rand->Rannyu(1.,31-step); 
	int a=path[r];
	int b=path[r+step];
	path[r]=b;
	path[r+step]=a;

return;
}



void Path :: Check(){

 int wrong=0;

	for(int i = 0; i < n_city-1; i++){
	  for(int j = i+1; j < n_city; j++){
	   if(path[i] == path[j]) wrong++;
	  }
	}

	if(wrong != 0) cout << "The salesman visits a city more than one time!" << endl; 
        
	if(path[0] != 1) cout << "The first city the salesman visits is not number 1! " << endl; 
}


double Path :: ComputeL(){

  double sum = 0; 
	double dx,dy;
	for(int i=0;i<n_city;i++){
		if(i==n_city-1){
			dx=x[path[i]-1]- x[path[0]-1];	
			dy=y[path[i]-1]- y[path[0]-1];
		}
		else{
			dx=x[path[i]-1]-x[path[i+1]-1];	
			dy=y[path[i]-1]-y[path[i+1]-1];	//idem per l'asse y					
		}
		sum+=pow((pow(dx,2.)+pow(dy,2.)),0.5);		//calcolo la distanza tra le due città e le sommo alle sistanze precedenti-->funzione di costo
	}

  return sum;
}


void Path :: Shift(){ 

  int m = int(_rand->Rannyu(2,31));
  int n = int(_rand->Rannyu(1,10));
  int value=0;

  if(n+m<n_city-1){ //per evitare di andare oltre a 32
	if(n>=m){
  		 for(int i=m; i>0;i--){
	  		value = path[i+n];
	  		path[i+n] = path[i];
	  		path[i] = value;
	 	 }
	}
	else{
  	 	 for(int i=m; i>1;i--){
	  		value = path[i+n];
	  		path[i+n] = path[i];
	  		path[i] = value;
	 	 }
	 	 for(int i=1; i<n+1;i++){ //una volta arrivata al 2 devo traslare tutto di uno 
	  		value = path[i+1];
	  		path[i+1] = path[i];
	  		path[i] = value;
	 	 }
	}
  }
	
 return;
}

void Path :: Permut (){

  	int m =(int)_rand->Rannyu(1.,16);
	int value = 0;
	for(int j=1;j<=m;j++){
		value=path[m+j];
		path[m+j]=path[j];
		path[j]=value;
	}
 return;
}

void Path :: Inversion(){

	int m=(int)_rand->Rannyu(1.,n_city);
	int value=0;
	if(m%2==0){
		int cont=1;	
		for(int j=m;j>m/2;j--){
			value=path[j];
			path[j]=path[cont];
			path[cont]=value;
			cont=cont+1;
		}
	}
	else{
		int cont=1;
		for(int j=m;j>(m+1)/2;j--){
			value=path[j];
			path[j]=path[cont];
			path[cont]=value;
			cont=cont+1;
		}
	}
 return;
}

void Path :: Mutazione(double temp){

    for(int k=0; k<4;k++){

	for(int i=0; i< n_city;i++){
		old_path[i] = path[i];
	}

	double old_L = L; 

	if(k==0) RandomPath();
	if(k==1) Shift();
	if(k==2) Permut();
	if(k==3) Inversion();

	double new_L = ComputeL(); 

	double p = exp(-(new_L-old_L)/temp);
	if(p < _rand->Rannyu()){
		for(int i=0; i<n_city;i++){	
			path[i] = old_path[i];
		}
	} 

	L = ComputeL(); 
    }

    Check();

 return;
 			
}


void Path:: BestPath(int which){
	ofstream Out;
	if(which==0){
		Out.open("bestpath_circ.out");
	}
	else{
		Out.open("bestpath_quadr.out");
	}
	for(int i=0;i<n_city;i++){
		Out<< path[i] << "  " << x[path[i]-1]<<" "<<y[path[i]-1]<<endl;
	}
	Out<< path[0] << "  " << x[path[0]-1]<<" "<<y[path[0]-1]<<endl;
	Out.close();
	return;
}
 
void Path:: outL(int which, double temp){
	ofstream Out;
	if(which==0){
		Out.open("L_circ.out",ios::app);
	}
	else{
		Out.open("L_quadr.out",ios::app);
	}
	Out << temp << "  " << L << "  " << endl;

	Out.close();
	return;
}

		

