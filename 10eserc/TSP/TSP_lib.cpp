#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>
#include "TSP_lib.h"
#include "random.h"

using namespace std;

//Path :: Path(){}

Path :: Path(Random *rnd,int rank){

   n_city = 32;
  _rand = rnd;

  //InitializeRank(rank);
  for(int i=0; i < 100; i++){
  	path[i] = new int [n_city];
	old_path[i] = new int[n_city]; 
  }

  x = new double[n_city];
  y = new double[n_city];
  L = new double[100];

 // InitialPathRank(); //così è uguale per tutti i rank
  ifstream in("config.in"); 
  if(in.is_open()){
	for(int i=0; i<n_city;i++){
		in >> x[i] >> y[i];
	}
  } else cerr << "Problem in opening the initial path" << endl; 
  in.close();

  for(int k=0;k<100;k++){
	for(int i=0;i<n_city;i++){
		path[k][i]=i+1;     //crea 100 percorsi tutti uguali
	}
  }

 for(int j=0; j<100;j++){
	RandomPath(j);				//modifica i percorsi creati
 } 

  Check();

  for(int i=0;i<100;i++){
	L[i]=ComputeL(i); //calcolo distanza 
  }
	
  Ordina(); //ordino da L minore a L maggiore		
		
}

/*void Path::InitializeRank(int rank) {

   int p1, p2;
   ifstream Primes("Primes");
        if (Primes.is_open()){
                Primes >> p1>> p2;
        } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   int seed[4]; //questa prima parte va fatto sempre per implementare i semi di completamento
   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            _rand.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
//cout << seed[0] << " " << seed[1] << " " << seed[2] << endl; 
//cout << "IO sono il nodo" << rank << "e il mio numero casuale è" << (int)rnd->Rannyu(1.,31.) << endl;
 return;

}*/

void Path :: RandomPath(int pa) {


	/*int o = int(_rand->Rannyu(1,31)); 
	int value = path[pa][o]; 
	int step = int(_rand->Rannyu(1,30)); 

	path[pa][o] = path[pa][Pbc(o+step)];
	path[pa][Pbc(o+step)] = value;*/	
 
	int r=int(abs(_rand->Rannyu())*(29)+1);	//genero casualmente il punto del vettore in cui devo scambiare il valore selezionato è quello successivo.
	int a=path[pa][r];
	int b=path[pa][r+1];
	path[pa][r]=b;
	path[pa][r+1]=a;

return;
}

/*void Path :: SetPath(int * value) {

  	for(int i=0; i < n_city; i++){
        path[i] = value[i]; 
	}	
 
return;
}*/


void Path :: Check(){

 int wrong=0;

     for(int k=0; k <100; k++){
	wrong=0;
	for(int i = 0; i < n_city-1; i++){
	  for(int j = i+1; j < n_city; j++){
	   if(path[k][i] == path[k][j]) wrong++;
	  }
	}
     }

	if(wrong != 0) cout << "The salesman visits more than one city at a time!" << endl; 

     for(int k=0; k<100;k++){
        if(path[k][0] != 1){
		 cout << "The first city the salesman visits is not number 1! " << endl;
		for(int i=0; i <n_city;i++){
			cout << path[k][i] << "  ";
		}
		cout << endl;  
	}
     }

return; 
}


double Path :: ComputeL(int pa){

  double sum = 0; 
	double dx,dy;
	for(int i=0;i<n_city;i++){
		if(i==n_city-1){
			dx=x[path[pa][i]-1]- x[path[pa][0]-1];	
			dy=y[path[pa][i]-1]- y[path[pa][0]-1];
		}
		else{
			dx=x[path[pa][i]-1]-x[path[pa][i+1]-1];	//calcolo la distanza sull'asse x dei due indici consecutivi del percorso numero "ind"
			dy=y[path[pa][i]-1]-y[path[pa][i+1]-1];	//idem per l'asse y					
		}
		sum+=pow((pow(dx,2.)+pow(dy,2.)),0.5);		//calcolo la distanza tra le due città e le sommo alle sistanze precedenti-->funzione di costo
	}

  return sum;
}

void Path :: Ordina (){

 int scambio;
 double value=0;
 int ordin[100];
 double disordineL[100];
 int ordinepath[100][n_city];
 scambio=1;

	for(int i=0;i<100;i++){
		disordineL[i]=L[i];   
	}

	while(scambio) {
		scambio = 0;
		for(int i=0; i<99; i++) {
			if( L[i] > L[i+1]) {
				value = L[i];
				L[i] = L[i+1];
				L[i+1] = value;
				scambio = 1;
			}
		}
	}


	for(int i=0; i<100;i++){
		for(int k =0; k<100;k++){
			if(L[i] == disordineL[k]){
				ordin[i] = k; 
			}
		}
	}
	

	for(int i =0; i<100; i++){
		int ind = ordin[i];
		for(int j=0; j<n_city;j++){
			ordinepath[i][j] = path[ind][j];
		}
	}

	for(int i =0; i<100; i++){
		for(int j=0; j<n_city;j++){
			path[i][j] = ordinepath[i][j];
		}
	}
  	
	for(int i =0; i<100;i++){
	 L[i] = ComputeL(i);
	}
 return;
}


void Path :: Shift(int pa){ 

  int m = int(abs(_rand->Rannyu())*(31-2)+2);
  int n = int(abs(_rand->Rannyu())*(10-1)+1);
  int value=0;

  if(n+m<n_city-1){ //per evitare di andare oltre a 32
	if(n>=m){
  		 for(int i=m; i>0;i--){
	  		value = path[pa][i+n];
	  		path[pa][i+n] = path[pa][i];
	  		path[pa][i] = value;
	 	 }
	}
	else{
  	 	 for(int i=m; i>1;i--){
	  		value = path[pa][i+n];
	  		path[pa][i+n] = path[pa][i];
	  		path[pa][i] = value;
	 	 }
	 	 for(int i=1; i<n+1;i++){ //una volta arrivata al 2 devo traslare tutto di uno 
	  		value = path[pa][i+1];
	  		path[pa][i+1] = path[pa][i];
	  		path[pa][i] = value;
	 	 }
	}
  }
	
 return;
}

void Path :: Permut (int pa){
/*
  int o = int(_rand->Rannyu(1,16-m)); // m = 14, pesco 1
  int u = int(_rand->Rannyu(15,33-m)); // pesco 14

	for(int i = 0; i < m; i++){
	  int value = path[tap][1+i];
	  path[tap][1+i] = path[tap][m+1+i]; //1 con 14
	  path[tap][(1+m)+i] = value;
	}*/

  	int m =int(abs(_rand->Rannyu())*(16-1)+1);
	int value = 0;
	for(int j=1;j<=m;j++){
		value=path[pa][m+j];
		path[pa][m+j]=path[pa][j];
		path[pa][j]=value;
	}
 return;
}

void Path :: Inversion(int pa){

  //int o = int(_rand->Rannyu(1,33-m)); //m =31, o=(1,2), quindi 1
	/*int m=(int)_rand->Rannyu(2,n_city-1);
	int value=0;
	for(int i = 0; i < m; i++){
	  int value = path[pa][1+i];
	  path[pa][1+i] = path[pa][1+(m-1)-i]; //1+30
	  path[pa][1+(m-1)-i] = value;
	}*/

	int m=int(abs(_rand->Rannyu())*(32-1)+1);
	int value=0;
	if(m%2==0){
		int cont=1;	
		for(int j=m;j>m/2;j--){
			value=path[pa][j];
			path[pa][j]=path[pa][cont];
			path[pa][cont]=value;
			cont=cont+1;
		}
	}
	else{
		int cont=1;
		for(int j=m;j>(m+1)/2;j--){
			value=path[pa][j];
			path[pa][j]=path[pa][cont];
			path[pa][cont]=value;
			cont=cont+1;
		}
	}
 return;
}

void Path :: Mutazione(){

	for(int i=0;i<100;i++){
		if(abs(_rand->Rannyu())<0.1){ RandomPath(i);} //faccio la mutazione Pari sul percorso i-esimo
		if(abs(_rand->Rannyu())<0.1){ Shift(i);}
		if(abs(_rand->Rannyu())<0.1){ Permut(i);}
		if(abs(_rand->Rannyu())<0.1){ Inversion(i);}
	}
	
	Check();

	for(int i=0;i<100;i++){
		L[i]=ComputeL(i);
	}
	
	Ordina();	

}

void Path::SetBestPath(int * bp){

	for(int i=0; i<n_city;i++){
		path[0][i] = bp[i]; 
	}
 return; 
}

void Path::BestL(int rank){
	ofstream BL;

	BL.open("bestL_quadr"+to_string(rank)+".out",ios::app);
	BL << L[0] << "  ";
	for(int i=0;i<n_city;i++){
		BL<<path[0][i]<< "  ";
	}
	BL << endl;
	BL.close();
	return;
}

int * Path::ReadBestPath(){
 return path[0];
}


void Path:: InitialPath(int rank){
	ofstream Out;
	Out.open("initialpath_quadr"+to_string(rank)+".out");
	for(int i=0;i<n_city;i++){
		Out<< path[0][i] << "  " << x[path[0][i]-1]<<" "<<y[path[0][i]-1]<<endl;
	}
	Out<< path[0][0] << "  " << x[path[0][0]-1]<<" "<<y[path[0][0]-1]<<endl;
	Out.close();
	return;
}
 

void Path:: BestPath(int rank){
	ofstream Out;
	Out.open("bestpath_quadr"+to_string(rank)+".out");
	for(int i=0;i<n_city;i++){
		Out<< path[0][i] << "  " << x[path[0][i]-1]<<" "<<y[path[0][i]-1]<<endl;
	}
	Out<< path[0][0] << "  " << x[path[0][0]-1]<<" "<<y[path[0][0]-1]<<endl;
	Out.close();
	return;
}
 

int Path :: SelectOperator(int p){

 int j = 0; 
 double r = abs(_rand->Rannyu()); 
 
  j = int(100*pow(r,p)); 
 //cout << "roulette truccata: " << j << endl;
 return j; 
}

void Path :: Crossover(int mother, int father, int ind){//int idx, int idx2,int cont, int cont2){

  int cut = int(abs(_rand->Rannyu())*(30-1)+1); //punto fino al quale viene salvato il vettore path

  int m[n_city];
  int f[n_city];

  int assent = 0;
  int indice = 0;
  //cout << " numero: " << exchange << endl;
	for(int i = 0; i < n_city; i++){
	  if(i < cut){
	  m[i] = path[mother][i]; 
	  f[i] = path[father][i];
	  }
	  else{
	  m[i] = -9;
          f[i] = -9;
          }
	}


	for(int j = 0; j < n_city; j++){
     	  assent = 0;
		for(int i = 0; i < cut; i++){ //se o = 1 non fa questo for
			if(m[i] == path[father][j]) break; 
			else assent++;
		}
		
		if(assent == cut) {		
		  m[cut+indice] = path[father][j];
		  indice++;
		}
	}

  indice = 0;

	for(int j = 0; j < n_city; j++){
     	  assent = 0;
		for(int i = 0; i < cut; i++){
			if(f[i] == path[mother][j]) break; 
			else assent++;
		}
		
		if(assent == cut) {		
		  f[cut+indice] = path[mother][j];
		  indice++;
		}
	}

	for(int i =0; i<n_city;i++){
	  old_path[ind][i] = f[i];
	  old_path[ind+1][i] = m[i];
	}

 return;
}
 

void Path :: NuovaPop (){

	//int save[2][n_city];

	for(int i=0; i<100;i++){
		for(int k=0; k<n_city;k++){
			old_path[i][k] = path[i][k];
		}
	}

	for(int i=0; i < 99;i=i+2){
		int father = SelectOperator(2); //mi restituisce un certo numero j che seleziona il percorso 
		int mother = SelectOperator(2);

		if(abs(_rand->Rannyu()) < 0.7) Crossover(mother,father,i); 

	 }

	for(int i=0; i<100;i++){
		for(int k=0; k<n_city;k++){
			path[i][k] = old_path[i][k];
		}
	}

	Check();

	for(int i=0; i < 100; i++){
  	  L[i] = ComputeL(i);
	}

  Ordina(); 

 return;

}

