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
  
  	for(int i=0; i < 100; i++){
  	path[i] = new int [n_city];
	old_path[i] = new int[n_city]; 
	}

	x = new double[n_city];
	y = new double[n_city];
	L = new double[100];


	if(which==0){ //circo di raggio 1
		for(int i=0;i<n_city;i++){
			double theta=_rand->Rannyu(0,2*M_PI);
				x[i]=cos(theta);
				y[i]=sin(theta);
		}
	}

	else{ //quadrato di lato 2
		for(int i=0;i<n_city;i++){
			x[i]=_rand->Rannyu(0.,2.);
			y[i]=_rand->Rannyu(0.,2.);
		}
	}


	for(int k=0;k<100;k++){
		for(int i=0;i<n_city;i++){
			path[k][i]=i+1;     //crea 100 percorsi tutti uguali
		}
	}

	for(int j=1;j<100;j++){
		RandomPath(j);				//modifica i percorsi creati
	}

	Check();
	for(int i=0;i<100;i++){
		L[i]=ComputeL(i); //calcolo distanza 
	}
	
	Ordina(); //ordino da L minore a L maggiore		
		
}



void Path :: RandomPath(int pa) {

	int step=(int)_rand->Rannyu(1.,30.);	//genero casualmente il punto del vettore in cui devo scambiare il valore selezionato è quello successivo.
	int r=(int)_rand->Rannyu(1.,31-step); //così che r + step non vada oltre 31
	int a=path[pa][r];
	int b=path[pa][r+step];
	path[pa][r]=b;
	path[pa][r+step]=a;

return;
}


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
        if(path[k][0] != 1) cout << "The first city the salesman visits is not number 1! " << endl; 
     }
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

	while(scambio) { //ordino il vettore L dei 100 individui
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
			if(L[i] == disordineL[k]){ //assegno agli L disordinato il numero del loro ordine
				ordin[i] = k; 
			}
		}
	}
	

	for(int i =0; i<100; i++){ //metto in un vettore i path ordinati
		int ind = ordin[i];
		for(int j=0; j<n_city;j++){ 
			ordinepath[i][j] = path[ind][j];
		}
	}

	for(int i =0; i<100; i++){ //riordino i path nella classe secondo l'ordine di L
		for(int j=0; j<n_city;j++){
			path[i][j] = ordinepath[i][j];
		}
	}
  	

	 SetL(L); 		//assegno a ogni path ordinato la sua L ordinata

 return;
}


void Path :: Shift(int pa){  

  int m = int(_rand->Rannyu(2,31)); //parto dalla numero due e scelgo m città contigue 
  int n = int(_rand->Rannyu(1,10));
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
  	 	 for(int i=m; i>1;i--){ // mi fermo alla seconda città
	  		value = path[pa][i+n];
	  		path[pa][i+n] = path[pa][i];
	  		path[pa][i] = value;
	 	 }
	 	 for(int i=1; i<n+1;i++){ //una volta arrivata al 2 devo traslare tutto di uno fino a un numero di città pari a n
	  		value = path[pa][i+1];
	  		path[pa][i+1] = path[pa][i];
	  		path[pa][i] = value;
	 	 }
	}
  }
	
 return;
}

void Path :: Permut (int pa){

  	int m =(int)_rand->Rannyu(1.,16);
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

	int m=(int)_rand->Rannyu(1.,n_city);
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
		if(_rand->Rannyu()<0.1){ RandomPath(i);} //faccio la mutazione Pari sul percorso i-esimo
		if(_rand->Rannyu()<0.1){ Shift(i);}
		if(_rand->Rannyu()<0.1){ Permut(i);}
		if(_rand->Rannyu()<0.1){ Inversion(i);}
	}
	
	Check();

	for(int i=0;i<100;i++){
		L[i]=ComputeL(i);
	}
	
	Ordina();	

}

void Path::BestL(int which){
	ofstream BL;
	if(which==0){
		BL.open("bestL_circ.out",ios::app);
	}
	else{
		BL.open("bestL_quadr.out",ios::app);
	}
	BL << L[0] << "  ";
	for(int i=0;i<n_city;i++){
		BL<<path[0][i]<< "  ";
	}
	BL << endl;
	BL.close();
	return;
}

void Path::AverageL(int which){

	double sum=0;
	ofstream AL;
	if(which==0){
		AL.open("averageL_circ.out",ios::app);
	}
	else{
		AL.open("averageL_quadr.out",ios::app);
	}

	for(int i=0; i<50;i++){
		sum += L[i];
	}
	AL << sum/50 << endl;
	AL.close();
  return;
}

void Path:: InitialPath(int which){
	ofstream Out;
	if(which==0){
		Out.open("initialpath_circ.out");
	}
	else{
		Out.open("initialpath_quadr.out");
	}
	for(int i=0;i<n_city;i++){
		Out<< path[0][i] << "  " << x[path[0][i]-1]<<" "<<y[path[0][i]-1]<<endl;
	}
	Out<< path[0][0] << "  " << x[path[0][0]-1]<<" "<<y[path[0][0]-1]<<endl;
	Out.close();
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
		Out<< path[0][i] << "  " << x[path[0][i]-1]<<" "<<y[path[0][i]-1]<<endl;
	}
	Out<< path[0][0] << "  " << x[path[0][0]-1]<<" "<<y[path[0][0]-1]<<endl;
	Out.close();
	return;
}
 

int Path :: SelectOperator(int p){

 int j = 0; 
 double r = _rand->Rannyu(); 
 
  j = int(100*pow(r,p)); 

 return j; 
}

void Path :: Crossover(int mother, int father, int ind){//int idx, int idx2,int cont, int cont2){

  int cut = int(_rand->Rannyu(1,30)); //punto fino al quale viene salvato il vettore path

  int m[n_city] = {0};
  int f[n_city] = {0};

  int assent = 0;
  int indice = 0;

	for(int i = 0; i < n_city; i++){
	  if(i < cut){ // m e f sono vettori che contengono solo la prima parte dei due vettori mother e father scelti
	  m[i] = path[mother][i]; 
	  f[i] = path[father][i];
	  }
	  else{
	  m[i] = -9;
          f[i] = -9;
          }
	}

//per ricostruire l'altro pezzo del vettore m dal vettore father
	for(int j = 0; j < n_city; j++){
     	  assent = 0;
		for(int i = 0; i < cut; i++){ //così capisco in che ordine sono presenti le città nel vettore father e nello stesso ordine le mette nel pezzo mancante del vettore m
			if(m[i] == path[father][j]) break; //se una città del vettore father è già presente nell'inizio del vettore m la escludo 
			else assent++;
		}
		
		if(assent == cut) { // se quella città j-esima del vettore father è assente dal vettore iniziale m che ha cut elementi allora la aggiungo
		  m[cut+indice] = path[father][j];
		  indice++; //così tengo conto a che punto sono arrivata a riempire la seconda parte del vettore m
		}
	}

//stesso ragionamento per il riempimento della seconda parte del vettore f dal vettore mother
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

//metto i nuovi vettori così creati in old_path
	for(int i =0; i<n_city;i++){
	  old_path[ind][i] = f[i];
	  old_path[ind+1][i] = m[i];
	}

 return;
}
 

void Path :: NuovaPop (){


	for(int i=0; i<100;i++){
		for(int k=0; k<n_city;k++){
			old_path[i][k] = path[i][k];
		}
	}

	for(int i=0; i < 99;i=i+2){
		int father = SelectOperator(2); //mi restituisce un certo numero j che seleziona il percorso 
		int mother = SelectOperator(2);

		if(_rand->Rannyu() < 0.7) Crossover(mother,father,i); 

	 }

//ho creato i nuovi vettori e li rimetto in path
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





