/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main() //dobbiamo fare più simulazioni a diverse T (15/20 punti tra 0.5 e 2) -> valore finale del data blocking ti dà il punto per fare il grafico delle varie quantità in funzione di T 
{ 
  Input(); //Inizialization

   for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
   {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro); //Newman-Barckman sono citati altri metodi un po' più sofisticati ma non troppo
      Measure();
      Accumulate(); //Update block averages
    }
      Averages(iblk); //Print results for current block (iblk), mentre per ogni step (nstep)       
   }
   ConfFinal(); 

   if(irestart == 1){
   cout << "Printing results..." << endl;
   DatiTemp();
   }

  return 0;
}


void Input(void)
{
  ifstream ReadInput;
  ifstream ReadConf;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> num; //questo non prende la vera temperatura ma a che punto dell'intervallo [0.5,2] si è
  temp = 0.5 + num * (1.5/20.); //così viene calcolata la temperatura
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep; //nstep per BLOCCO! quind nstep tot= nblk*nstep

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;

  ReadInput >> irestart;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
//uso questi indici per riempire di quantità misurate il vettore walker e ogni quantità andrà a riempire l'indice rispettivo -> questo facilita la misura delle medie e delle incertezze
 
  n_props = 4; //Number of observables

//initial configuration


  if(irestart == 0){
   for (int i=0; i<nspin; ++i) //questa corrisponde a T= inf, spin totalmente casuali
   {
     if(rnd.Rannyu() >= 0.5) s[i] = 1;
     else s[i] = -1;
   }
  }
  else{
   ReadConf.open("config.final"); 
    cout << "Reading from config.final..." << endl << endl;
    for (int i=0; i<nspin; ++i){ //questa corrisponde a T= inf, spin totalmente casuali

     ReadConf >> s[i];
    }
   ReadConf.close();
  }

//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro) //è da completare 
{
  int o, random;
  double sm, acc, alpha, r;
  double prob_down, prob_up; 

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle
    o = (int)(rnd.Rannyu()*nspin); //c'è rischio di non chiamare alcuni indici, ma non è un problema perchè gli altri saranno chiamati (è cmq una funzione random)

    if(metro==1){ //Metropolis avrà bisogno della funzione Boltzmann (già implementata) che dà energia associata allo spin sm grazie agli spin primi vicini ip-1 e ip+1 con le adeguate PBC
    random = int(rnd.Rannyu(0,2)); //per decidere come viene flippato sm
      if(random == 0) sm = -1;
      else sm = 1;
    acc = exp(-(Boltzmann(sm,o)-Boltzmann(s[o],o))*pow(temp,-1));
      if(acc < 1) alpha = acc;
      else alpha = 1; 
    r = rnd.Rannyu(); 
      if(r <= alpha){
      s[o] = sm;  
      accepted++; 
      }
    attempted++;
    }

    else //Gibbs sampling: due prob (prob sia up che down, la cui somma fa 1 chiaramente) per decidere quale delle due: genero numero casuale tra 0 e 1, e vedi dove cade r, e poi accetti la config in cui il tuo r è caduto
      {
       prob_up = 1./(1 + exp(-pow(temp,-1)*(Boltzmann(-1,o)-Boltzmann(1,o)))); //prob che dopo lo switch sia down
       prob_down = 1./(1 + exp(-pow(temp,-1)*(Boltzmann(1,o)-Boltzmann(-1,o)))); //prob che dopo lo switch sia up
       r = rnd.Rannyu(); 
       if(prob_down <= prob_up){
	  if(r < prob_down) s[o] = -1; //ricade nella prob che sia down dopo lo switch
	  if(r>= prob_down) s[o] = 1; //ricade nella prob che sia up dopo lo switch
       }
       if(prob_down >= prob_up){
	  if(r < prob_up) s[o] = 1;
	  if(r>= prob_up) s[o] = -1;
       }

    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm; //sm è valore, ip è la posizione di quello spin sm
  return ene;
}

void Measure()
{
  //int bin?
  double u = 0.0, m = 0.0, u_2 = 0.0, x = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  } //consiglio: stimare in Averages la capacità termica e non qui, qui calcoliamo solo u^2
  x = pow(temp,-1)*pow(m,2);
  u_2 = pow(u,2);
  walker[iu] = u; //dobbiamo calcolare le altre componenti di walker
  walker[im] = m;
  walker[ic] = u_2;  
  walker[ix] = x; 
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i]; //valori dell'energia ai vari step che vengono sommati tra di loro     
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
    if(metro == 1){
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    }

    Ene.open("output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //valore medio dell'energy nel BLOCCO -> questi sono gli Ai (se voglio che sia energia per gdl devo dividere per nspin)
    glob_av[iu]  += stima_u; //qui si sommano i vari Ai 
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    Heat.open("output.heat.0",ios::app);
    stima_u2 = blk_av[ic]/blk_norm; //valore medio dell'energy per grado di libertà nel BLOCCO -> questi sono gli Ai
    stima_u = pow(stima_u*nspin,2); 
    stima_c = pow(temp,-2)*(stima_u2 - stima_u)/(double)nspin;
    glob_av[ic]  += stima_c; //qui si sommano i vari Ai 
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();

    Mag.open("output.mag.0",ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin; //valore medio dell'energy per grado di libertà nel BLOCCO -> questi sono gli Ai
    glob_av[im]  += stima_m; //qui si sommano i vari Ai 
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();

    Chi.open("output.chi.0",ios::app);
    stima_x = blk_av[ix]/blk_norm/(double)nspin; //valore medio dell'energy per grado di libertà nel BLOCCO -> questi sono gli Ai
    glob_av[ix]  += stima_x; //qui si sommano i vari Ai 
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();


    //cout << "----------------------------" << endl << endl;
}

void DatiTemp(){
 ifstream OpenDati;
 ofstream StampaDati;
 const int wd=12;
 double x[20], medie[20], glob_medie[20], errori[20];

 OpenDati.open("output.ene.0");
  for(int i = 0; i <20; i++){
  OpenDati >> x[i] >> medie[i] >> glob_medie[i] >> errori[i];
  }
 OpenDati.close();

 if(metro==1) StampaDati.open("ene_temp_metrop.txt",ios::app);
 else StampaDati.open("ene_temp_gibbs.txt",ios::app);
  StampaDati << temp << "    " << glob_medie[19] << "    " << errori[19] << endl;
 StampaDati.close();

 OpenDati.open("output.chi.0");
  for(int i = 0; i <20; i++){
  OpenDati >> x[i] >> medie[i] >> glob_medie[i] >> errori[i];
  }
 OpenDati.close();

 if(metro==1) StampaDati.open("chi_temp_metrop.txt",ios::app);
 else StampaDati.open("chi_temp_gibbs.txt",ios::app);
  StampaDati << temp << "    " << glob_medie[19] << "    " << errori[19] << endl;
 StampaDati.close();


 OpenDati.open("output.mag.0");
  for(int i = 0; i <20; i++){
  OpenDati >> x[i] >> medie[i] >> glob_medie[i] >> errori[i];
  }
 OpenDati.close();

 if(metro==1) StampaDati.open("mag_temp_metrop.txt",ios::app);
 else StampaDati.open("mag_temp_gibbs.txt",ios::app);
  StampaDati << temp << "    " << glob_medie[19] << "    " << errori[19] << endl;
 StampaDati.close();


 OpenDati.open("output.heat.0");
  for(int i = 0; i <20; i++){
  OpenDati >> x[i] >> medie[i] >> glob_medie[i] >> errori[i];
  }
 OpenDati.close();

 if(metro==1) StampaDati.open("heat_temp_metrop.txt",ios::app);
 else StampaDati.open("heat_temp_gibbs.txt",ios::app);
  StampaDati << temp << "    " << glob_medie[19] << "    " << errori[19] << endl;
 StampaDati.close();

}



void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
    if(irestart==0) cout << s[i] << endl; 
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk) 
{ 
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk); 
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
