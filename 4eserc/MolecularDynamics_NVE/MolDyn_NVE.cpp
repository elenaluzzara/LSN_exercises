/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"


using namespace std;

int main(){ 
  Input();             //Inizialization
  int nconf = 1;

  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
        Move();  
      if(istep%10 == 0){
        Measure();
        Accumulate(); //Update block averages
        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
	nconf += 1;
      }
    }
    Averages(iblk);  
  }
 
  OldConf();
  ConfFinal();

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input
 
  ReadInput >> irestart;
  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl; //larghezza del box in unità di LJ

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> iprint;


  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

  igofr = 4;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;
  cout << "Dimension of the bin for g(r): " << bin_size << endl;

//Read initial configuration

 if(irestart == 1){
  cout << "Read initial configuration from file old.0" << endl << endl;
  ReadConf.open("old.0");

  if (ReadConf.is_open()){

    for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i]; 
    x[i] = x[i] * box; // in unità di LJ
    y[i] = y[i] * box;
    z[i] = z[i] * box;
    }
  } else cerr << "PROBLEM: Unable to open old.0" << endl;

  ReadConf.close();
 } 

 else{
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box; // in unità di LJ
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();
 }


//Prepare old position
 if(irestart == 1){
  
  cout << "Prepare old position from the previous simulation" << endl << endl;
  ReadConf.open("old.final");

  if (ReadConf.is_open()){

    for (int i=0; i<npart; ++i){
    ReadConf >> xold[i] >> yold[i] >> zold[i]; 
    xold[i] = xold[i] * box; // in unità di LJ
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
    }
  } else cerr << "PROBLEM: Unable to open old.final" << endl;

  ReadConf.close();

  Move();   //faccio uno step di Verlet per trovare r(t+dt)

  double sum_v = 0;

  for(int i = 0; i < npart; ++i){
 
    vx_2[i] = Pbc(x[i] - xold[i])/(delta); //vx_2 ecc servono solo per calcolare T(t + dt/2) e trovare così il fattore con il quale riscalare le velocità vere al tempo t con cui calcora Ekin ecc
    vy_2[i] = Pbc(y[i] - yold[i])/(delta);
    vz_2[i] = Pbc(z[i] - zold[i])/(delta);

  sum_v += vx_2[i]*vx_2[i] + vy_2[i]*vy_2[i] + vz_2[i]*vz_2[i];
  }

   sum_v /= (double)npart;
   double T_2 = sum_v/3;  

   double fs = sqrt(temp/T_2);   // fs = velocity scale factor 

   for (int i=0; i<npart; ++i){
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;

    xold[i] = Pbc(x[i] - vx[i] * delta); 
    yold[i] = Pbc(y[i] - vy[i] * delta);
    zold[i] = Pbc(z[i] - vz[i] * delta);

   } 
 }
 
//Prepare initial velocities
 else{
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){ //veloctà scelte in modo casuale tra -0.5 e 0.5
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0]; //così non c'è nessun drift del momento iniziale del sistema -> il momento totale di tutte le direzioni è zero, non deve esserci una preferenza iniziale 
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta); //ottengo così le posizioni allo step precedente -> non è proprio soddisfacente e crea all'inizio un salto nell'Etot -> noi non vogliamo ripartire da zero ogni volta ma dalle posizioni precedentemente ottenute con la simulazione prima
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);

   } 
 }
 
   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) ); //non c'è divisione per massa perchè è in unità di misura della massa stessa quindi m = 1
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta); 
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  
  return;
}


void OldConf(void){

  ofstream OldConf;

  OldConf.open("old.0");

  for (int i=0; i<npart; ++i){
    OldConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl; 
  }
  OldConf.close();

  OldConf.open("old.final");

  for (int i=0; i<npart; ++i){
    OldConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl; 
  }
  OldConf.close();

  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){ //ciclo su tutte le particelle
    if(i != ip){ //esclusa la particella stessa
      dvec[0] = Pbc( x[ip] - x[i] );  // distanza tra la particella e le altre 
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){ // se il modulo della distanza è minore di cutoff calcolo il contributo dell'interazione tra la particella stessa e una ad essa vicina 
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}


void Measure(){ //Properties measurement
  int bin = 0;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream ValoriIstantanei;

//reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){ //sto calcolando l'energia potenziale 
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] );
     dy = Pbc( yold[i] - yold[j] );
     dz = Pbc( zold[i] - zold[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

//update of the histogram of g(r) 
   bin=dr/bin_size; // trovo numero bin corrispondente
   walker[igofr+bin]+=2; //aggiungo 2


     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    walker[iv] = v; //Potential energy
    walker[ik] = t; //Kinetic energy 
    walker[it] = (2.0 / 3.0) * t; //Temperature
    walker[ie] = (t+v); //Total energy

  ValoriIstantanei.open("instant_etot.out",ios::app);
  ValoriIstantanei << walker[ie]/(double)npart << endl; 
  ValoriIstantanei.close(); 

  ValoriIstantanei.open("instant_ekin.out",ios::app);
  ValoriIstantanei <<  walker[ik]/(double)npart << endl; 
  ValoriIstantanei.close(); 

  ValoriIstantanei.open("instant_epot.out",ios::app);
  ValoriIstantanei <<  walker[iv]/(double)npart << endl; 
  ValoriIstantanei.close(); 

  ValoriIstantanei.open("instant_temp.out",ios::app);
  ValoriIstantanei <<  walker[it]/(double)npart << endl; 
  ValoriIstantanei.close();
 
    return;
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
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{

   double _r, DV;
   ofstream Gofr, Gave, Ekin, Epot, Temp, Etot;
   const int wd=12;
    

    cout << "Block number " << iblk << endl;
    
    Etot.open("ave_etot.out",ios::app);
    Ekin.open("ave_ekin.out",ios::app);
    Epot.open("ave_epot.out",ios::app);
    Temp.open("ave_temp.out",ios::app);
    Gofr.open("output.gofr.0",ios::app);
    Gave.open("output.gave.0",ios::app);
    
    stima_etot = blk_av[ie]/blk_norm/(double)npart; //Total energy per particle
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

    stima_ekin = blk_av[ik]/blk_norm/(double)npart; //Kinetic energy per particle
    glob_av[ik] += stima_ekin;
    glob_av2[ik] += stima_ekin*stima_ekin;
    err_ekin=Error(glob_av[ik],glob_av2[ik],iblk);

    stima_epot = blk_av[iv]/blk_norm/(double)npart; //Potential energy per particle
    glob_av[iv] += stima_epot;
    glob_av2[iv] += stima_epot*stima_epot;
    err_epot=Error(glob_av[iv],glob_av2[iv],iblk);

    stima_temp = blk_av[it]/blk_norm/double(npart); //Temperature per particle
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);
    
    Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;

    Ekin << setw(wd) << iblk <<  setw(wd) << stima_ekin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_ekin << endl;

    Epot << setw(wd) << iblk <<  setw(wd) << stima_epot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_epot << endl;

    Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;

    cout << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;

   Gofr << setw(wd) << iblk;
    for(int i = igofr ; i < igofr + nbins; i++){
      _r = (i - 4)*bin_size;
      DV = (4./3.)*M_PI*(pow(_r+bin_size,3)-pow(_r,3));
      stima_gofr = blk_av[i]/blk_norm/(rho*DV*(double)npart); //Pressure
      glob_av[i] += stima_gofr;
      glob_av2[i] += stima_gofr*stima_gofr;
      err_gofr[i]=Error(glob_av[i],glob_av2[i],iblk);
      Gofr << setw(wd) << stima_gofr;
    }
    Gofr << endl;

//g(r) per blocco finale per ogni bin
    if(iblk == nblk){ 
     for(int i = igofr; i < igofr + nbins; i++){
     float _r = (i-4)*bin_size;
     Gave << setw(wd) << i - 4 << setw(wd) << _r+bin_size*0.5 << setw(wd) << glob_av[i]/(double)iblk << setw(wd) << err_gofr[i] << endl;
     }
    }

    cout << "----------------------------" << endl << endl;

    Etot.close();
    Ekin.close();
    Epot.close();
    Temp.close();
    Gofr.close();
    Gave.close();
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl; //sono scritte in unità di misura del lato e non di LJ in modo che possa utilizzarlo indep da densità del sistema
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz"); //così  si scrive nella cartella frames
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

double * ReadFile(const char* in_filename, int dime){

 ifstream ReadFile;

  double *valore = new double[dime];

  ReadFile.open(in_filename);
  if (ReadFile.is_open()){
    for (int i=0; i < dime; ++i){
    ReadFile >> valore[i];
    }
  } else cerr << "PROBLEM: Unable to open file" << endl;
  ReadFile.close();

 return valore;

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
