/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=1000;
int n_props;
int iv,ik,it,ie,igofr;
double walker[m_props];
double bin_size,nbins;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];
double vx_2[m_part],vy_2[m_part],vz_2[m_part]; //per riscalare le velocità

// averages
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double stima_epot, stima_ekin, stima_etot, stima_temp, stima_gofr,stima_gave,err_epot,err_ekin, err_etot,err_temp;
double err_gofr[100];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, nblk, iprint, seed, irestart;
double delta;

//functions
void Input(void);
void Move(void);
void OldConf(void);
void ConfFinal(void);
void ConfXYZ(int);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Error(double,double ,int);
void Recupera(void);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
