#ifndef _Pat_h_
#define _Pat_h_

#include "random.h"

class Path{

protected:
  int n_city;
  int *path;
  int *old_path;
  double L; 
  double * x; 
  double * y;
  Random *_rand; 
  double walker;
  double glob_av;
  double glob_av2;
  double blk_av;
  int blk_norm;

public:
  // constructors
  //Path();
  Path(Random *rnd, int which);
  // destructor
  //~Path();
  // methods
  //void SetInitPath(void);
  //void Which(int which, int p1, int p2, int seed[]);
  void SetCities(int cities) { n_city = cities; }
  void RandomPath();
  //void SetPath(int * value); //per inizializzarlo 
  //void SetRndom(Random *rnd) { _rand = rnd; }
  void Check();
  void outL(int which, double temp);
  //int * GetPath() { return path; }
  int GetCities() { return n_city; }
  double * GetX() { return x; }
  double * GetY() { return y; }
  //void SetX(double * a) { x = a; }
  //void SetY(double * b) { y = b; }
  double ComputeL(); 
  double GetL() { return L; }
  void Mutazione(double temp);
  void Shift();
  void Permut();
  void Inversion();
  void Ordina ();
  void BestPath(int which);
  double Error(double sum, double sum2, int iblk);
  void Reset(int iblk);
  void Accumulate();
  void Averages(double * Lx, int nstep, int nblk,int which); //int temp);
  int Pbc(int r);

};

#endif 
int SelectOperator(double p, int pop, Random * rnd);
void Crossover(Path mpath, Path fpath, Random * rnd);

