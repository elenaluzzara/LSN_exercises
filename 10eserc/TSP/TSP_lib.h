#ifndef _Pat_h_
#define _Pat_h_

#include "random.h"

class Path{

protected:
  int n_city;
  int *path[100];
  int *old_path[100];
  int * tappe;
  double * L; 
  double * x; 
  double * y;
  Random * _rand; 
  double *walker;
  double *glob_av;
  double *glob_av2;
  double *blk_av;
  int blk_norm;

public:
  // constructors
  //Path();
  Path(Random * rnd, int rank);
  void InitialPathRank(); 
  void SetCities(int cities) { n_city = cities; }
  void RandomPath(int tap);
  void ResetPath();
  void BestSetPath(int * value); 
  void Check();
  void NuovaPop();
  void InitialPath(int which);
  int SelectOperator(int p); 
  void Crossover(int mother,int father, int ind); //int idx, int idx2,int cont, int cont2);
  int GetCities() { return n_city; }
  double * GetX() { return x; }
  double * GetY() { return y; }
  double ComputeL(int tap); 
  double * GetL() { return L; }
  void Mutazione();
  void Shift(int tap);
  void Permut(int tap);
  void Inversion(int tap);
  void Ordina ();
  void BestPath(int rank);
  void SetBestPath(int*value);
  int * ReadBestPath();
  void BestL(int rank);
  double Error(double sum, double sum2, int iblk);
  void Reset(int iblk);
  void Accumulate();
  void Averages(int iblk, int which,int nstep);
  void InitializeRank(int rank);

};

#endif 
