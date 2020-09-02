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
  Path(Random *rnd, int which);
  void SetCities(int cities) { n_city = cities; }
  void RandomPath();
  void Check();
  void outL(int which, double temp);
  int GetCities() { return n_city; }
  double * GetX() { return x; }
  double * GetY() { return y; }
  double ComputeL(); 
  double GetL() { return L; }
  void Mutazione(double temp);
  void Shift();
  void Permut();
  void Inversion();
  void Ordina ();
  void BestPath(int which);

};

#endif

