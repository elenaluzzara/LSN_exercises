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
  Random *_rand; 

public:
  // constructors
  Path(Random *rnd, int which);
  void SetCities(int cities) { n_city = cities; }
  void RandomPath(int tap);
  void ResetPath();
  void SetL(double *value) {L=value;} //per inizializzarlo 
  void Check();
  void NuovaPop();
  void InitialPath(int which);
  int SelectOperator(int p); 
  void Crossover(int mother,int father, int ind);
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
  void BestPath(int which);
  void BestL(int which);
  void AverageL(int which);

};

#endif 

