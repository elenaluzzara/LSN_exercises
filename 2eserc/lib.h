#ifndef _Pat_h_
#define _Pat_h_

class BlockingMethod {

private:

  double sum;
  double *sum_prog;
  double *su2_prog;
  double *error_prog;
  double *ave;
  double *av2;
  int k;
  int nblk; //numero di block
  int nstep; //numero di step per blocchi

public:
  // constructors
  BlockingMethod();
  BlockingMethod(int N, int L);
  // destructor
  ~BlockingMethod();
  // methods
  void Set_nblk(int N) { nblk = N;}
  void Set_nstep(int L) { nstep = L;}
  double * Get_SumProg() { return sum_prog; }
  double * Get_ErrorProg() { return error_prog; }
  void BlkMethod(double *casual);
  void Stampa(const char* nomefile);
  void SetAve(double *av) { ave = av; }
  void BlkMethod_pigreco(const char* nomefile);
  double ErrBlk(int iblk);
};

#endif 

void data_chi2(const char* nomefile,  int x[], double dati[], int dime);

void data_histo(const char* nomefile,  int x[], double dati0[], double dati1[], double dati2[], double dati3[], int dime);
