#ifndef INONEWFITALG_H
#define INONEWFITALG_H
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include "vect_manager.h"
#include "DetectorParameterDef.hh"
#include "MultiSimAnalysisDigi.hh"
#include "FieldPropagator.hh"
#include "SwimSwimmer.h"
#include "SwimParticle.h"
#include <sys/time.h>
#include "G4SystemOfUnits.hh"
#include <cassert>
#include <cmath>
#include "TMath.h"
#include "TVector3.h"
#include "TMatrixD.h"
//#include "TVectorD.h"
//#include "TMatrixTBase.h"
#include "TMatrixDEigen.h"

using std::vector;
using namespace std;

const int nlayermx=10;
const int nvectormx=6;

class InoNewFitAlg {
 public:
  InoNewFitAlg();
  InoNewFitAlg(double* psvVtx, double* mpts,double* mptserr, double* mptsz, int* occu, int nmiss, int vtxp, bool TrkDir);
  bool DoIterations();
  double GetFinalStateVectorElement(int ix) {return finalStateVector[ix];};
  double GetFinalStateVectorError(int ix);
  ~InoNewFitAlg();
  bool Swim(double* StateVector, double* Output, const int Plane, int& NewPlane, const bool GoForward);
  void GetNoiseMatrix(const int Plane, const int NewPlane);

 private:
  int VtxPlane;
  double input_x_k[nvectormx];
  double input_x_k_err[nvectormx][nvectormx];
  double vtx_x_k[nvectormx];
  double vtx_x_k_err[nvectormx][nvectormx];
  double posin[2*nlayermx];
  double posinerr2[2*nlayermx];
  double posinz[nlayermx];
  int fgood;
  double datapts[2*nlayermx];
  double datavec[2*nlayermx];
  vector<int> occulyr;
  double CovMatrix[2*nlayermx][2*nlayermx];
  double AMatrix[2*nlayermx][nvectormx];
  double exppts[2*nlayermx];
  double finalStateVector[nvectormx];
  double finalSVerror[nvectormx];
  bool ZIncreasesWithTime;
  double LayerThickness;
  bool debug_fit;
  double x_k_minus[6];
  double C_k_intermediate[5][5];
  double Q_k_minus[2][2];
  double F_k_minus[5][5];
  double ZPosLayer[255];
  int nLayer;
  bool GoFrd;
  double ShiftInX;
  double ShiftInY;
  double ShiftInZ;
};
#endif
