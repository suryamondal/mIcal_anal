#ifndef INOTRACKFITALG_H
#define INOTRACKFITALG_H
//AlgFitTrackCam

#include <vector>
using std::vector;
#include "TVector3.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "vect_manager.h"
#include "DetectorParameterDef.hh"
#include "MultiSimAnalysisDigi.hh"
#include "FieldPropagator.hh"
using namespace CLHEP;
const unsigned int doubleLa=500;
const unsigned int shiftLa=250;
class InoTrackCand;

//typedef struct {
//  InoHit* csh;
//}HitStruct;

typedef struct {
  InoCluster* csh;
}ClustStruct;


typedef struct{
  bool   Straight;
  double XPos;
  double YPos;
  double ZPos;
  int    PlaneView;
  double XPosErrSq;
  double YPosErrSq;
  int    numInList;
  double cltime;
}TrkDataStruct;

typedef struct{
  double x_k0;
  double x_k1;
  double x_k2;
  double x_k3;
  double x_k4;
  int    x_k5; //withcls;
  double C_k[5][5];
  bool   x_k6; //true for cluster in linear part, false : bended back
}FiltDataStruct;


class InoTrackFitAlg 
{
  
 public:
  InoTrackFitAlg();
  virtual ~InoTrackFitAlg();

  virtual void RunAlg( );// AlgConfig &ac, CandHandle &ch, CandContext &cx);

  //  void InitialFramework(); //const CandSliceHandle* slice, CandContext &cx);
  void InitialFramework_new();
  //  void RunTheFitter(); // CandFitTrackCamHandle &cth);
  void RunTheFitter_new();

  void StoreFilteredData(const int NewPlane);
  void StoreFilteredData_sr(const int NewPlane, double*, bool);
  void FillGapsInTrack();
  //  bool FindTheStrips(CandFitTrackCamHandle &cth, bool MakeTheTrack);
  //  bool FindTheStrips(bool MakeTheTrack);
  //  void GetFitData(int Plane1, int Plane2);
  void GetFitData_new(int& Plane1, int& Plane2);

  void ShowerStrips();
  void RemoveTrkHitsInShw();
  void ShowerSwim();

  //  void GoBackwards();
  //  void GoForwards();

  void GoBackwards_new(const bool first);
  void GoForwards_new(const bool first);

  void StraightLineFit(vector<Hep3Vector>& finder, vector<int>& loczlay, int eplnstrt, double* locslope, double* locintercpt, double* locchi2, int* locnhit, double* locexpecpos);

  void simple_track_fit(vector<Hep3Vector>& finder, int elpnsimple, double& curve, double& radii, double* x0, double* locchisq, double* xavg, int& locnhits, double& locexpecpos);
  bool GetCombiPropagator(const int Plane, const int NewPlane, const bool GoForward, double* ddS, double* ddRange);

  bool GetCombiPropagator_new(const int Plane, const int NewPlane, const bool GoForward);

  void SetBsmrGaussParameter(double bsmr){BsmrGaussPar = bsmr;};
  double GetBsmrGaussParameter() {return BsmrGaussPar;};
  //  bool Swim(double* StateVector, double* Output, const int Plane, const int NewPlane, const bool GoForward); 

  bool Swim(double* StateVector, double* Output, const int Plane, const int NewPlane, const bool GoForward, double* dS=0, double* Range=0, double* dE=0);


  bool Swim(double* StateVector, double* Output, const double zbeg, const int NewPlane, const bool GoForward, double* dS=0, double* Range=0, double* dE=0);
  bool Swim(double* StateVector, double* Output, const int Plane, const double zend, const bool GoForward, double* dS=0, double* Range=0, double* dE=0);


  bool Swim_new(double* StateVector, double* Output, const int Plane,  int& NewPlane, const bool GoForward, double* dS=0, double* Range=0, double* dE=0);


  void GetInitialCovarianceMatrix(const bool FirstIteration);
  void GetNoiseMatrix(const int Plane, const int NewPlane);
  void ExtrapCovMatrix();
  void CalcKalmanGain(const int NewPlane);
  void UpdateStateVector(const int Plane, const int NewPlane, const bool GoForward);

  int CheckFCPCUpOrDn(double *ax_k, bool DirExtraPol, int MaxMinPlane, bool GoDir);


  void UpdateStateVector_new(const int Plane, const int NewPlane, double* Output, const bool GoForward);



  
  void UpdateCovMatrix();
  void MoveArrays();
  void CheckValues(double* Input, const int NewPlane);

  void SetTrackProperties( double* Input); //CandFitTrackCamHandle &cth);  
  /*
  void SetPropertiesFromFinderTrack(InotTrackCand &cth);
  void SpectrometerSwim(CandFitTrackCamHandle &cth);
  */
  void SetRangeAnddS(); //CandFitTrackCamHandle& cth);
  void TimingFit(); //CandFitTrackCamHandle &cth);

  //  bool NDPlaneIsActive(int plane, float u, float v, float projErr);
  virtual void Trace(const char *c) const;

  void ResetCovarianceMatrix();
  //  double NDStripBegTime(CandStripHandle* Strip, double U=0, double V=0);


  void SetT( );
  void CalculateTrace(){;};
  
  bool DirectionFromFinderHits(InoTrack *trk );
  bool DirectionFromFitterHits(InoTrackCand *trk, int epln, double& xslope, double& xintercept, double& xexp);

 private:
  //  vector<StripStruct> SlcStripData[490];
  //  vector<StripStruct> InitTrkStripData[490];

  //  vector<HitStruct> SlcHitData[doubleLa]; //GMA put a very large value, but need to be put from database 
  //  vector<HitStruct> InitTrkHitData[doubleLa];

  vector<ClustStruct> SlcClustData[doubleLa]; //GMA put a very large value, but need to be put from database 
  vector<ClustStruct> InitTrkClustData[doubleLa]; //Only Finder track cluster

  vector<TrkDataStruct> TrkClustsData[doubleLa]; //TrkHitsData[doubleLa]; // TrkStripData[150];
  vector<FiltDataStruct> FilteredData[doubleLa];
  vector<FiltDataStruct> ExtraPolData[doubleLa];

  double ZPosLayer[doubleLa];
  int FCorPC; //=0;
  int FCorPCForward; // FCorPCUp; //=0;
  int FCorPCBackward; //FCorPCDn; // = 0;
  double BsmrGaussPar;

  Int_t nbfield;
  Double_t bave;
  Double_t t_bave;
  Bool_t EndofRange;
  Int_t EndofRangePlane;
  Bool_t LastIteration;
  double x_k4_biased;
  int UseGeoSwimmer;
  bool iternew;
  double x_k[6]; // x_k5 for tagging used hit or not 5];
  double predicted_k[6]; // x_k5 for tagging used hit or not 5];
  double x_k_minus[6]; // 5];
  double C_k[5][5];
  double C_k_minus[5][5];
  double C_k_intermediate[5][5];
  double F_k[5][5];
  double F_k_minus[5][5];
  double Q_k[5][5];
  double Q_k_minus[5][5];
  double K_k[5][2];
  double newfit_x_k[6];
  int H_k[2][5];
  int Identity[5][5];

  double VtxCov[5];
  double EndCov[5];

  int MaxPlane;
  int MinPlane;

  double DeltaZ;
  double DeltaPlane;
  
  //  VldContext* vldc;
  InoTrackCand* fTrackCand;

  bool timecheck;
  bool debug_fit;
  bool debug_mical;
  bool debug_fcpc;
  bool debug_new;
  int TrkFitterDebug;
  bool ZIncreasesWithTime;
  bool PassTrack;

  bool SaveData;
  bool SwimThroughShower;

  int ShowerEntryPlane;
  
  int NIter;
  int TotalNSwimFail;

  int NumFinderStrips;
  double MeanTrackTime;
  //  PlaneOutline fPL;

  double StripListTime;

  InoFittedTrack_Manager* inoFittedTrack_pointer; 
  InoTrackCand_Manager* inoTrackCand_pointer; 
  MultiSimAnalysisDigi *pAnalysis; 
  FieldPropagator *pFieldMap;
  const InoTrack* fFinderTrack;

  double StripXWidth;
  double StripYWidth;
  double LayerThickness;
  unsigned int    nLayer;

  TVector3 shiftvector; //Shift in position, while track bend back to the same layer  
  double ShiftInX;
  double ShiftInY;
  double ShiftInZ;

  TGeoManager* icalGeometry;

};
    
#endif   // ALGFITTRACKCAM_H 
