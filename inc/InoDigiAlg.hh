#ifndef InoDigiAlg_h
#define InoDigiAlg_h 1

#include "DetectorParameterDef.hh"
#include "MultiSimAnalysisDigi.hh"
#include "GeneralRecoInfo.hh"
#include "vect_manager.h"
#include "InoStrip.h"
#include "globals.hh"
#include "TH2.h"
#include "TH2D.h"
#include "TRandom.h"
#include "TMath.h"
using namespace std;

class InoDigiAlg {
public:
  InoDigiAlg();
  ~InoDigiAlg();

  void ReadEvent(int ixt);
  void DigitiseSimData();
  void NoiseGenLoop();
  void CalculateTrigger();
  void SaveDigiData();
  int GetRandomXY(double& GapX, TH2D* tmphistx);
  
  void SetCorrTimeSmear(G4double val) {TimeCorrSmr=val;};
  void SetUnCorrTimeSmear(G4double val) {TimeUnCorrSmr=val;};
  void SetCorrInefficiency(G4double val) {CorrIneffiPar=val;};
  void SetUnCorrXInefficiency(G4double val) {UnCorrXIneffiPar=val;};
  void SetUnCorrYInefficiency(G4double val) {UnCorrYIneffiPar=val;};
  void SetTimeToDigiConv(G4double val) {TimeToDigiConv=val;};
  void SetSignalSpeed(G4double val) {SignalSpeed=val;};
  
private:
  DetectorParameterDef* paradef;
  MultiSimAnalysisDigi* pAnalysis;
  InoStripX_Manager* inoStripX_pointer;
  InoStripY_Manager* inoStripY_pointer;
  InoCal0Hit_Manager* inocal0hit_pointer;
  double LayerThickness;
  int nLayer;
  int nIRLayer;
  int numberInMO;
  int numberInCH;
  int numberInX;
  int numberInY;
  double parino[3];
  double parlay[3];
  double parirlay[3];
  double gapino;
  double parmod[3];
  double ShiftInX;
  double ShiftInY;
  double ShiftInZ;
  double parchm[3];
  double pargas[3];
  double Xstrwd;
  double Ystrwd;
  int DetectorType;
  double ZLayerPos[300];

  int TrgLayer[10]; // = {6,7,8,9};
  int ntriglay; // = 4;
  int TrgDataX[10]; // = {0,0,0,0,0,0,0,0,0,0};
  int TrgDataY[10]; // = {0,0,0,0,0,0,0,0,0,0};
  int trigStoreX; // = 0;
  int trigStoreY; // = 0;

  int NewMultiplicity;
  double TimeCorrSmr;
  double TimeUnCorrSmr;
  double CorrIneffiPar;
  double UnCorrXIneffiPar;
  double UnCorrYIneffiPar;
  double TimeToDigiConv;
  double SignalSpeed;

  
};

#endif
