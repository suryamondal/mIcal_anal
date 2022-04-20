#ifndef MULTISIMDIGI_H
#define MULTISIMDIGI_H 1
#include <vector>
using std::vector;

#include "ParameterMessenger.hh"
#include "DetectorParameterDef.hh"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TTree.h"
#include "TFile.h"
// #include "evetree.h"
#include "RPCEve.h"
#include "globals.hh"
#include "Hits.h"
#include "HitPos.h"
#include "TProfile.h"
#include <iostream>
#include <fstream>
using namespace std;

struct vectGr{
  float x;
  float y;
  float z;
  float dx;
  float dy;
  float dz;
};

class MultiSimAnalysisDigi {
public:
  MultiSimAnalysisDigi();
  ~MultiSimAnalysisDigi();
  static MultiSimAnalysisDigi* AnPointer;

  void OpenInputRootFiles(char* inf);
  void OpenOutputRootFiles(char* outf);
  void OpenCollatedRootFile();
  void CloseInputRootFiles();
  void CloseOutputRootFiles();
  void SaveGenVisFile();

  void SetCorrTimeError(G4double val);
  void SetUnCorrTimeError(G4double val);
  void SetTimeToDigiConvVal(G4double val);
  void SetSignalSpeedVal(G4double val);
  double GetCorrTimeError() {return CorrTimeError;}
  double GetUnCorrTimeError() {return UnCorrTimeError;}
  double GetTimeToDigiConvVal() {return TimeToDigiConv;}
  double GetSignalSpeedVal() {return SignalSpeed;}

  int isInOut;
  int isVisOut;
  int isXtermOut;
  int collatedIn;

  // Collated Histograms
  TH2D* inefficiency_corx[20];
  TH2D* inefficiency_uncx[20];
  TH2D* inefficiency_uncy[20];
  TH2D* triggereffi_xevt[20];
  TH2D* triggereffi_yevt[20];
  TH2D* strp_xmulsim_cor[20];
  TH2D* strp_ymulsim_cor[20];
  TH2D* block_xmulsim[20][16][16];
  TH2D* block_ymulsim[20][16][16];

  TH1D* hdifftime1[20];
  TH1D* hdifftime2[20];
  TH1D* hxtime_ext[20];
  TH1D* hytime_ext[20];

  TH1D* hxpos_ext[20];
  TH1D* hypos_ext[20];
  TH1D* hxpos_ext_kalman[20];
  TH1D* hypos_ext_kalman[20];
  TH1D* h_hit_time_ext[20];

  TH1D* xtdc_minus_ref[20][8];
  TH1D* ytdc_minus_ref[20][8];

  TH1D* tshift_xtdc_minus_ref[20][8];
  TH1D* tshift_ytdc_minus_ref[20][8];

  Hits *H;
  HitPos *Hp;
  int EveCnt;
  int nloops;

  TFile *pRootFile;
  TFile *inputRootFile;  
  TFile *pVisFile;
  TFile* collatedRootFile;

  TTree *pEventTree;
  TTree *inputEventTree;
  TTree *visTree;
  RPCEve *data_event;
  // evetree *data_event;
  TH1D* DGap;
  TH1D* ShwXw;
  TH1D* ShwYw;
  TH1D* trk_gap;
  TH2D* trk_edge;
  TH1D* pPosX;
  TH1D* pPosY;  
  TH1D* pPosZ;
  TH2D* pPosXX;
  TH2D* pPosYY;  
  TH2D* pPosZZ;

  
  static const int nhistmx=1000;
  int   ihist;
  TH3F* gens_list[6][nhistmx]; // Not used
  vector<vectGr> gens_vect[6];
  
  // Common input read and output store for ICALsim and miniICALsim
  static const unsigned int ngenmx=50;
  UInt_t          irun;
  UInt_t          ievt;
  UInt_t          ievt2;
  UInt_t          ievt3;
  UInt_t          ngent;
  Int_t           pidin[ngenmx];   //[ngent]
  Float_t         ievt_wt;
  Int_t           intxn_id;
  Float_t         momin[ngenmx];   //[ngent]
  Float_t         thein[ngenmx];   //[ngent]
  Float_t         phiin[ngenmx];   //[ngent]
  Float_t         posxin[ngenmx];   //[ngent]
  Float_t         posyin[ngenmx];   //[ngent]
  Float_t         poszin[ngenmx];   //[ngent]

  // ICALsim Root Files Input Data Read (SIM)
  static const unsigned int nsimhtmx=4000;
  UInt_t          nsimht;
  UInt_t          detid[nsimhtmx];   //[nsimht]
  Int_t           simpdgid[nsimhtmx];   //[nsimht]
  Float_t         simtime[nsimhtmx];   //[nsimht]
  Float_t         simenr[nsimhtmx];   //[nsimht]
  Float_t         simvx[nsimhtmx];   //[nsimht]
  Float_t         simvy[nsimhtmx];   //[nsimht]
  Float_t         simvz[nsimhtmx];   //[nsimht]
  Float_t         simpx[nsimhtmx];   //[nsimht]
  Float_t         simpy[nsimhtmx];   //[nsimht]
  Float_t         simpz[nsimhtmx];   //[nsimht]
  Float_t         simlocvx[nsimhtmx];   //[nsimht]
  Float_t         simlocvy[nsimhtmx];   //[nsimht]
    
  // ICALsim Root Files Input Data Read
  static const unsigned int ndigihtmx=5000;
  UInt_t          ndigiht;
  Int_t           trigx;
  Int_t           trigy;
  UInt_t          stripid[ndigihtmx];   //[ndigiht]
  Int_t           digipdgid[ndigihtmx];   //[ndigiht]
  Int_t           digitime[ndigihtmx];   //[ndigiht]
  Int_t           digitruetime[ndigihtmx];   //[ndigiht]
  Float_t         digienr[ndigihtmx];   //[ndigiht]
  Float_t         digivx[ndigihtmx];   //[ndigiht]
  Float_t         digivy[ndigihtmx];   //[ndigiht]
  Float_t         digivz[ndigihtmx];   //[ndigiht]
  Float_t         digipx[ndigihtmx];   //[ndigiht]
  Float_t         digipy[ndigihtmx];   //[ndigiht]
  Float_t         digipz[ndigihtmx];   //[ndigiht]
  int diginoise[ndigihtmx];

  // Reco Output store for both sim.
  static const unsigned int nvishtmx=5000;
  G4float fitposxx[nvishtmx];
  G4float fitposyy[nvishtmx];
  G4float fitposzz[nvishtmx];
  G4float fitlayzz[nvishtmx];
  G4float fitlayx2[nvishtmx];
  G4float fitlayx3[nvishtmx];
  G4float fitlayx4[nvishtmx];
  G4float fitlaymom[nvishtmx];
  G4float fitlaythe[nvishtmx];
  G4float fitlayphi[nvishtmx];
  G4float extrapolxx[nvishtmx];
  G4float extrapolyy[nvishtmx];
  G4float extrapolmom[nvishtmx];

  G4float momdiff1;
  G4float radialdiff1;

  unsigned int nvisht;
  G4float clstposxx[nvishtmx];
  G4float clstposyy[nvishtmx];
  G4float clstposzz[nvishtmx];
  G4int clstposzpln[nvishtmx];
  unsigned int nvisclst;

  static  const unsigned int nthtmx=100;

  Int_t ntrecord1x;
  Int_t ntrecord1y;
  Int_t ntrecord2x;
  Int_t ntrecord2y;
  Int_t striprec1x[nthtmx];
  Int_t striprec1y[nthtmx];
  Int_t striprec2x[nthtmx];
  Int_t striprec2y[nthtmx];
  Float_t tdcrec1x[nthtmx];
  Float_t tdcrec1y[nthtmx];
  Float_t tdcrec2x[nthtmx];
  Float_t tdcrec2y[nthtmx];

  Int_t nhits_last;
  Int_t nhits_last_m1;

  Int_t strtnhitsx;
  Int_t strtnhitsy;
  Float_t strtchisqx;
  Float_t strtchisqy;
  Float_t strtintercptx;
  Float_t strtintercpty;
  Float_t strtslopex;
  Float_t strtslopey;
  Int_t strtxexpec;
  Int_t strtyexpec;

  Float_t simpleradii;
  Float_t simplecurv;
  Float_t simplex0;
  Float_t simplez0;
  Float_t simplechisqpos;
  Float_t simplechisqneg;
  Float_t simplechisqcndn;
  Float_t simpleavgxpos;
  Float_t simpleavgxneg;
  Float_t simpleavgxcndn;
  Float_t simpleavgxmeas;
  Int_t simplenhits;
  Int_t simplexexpec;
  
  Int_t ntdc1x;
  Int_t ntstrp1x;
  Int_t tdcID1x[nthtmx];
  Int_t StrpID1x[nthtmx];
  Float_t TDCval1x[nthtmx];

  Int_t ntdc2x;
  Int_t ntstrp2x;
  Int_t tdcID2x[nthtmx];
  Int_t StrpID2x[nthtmx];
  Float_t TDCval2x[nthtmx];

  Int_t ntdc1y;
  Int_t ntstrp1y;
  Int_t tdcID1y[nthtmx];
  Int_t StrpID1y[nthtmx];
  Float_t TDCval1y[nthtmx];

  Int_t ntdc2y;
  Int_t ntstrp2y;
  Int_t tdcID2y[nthtmx];
  Int_t StrpID2y[nthtmx];
  Float_t TDCval2y[nthtmx];

  Int_t nhits_below;
  Float_t ftime_last;

  static  const unsigned int ntrkmx=20;
  // Int_t           hw_trigx;
  UInt_t          ngenerated;
  UInt_t          naperture;
  UInt_t          triggeracceptance;
  Int_t           hw_trig;
  Int_t           sw_trigx;
  Int_t           sw_trigy;
  UInt_t          ntrkt;
  Int_t           itype[ntrkmx];   //[ntrkt]
  Int_t           nLayer;
  Int_t           nhits[ntrkmx];   //[ntrkt]
  Int_t           nhits_finder[ntrkmx];   //[ntrkt]
  Float_t         chisq[ntrkmx];   //[ntrkt]
  Float_t         cvalue[ntrkmx];   //[ntrkt]
  Int_t           fc_or_pc[ntrkmx];   //[ntrkt]
  Float_t         trkmm[ntrkmx];   //[ntrkt]
  Float_t         trkth[ntrkmx];   //[ntrkt]
  Float_t         trkph[ntrkmx];   //[ntrkt]
  Float_t         momvx[ntrkmx];   //[ntrkt]
  Float_t         thevx[ntrkmx];   //[ntrkt]
  Float_t         phivx[ntrkmx];   //[ntrkt]
  Float_t         posxvx[ntrkmx];   //[ntrkt]
  Float_t         posyvx[ntrkmx];   //[ntrkt]
  Float_t         poszvx[ntrkmx];   //[ntrkt]
  Float_t         momend[ntrkmx];   //[ntrkt]
  Float_t         theend[ntrkmx];   //[ntrkt]
  Float_t         phiend[ntrkmx];   //[ntrkt]
  Float_t         posxend[ntrkmx];   //[ntrkt]
  Float_t         posyend[ntrkmx];   //[ntrkt]
  Int_t            strpxend[ntrkmx];   //[ntrkt]
  Int_t            strpyend[ntrkmx];   //[ntrkt]
  Float_t         poszend[ntrkmx];   //[ntrkt]
  Float_t         tx_end[ntrkmx];   //[ntrkt]
  Float_t         ty_end[ntrkmx];   //[ntrkt]
  Float_t         momds[ntrkmx];   //[ntrkt]
  Float_t         momrg[ntrkmx];   //[ntrkt]
  Float_t         mcxgnvx[ntrkmx];   //[ntrkt]
  Float_t         mcygnvx[ntrkmx];   //[ntrkt]
  Float_t         momgnvx[ntrkmx];   //[ntrkt]
  Float_t         thegnvx[ntrkmx];   //[ntrkt]
  Float_t         phignvx[ntrkmx];   //[ntrkt]
  Float_t         momgnend[ntrkmx];   //[ntrkt]
  Float_t         thegnend[ntrkmx];   //[ntrkt]
  Float_t         phignend[ntrkmx];   //[ntrkt]
  Int_t           vtxzplane[ntrkmx];   //[ntrkt]
  Int_t           endzplane[ntrkmx];   //[ntrkt]
  Int_t           ntrkcl[ntrkmx];   //[ntrkt]
  Int_t           ntrkst[ntrkmx];   //[ntrkt]
  Int_t           ntotcl;
  Int_t           ntotst;
  Int_t           inohits;
  Int_t           orighits;
  Int_t           inoclust;
  Int_t           origclust;
  Float_t         hPathlength;
  Int_t           x_hits;
  Int_t           y_hits;
  Int_t           inohits_old;
  Int_t           orighits_old;
  Int_t           x_hits_old;
  Int_t           y_hits_old;
  Int_t           hit_wo_ghst;
  Float_t         e_hadron;
  Int_t           nhits_largest_cluster;
  Int_t           orighits_trape;
  Int_t           orighits_cluster;
  Int_t           hit_wogh_orighits;
  Float_t         theta_hadron_shw;
  Float_t         had_eigen_val[3];
  Float_t         phi_hadron_shw;
  Float_t         theta_hadron_in;
  Float_t         phi_hadron_in;
  Float_t         dot_angle_had_shw;
  Int_t           nhits_largest_cluster_selected;
  Float_t         range;
  Float_t         tx[ntrkmx];
  Float_t         ty[ntrkmx];
  Float_t         xxin[ntrkmx];
  Float_t         yyin[ntrkmx];
  Float_t         txin[ntrkmx];
  Float_t         tyin[ntrkmx];
  Float_t         therr[ntrkmx];
  Float_t         pherr[ntrkmx];
  Float_t         xxerr[ntrkmx];
  Float_t         yyerr[ntrkmx];
  Float_t         txerr[ntrkmx];
  Float_t         tyerr[ntrkmx];
  Float_t         qperr[ntrkmx];
  Float_t         xxenderr[ntrkmx];
  Float_t         yyenderr[ntrkmx];
  Float_t         txenderr[ntrkmx];
  Float_t         tyenderr[ntrkmx];
  Float_t         qpenderr[ntrkmx];
  Int_t nmxhit;

private:
  DetectorParameterDef *paradef;
  ParameterMessenger *CardFile;

  double CorrTimeError;
  double UnCorrTimeError;
  double TimeToDigiConv;
  double SignalSpeed;

  int numberInLA;
  double parirlay[3];
  double parlay[3];

  int nUp;
  int nDown;
  int DetectorType;

};

#endif


// exp([0]*x+[1])+[2]
