#ifndef GENERALRECOINFO_H
#define GENERALRECOINFO_H 1

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
#include "globals.hh"
#include "TProfile.h"
#include <iostream>
#include <fstream>
using namespace std;

class GeneralRecoInfo {
public:
  GeneralRecoInfo();
  GeneralRecoInfo(char* fileInName);
  ~GeneralRecoInfo();
  void OpenRootFiles(char* outf);
  void CloseRootFiles();
  static GeneralRecoInfo* GnPointer;
  TFile* GeneralFileOut;

  TH2D* xyvsbxin;
  TH2D* xyvsbyin;
  TH2D* xyvsbxout;
  TH2D* xyvsbyout;
  
  TH1D* xlayer_mult[200];
  TH1D* ylayer_mult[200];
  TH2D* raw_occu[200]; // not
  TH1D* xlayer_occu[200]; 
  TH1D* ylayer_occu[200];
  TH1D* xlayer_muon_occu[200]; // not
  TH1D* ylayer_muon_occu[200]; // not
  TH1D* nlayer_fit;
  TH1D* fit_nLayer;

  TFile* InCorrFile;
  TH2D*	h_align_xstr_ydev;
  TH2D*	h_align_ystr_xdev;
  TH2D*	h_align_xstr_xdev;
  TH2D*	h_align_ystr_ydev;
  TH3D*	h_strpos_vs_time;
  TH2D*	h_xposerrsq;
  TH2D*	h_yposerrsq;
  TH1D*	h_timeserrx2;
  TH1D*	h_timeserry2;
  TH1D*	h_timeoffsetx;
  TH1D*	h_timeoffsety;
  TH2D*	h_xtoffset;
  TH2D*	h_ytoffset;
  TH2D*	h_xtoffystr;
  TH2D*	h_ytoffxstr;
  TH3D*	h_xt_slope_cor;
  TH3D*	h_yt_slope_cor;


  double align_xstr_ydev[200][3];
  double align_ystr_xdev[200][3];
  double align_xstr_xdev[200][3];
  double align_ystr_ydev[200][3];
  double timeoffsetx[200];
  double timeoffsety[200];
  double xtoffset[200][128];
  double ytoffset[200][128];
  double xt_slope_cor[200][128][128];
  double yt_slope_cor[200][128][128];
private:
  DetectorParameterDef *paradef;
  int nLayer;
  int nIRLayer;
  double parirlay[3];
  double parlay[3];
  int nXStrip;
  int nYStrip;
  int DetectorType;

};
#endif

