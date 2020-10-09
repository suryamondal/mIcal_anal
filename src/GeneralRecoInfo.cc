#include "GeneralRecoInfo.hh"

GeneralRecoInfo* GeneralRecoInfo::GnPointer;
GeneralRecoInfo::GeneralRecoInfo() {
  GnPointer = this;
  paradef = DetectorParameterDef::AnPointer;
  DetectorType = paradef->GetDetectorType();
  nLayer = paradef->GetnLayer();
  nIRLayer = paradef->GetnIRLayer();
  for(int ij=0; ij<3; ij++) {
    parirlay[ij] = paradef->GetParirlay(ij);
    parlay[ij] = paradef->GetParlay(ij);
  }
  nXStrip = paradef->GetnXStrip();
  nYStrip = paradef->GetnYStrip();

  for(int ij=0; ij<nLayer; ij++) {
    for(int jk=0; jk<3; jk++) {
      align_xstr_ydev[ij][jk] = 0.0;
      align_ystr_xdev[ij][jk] = 0.0;
      align_xstr_xdev[ij][jk] = 0.0;
      align_ystr_ydev[ij][jk] = 0.0;
    }
    timeoffsetx[ij] = 0.0;//h_timeoffsetx->GetBinContent(ij+1);
    timeoffsety[ij] = 0.0;//h_timeoffsety->GetBinContent(ij+1);
    for(int jk=0; jk<nYStrip; jk++) {
      xtoffset[ij][jk] = 0.0;//h_xtoffset->GetBinContent(ij+1,jk+1);
      ytoffset[ij][jk] = 0.0;//h_ytoffset->GetBinContent(ij+1,jk+1);
      for(int kl=0; kl<nYStrip; kl++) {
	xt_slope_cor[ij][jk][kl] = 0.0;//h_xt_slope_cor->GetBinContent(ij+1,jk+1,kl+1);
	yt_slope_cor[ij][jk][kl] = 0.0;//h_yt_slope_cor->GetBinContent(ij+1,jk+1,kl+1);
      }
    }
  }

}

GeneralRecoInfo::GeneralRecoInfo(char* fileInName) {
  GnPointer = this;
  paradef = DetectorParameterDef::AnPointer;
  DetectorType = paradef->GetDetectorType();
  nLayer = paradef->GetnLayer();
  nIRLayer = paradef->GetnIRLayer();
  for(int ij=0; ij<3; ij++) {
    parirlay[ij] = paradef->GetParirlay(ij);
    parlay[ij] = paradef->GetParlay(ij);
  }
  nXStrip = paradef->GetnXStrip();
  nYStrip = paradef->GetnYStrip();

  InCorrFile = new TFile(fileInName,"read");
  
  h_align_xstr_ydev = (TH2D*)InCorrFile->Get("h_align_xstr_ydev");
  h_align_ystr_xdev = (TH2D*)InCorrFile->Get("h_align_ystr_xdev");
  h_align_xstr_xdev = (TH2D*)InCorrFile->Get("h_align_xstr_xdev");
  h_align_ystr_ydev = (TH2D*)InCorrFile->Get("h_align_ystr_ydev");
  h_strpos_vs_time = (TH3D*)InCorrFile->Get("h_strpos_vs_time");
  h_xposerrsq = (TH2D*)InCorrFile->Get("h_xposerrsq");
  h_yposerrsq = (TH2D*)InCorrFile->Get("h_yposerrsq");
  h_timeserrx2 = (TH1D*)InCorrFile->Get("h_timeserrx2");
  h_timeserry2 = (TH1D*)InCorrFile->Get("h_timeserry2");
  h_xtoffset = (TH2D*)InCorrFile->Get("h_xtoffset");
  h_ytoffset = (TH2D*)InCorrFile->Get("h_ytoffset");
  h_xtoffystr = (TH2D*)InCorrFile->Get("h_xtoffystr");
  h_ytoffxstr = (TH2D*)InCorrFile->Get("h_ytoffxstr");
  h_xt_slope_cor = (TH3D*)InCorrFile->Get("h_xt_slope_cor");
  h_yt_slope_cor = (TH3D*)InCorrFile->Get("h_yt_slope_cor");
  h_timeoffsetx = (TH1D*)InCorrFile->Get("h_timeoffsetx");
  h_timeoffsety = (TH1D*)InCorrFile->Get("h_timeoffsety");

  for(int ij=0; ij<nLayer; ij++) {
    for(int jk=0; jk<3; jk++) {
      align_xstr_ydev[ij][jk] = h_align_xstr_ydev->GetBinContent(ij+1,jk+1);
      align_ystr_xdev[ij][jk] = h_align_ystr_xdev->GetBinContent(ij+1,jk+1);
      align_xstr_xdev[ij][jk] = h_align_xstr_xdev->GetBinContent(ij+1,jk+1);
      align_ystr_ydev[ij][jk] = h_align_ystr_ydev->GetBinContent(ij+1,jk+1);
    }
    timeoffsetx[ij] = h_timeoffsetx->GetBinContent(ij+1);
    timeoffsety[ij] = h_timeoffsety->GetBinContent(ij+1);
    for(int jk=0; jk<nYStrip; jk++) {
      xtoffset[ij][jk] = h_xtoffset->GetBinContent(ij+1,jk+1);
      ytoffset[ij][jk] = h_ytoffset->GetBinContent(ij+1,jk+1);
      for(int kl=0; kl<nYStrip; kl++) {
	xt_slope_cor[ij][jk][kl] = h_xt_slope_cor->GetBinContent(ij+1,jk+1,kl+1);
	yt_slope_cor[ij][jk][kl] = h_yt_slope_cor->GetBinContent(ij+1,jk+1,kl+1);
      }
    }
  }

}

GeneralRecoInfo::~GeneralRecoInfo() {
  cout<<"Deleting GeneralRecoInfo class pointer ..."<<endl;
}

void GeneralRecoInfo::OpenRootFiles(char* outfile) {
  sprintf(outfile,"%s.root",outfile);
  GeneralFileOut = new TFile(outfile, "RECREATE"); //VALGRIND
  if (!GeneralFileOut) {
    cout << "Error opening histogram file !" << endl;
    exit(-1);
  } else {
    cout<< "Output stored in root file: "<< outfile <<endl;
  }

  int nbinxMax, nbinyMax;
  double xrmin,xrmax,yrmin,yrmax;
  if(DetectorType) {
    nbinxMax = 80; nbinyMax = 80;
  } else {
    nbinxMax = 320; nbinyMax = 320;
  }  
  xrmin = -parirlay[0]-25;
  xrmax = parirlay[0]+25;
  yrmin = -parirlay[1]-25;
  yrmax = parirlay[1]+25;
  
  xyvsbxout = new TH2D("xyvsbxout", "xyvsbxout", nbinxMax+1, xrmin, xrmax, nbinyMax, yrmin, yrmax);
  xyvsbyout = new TH2D("xyvsbyout", "xyvsbyout", nbinxMax+1, xrmin, xrmax, nbinyMax, yrmin, yrmax);
  char title[300];
  for(int ij=0; ij<nLayer; ij++) {
    sprintf(title, "xlayer_mult_l%i", ij);
    xlayer_mult[ij] = new TH1D(title, title, nXStrip+1, -0.5, nXStrip+0.5);
    
    sprintf(title, "ylayer_mult_l%i", ij);
    ylayer_mult[ij] = new TH1D(title, title, nYStrip+1, -0.5, nYStrip+0.5);

    sprintf(title, "xlayer_occu_l%i", ij);
    xlayer_occu[ij] = new TH1D(title, title, nXStrip+1, -0.5, nXStrip+0.5);

    sprintf(title, "ylayer_occu_l%i", ij);
    ylayer_occu[ij] = new TH1D(title, title, nYStrip+1, -0.5, nYStrip+0.5);
    // sprintf(title, "raw_occu_l%i", ij);
    // raw_occu[ij] = new TH2D(title, title, nXStrip+1, -0.5, nXStrip+0.5, nYStrip+1, -0.5, nYStrip+0.5);

    // sprintf(title, "xlayer_muon_occu_l%i", ij);
    // xlayer_muon_occu[ij] = new TH1D(title, title, nXStrip+1, -0.5, nXStrip+0.5);
    
    // sprintf(title, "ylayer_muon_occu_l%i", ij);
    // ylayer_muon_occu[ij] = new TH1D(title, title, nYStrip+1, -0.5, nYStrip+0.5);
  }
  
  nlayer_fit = new TH1D("nlayer_fit", "nlayer_fit", nLayer+1, -0.5, nLayer+0.5);
  fit_nLayer = new TH1D("fit_nLayer", "fit_nLayer", nLayer+1, -0.5, nLayer+0.5);
  
  cout<<"void GeneralRecoInfo::OpenRootFiles(char* outfile) { complete..."<<endl;
  
}

void GeneralRecoInfo::CloseRootFiles() {
  if (GeneralFileOut) {
    GeneralFileOut->cd();
    // if(xyvsbxin) {xyvsbxin->Write(); delete xyvsbxin; xyvsbxin=0;}
    // if(xyvsbyin) {xyvsbyin->Write(); delete xyvsbyin; xyvsbyin=0;}
    if(xyvsbxout) {xyvsbxout->Write(); delete xyvsbxout; xyvsbxout=0;}
    if(xyvsbyout) {xyvsbyout->Write(); delete xyvsbyout; xyvsbyout=0;}
    for(int ij=0; ij<nLayer; ij++) {
      if(xlayer_mult[ij]) {xlayer_mult[ij]->Write(); delete xlayer_mult[ij]; xlayer_mult[ij]=0;}
      if(ylayer_mult[ij]) {ylayer_mult[ij]->Write(); delete ylayer_mult[ij]; ylayer_mult[ij]=0;}
      if(raw_occu[ij]) {raw_occu[ij]->Write(); delete raw_occu[ij]; raw_occu[ij]=0;}
      if(xlayer_occu[ij]) {xlayer_occu[ij]->Write(); delete xlayer_occu[ij]; xlayer_occu[ij]=0;}
      if(ylayer_occu[ij]) {ylayer_occu[ij]->Write(); delete ylayer_occu[ij]; ylayer_occu[ij]=0;}
      if(xlayer_muon_occu[ij]) {xlayer_muon_occu[ij]->Write(); delete xlayer_muon_occu[ij]; xlayer_muon_occu[ij]=0;}
      if(ylayer_muon_occu[ij]) {ylayer_muon_occu[ij]->Write(); delete ylayer_muon_occu[ij]; ylayer_muon_occu[ij]=0;}
    }
    if(nlayer_fit) {nlayer_fit->Write(); delete nlayer_fit; nlayer_fit=0;}
    if(fit_nLayer) {fit_nLayer->Write(); delete fit_nLayer; fit_nLayer=0;}
      
    GeneralFileOut->Close();
    delete GeneralFileOut; GeneralFileOut=0;
    cout << "Histogram file  made !" << endl;
  }  else {
    cout << "Histogram file not made !" << endl;
  } // if (GeneralFileOut) {
}


