#include "InoRecoAlg.hh"

InoRecoAlg::InoRecoAlg(int isInOut) {
  grecoi = GeneralRecoInfo::GnPointer;
  pAnalysis = MultiSimAnalysisDigi::AnPointer;
  paradef = DetectorParameterDef::AnPointer;
  DetectorType = paradef->GetDetectorType();
  isData = (isInOut==2) ? 1 : 0;
  nLayer = paradef->GetnLayer();
  nIRLayer = paradef->GetnIRLayer();
  Xstrwd = paradef->GetXStrwd();
  Ystrwd = paradef->GetYStrwd();
  gapino = paradef->GetGapino();
  ShiftInX = paradef->GetRPCShift(0) + paradef->GetStackShift(0);
  ShiftInY = paradef->GetRPCShift(1) + paradef->GetStackShift(1);
  ShiftInZ = paradef->GetStackShift(2);
  // cout<<"ShiftInZ "<<ShiftInZ<<endl;
  for(int ij=0; ij<3;ij++) {
    parino[ij] = paradef->GetParino(ij);
    parlay[ij] = paradef->GetParlay(ij);
    parmod[ij] = paradef->GetParmod(ij);
    pargas[ij] = paradef->GetPargas(ij);
    parchm[ij] = paradef->GetParchm(ij);
  }
  for(int jk=0; jk<nLayer; jk++) {
    ZLayerPos[jk] = paradef->GetRPCLayerPosZ(jk) + paradef->GetShiftInZ(jk) + ShiftInZ;
  }
  numberInMO = paradef->GetnModule();
  numberInCH = paradef->GetnChamber();
  numberInX = paradef->GetnXStrip();
  numberInY = paradef->GetnYStrip();
  debug_hits=false;
  // LayerThickness = paradef->GetLayerThickness()*(1/m);

  // cout<<"numberInMO "<<numberInMO<<endl;
  // cout<<"numberInCH "<<numberInCH<<endl;
  // cout<<"numberInX "<<numberInX<<endl;
  // cout<<"numberInY "<<numberInY<<endl;
  // cout<<"Shft "<<ShiftInX<<","<<ShiftInY<<","<<ShiftInZ<<endl;
  // cout <<"parino "<<parino[0]<<"*mm, "<<parino[1]<<"*mm, "<<parino[2]<<"*mm"<<endl;
  // cout <<"parlay "<<parlay[0]<<"*mm, "<<parlay[1]<<"*mm, "<<parlay[2]<<"*mm"<<endl;
  // cout <<"parmod "<<parmod[0]<<"*mm, "<<parmod[1]<<"*mm, "<<parmod[2]<<"*mm"<<endl;
  // cout <<"parchm "<<parchm[0]<<"*mm, "<<parchm[1]<<"*mm, "<<parchm[2]<<"*mm"<<endl;
  // cout <<"pargas "<<pargas[0]<<"*mm, "<<pargas[1]<<"*mm, "<<pargas[2]<<"*mm"<<endl;
  TimeToDigiConv = pAnalysis->GetTimeToDigiConvVal();
  // cout<<"TimeToDigiConv = pAnalysis->GetTimeToDigiConvVal() "<<TimeToDigiConv<<" "<<pAnalysis->GetTimeToDigiConvVal()<<endl;
  TrgLayer[0] = 6;
  TrgLayer[1] = 7;
  TrgLayer[2] = 8;
  TrgLayer[3] = 9;
  triglays = 4;
  for(int ij=0; ij<nLayer; ij++) {
    nTriggerX[ij] = 0; nTriggerY[ij]=0;
  }
}

InoRecoAlg::~InoRecoAlg() {
  if(inoStripX_pointer) {inoStripX_pointer->InoStripX_list.clear(); delete inoStripX_pointer; inoStripX_pointer=0;}
  if(inoStripX_pointer) {inoStripX_pointer->InoStripX_list.clear(); delete inoStripX_pointer; inoStripX_pointer=0;}
  if(inoTDCx_pointer) {
    for(int ij=0; ij<nLayer; ij++) {
      for(int jk=0; jk<8; jk++) {
	inoTDCx_pointer->xtdctiming[ij][jk].clear();
      }
    }
    delete inoTDCx_pointer; inoTDCx_pointer=0;
  }
  if(inoTDCy_pointer) {
    for(int ij=0; ij<nLayer; ij++) {
      for(int jk=0; jk<8; jk++) {
	inoTDCy_pointer->ytdctiming[ij][jk].clear();
      }
    }
    delete inoTDCy_pointer; inoTDCy_pointer=0;
  }
  // cout << "Deleting InoRecoAlg ..." << endl;
}

void InoRecoAlg::ReadEvent(int ixt) {

  int layXmult[200] = {0};
  int layYmult[200] = {0};
  
  const double posoffset[2][12]={0.0};
  inoStripX_pointer = new InoStripX_Manager();
  inoStripY_pointer = new InoStripY_Manager();
  inoTDCx_pointer = new InoTDCHitx_Manager();
  inoTDCy_pointer = new InoTDCHity_Manager();
  pAnalysis->inputRootFile->cd();
  pAnalysis->inputEventTree->GetEntry(ixt);


  if(isData) {

    Int_t tmpdigipdgid=13;
    Float_t tmpdigienr=1;
    Float_t tmpdigivx=1; 
    Float_t tmpdigivy=1; 
    Float_t tmpdigivz=1; 
    Float_t tmpdigipx=1; 
    Float_t tmpdigipy=1; 
    Float_t tmpdigipz=1; 
    Int_t tmpdiginoise=0; 
    int nInCH = 0;
    int nInMO = 0;
    int nInDT = 0;


    vector <int> xptsall[nLayerMx],yptsall[nLayerMx];   //number of hits after noise rejection for position fit
    vector <int> xptsalltdc[nLayerMx][nTDCpLayer],yptsalltdc[nLayerMx][nTDCpLayer]; // raw : for one hit, return strip number otherwise -10-multiplicity
    
    for(int jk=0;jk<nLayer;jk++) {
      for (int kl=0; kl<nTDCpLayer; kl++) {
	xptsalltdc[jk][kl].clear();  yptsalltdc[jk][kl].clear();
      }
      xptsall[jk].clear(); yptsall[jk].clear();
      //printf("layer\t");
      for(int kl=0; kl<numberInX; kl++) {
	if(pAnalysis->data_event->xLayer[jk]->TestBitNumber(kl)) {
	  xptsall[jk].push_back(kl);
	  xptsalltdc[jk][kl%8].push_back(kl);
	}
      }
      for(int kl=0; kl<numberInY; kl++) {
	if(pAnalysis->data_event->yLayer[jk]->TestBitNumber(kl)) {
	  yptsall[jk].push_back(kl);
	  yptsalltdc[jk][kl%8].push_back(kl);
	}
      }
      // cout<<"ixt "<<ixt<<" "<< jk<<" "<<xptsall[jk].size()<<" "<<yptsall[jk].size()<<endl;
    } // for(int jk=0;jk<nLayer;jk++) 
    
    for(int jk=0; jk<nLayer; jk++) {
      for(int kl=0; kl<nTDCpLayer; kl++) {
	if(xptsalltdc[jk][kl].size()) {
	  unsigned tpsize = pAnalysis->data_event->vxtdc_l[jk][kl]->size();
	  for(unsigned lm=0; lm<tpsize; lm++) {
	    int stimex = 2550 + int(pAnalysis->data_event->vxtdc_l[jk][kl]->at(lm) - pAnalysis->data_event->tdc_ref_l[jk]);
	    inoTDCx_pointer->xtdctiming[jk][kl].push_back(stimex);
	    // pAnalysis->tshift_xtdc_minus_ref[jk][kl]->Fill(stimex);
	    // pAnalysis->xtdc_minus_ref[jk][kl]->Fill(int(pAnalysis->data_event->vxtdc_l[jk][kl]->at(lm) - pAnalysis->data_event->tdc_ref_l[jk]));
	  }
	}

	if(yptsalltdc[jk][kl].size()) {
	  unsigned tpsize = pAnalysis->data_event->vytdc_l[jk][kl]->size();
	  for(unsigned lm=0; lm<tpsize; lm++) {
	    int stimey = 2550 + int(pAnalysis->data_event->vytdc_l[jk][kl]->at(lm) - pAnalysis->data_event->tdc_ref_l[jk]);
	    inoTDCy_pointer->ytdctiming[jk][kl].push_back(stimey);
	    // pAnalysis->tshift_ytdc_minus_ref[jk][kl]->Fill(stimey);
	    // pAnalysis->ytdc_minus_ref[jk][kl]->Fill(int(pAnalysis->data_event->vytdc_l[jk][kl]->at(lm) - pAnalysis->data_event->tdc_ref_l[jk]));
	  }
	}
      } // for(int kl=0; kl<nTDCpLayer; kl++) {
    } // for(int jk=0; jk<nLayer; jk++) {
    

    for(int iz=0; iz<nLayer; iz++) {
      for (unsigned ix=0; ix<xptsall[iz].size(); ix++) {
	int xstripid = 0;
	xstripid<<=2;
	xstripid +=nInDT;
	xstripid<<=8;
	xstripid +=iz;
	xstripid<<=3;
	xstripid +=nInMO;
	xstripid<<=3;
	xstripid +=nInCH;
	xstripid<<=7;
	xstripid +=xptsall[iz][ix];
	xstripid<<=5;
	xstripid +=0;
	xstripid<<=3;
	xstripid +=TMath::Min(3,7);
	// cout<<"ix "<<ix<<" "<<iz<<" "<<xptsall[iz][ix]<<endl;
	InoStrip*  Xstrip = new InoStrip(); //VALGRIND
	Xstrip->SetId(xstripid);
	Xstrip->SetpdgId(tmpdigipdgid);
	int itdc = xptsall[iz][ix]%8;
	if(inoTDCx_pointer->xtdctiming[iz][itdc].size()) {
	  Xstrip->SetSmrTime(inoTDCx_pointer->xtdctiming[iz][itdc][0]);
	} else {
	  Xstrip->SetSmrTime(-10000000);
	}
	Xstrip->SetStripNumLoc(xptsall[iz][ix]);
	Xstrip->SetTrueTime(iz);
	Xstrip->SetPulse(tmpdigienr);
	G4ThreeVector tmp3v(tmpdigipx, tmpdigipy, tmpdigipz);
	Xstrip->SetMomentum(tmp3v.mag());
	Xstrip->SetTheta(tmp3v.theta());
	Xstrip->SetPhi(tmp3v.phi());
	Xstrip->SetfNoise(tmpdiginoise);
	Xstrip->SetGenPosX(tmpdigivx);
	Xstrip->SetGenPosY(tmpdigivy);
	Xstrip->SetGenPosZ(tmpdigivz);
	inoStripX_pointer->InoStripX_list.push_back(Xstrip);
      } // for (int ix=0; ix<xptsall[iz].size(); ix++) {

      for (unsigned iy=0; iy<yptsall[iz].size(); iy++) {
	int ystripid = 1;
	ystripid<<=2;
	ystripid +=nInDT;
	ystripid<<=8;
	ystripid +=iz;
	ystripid<<=3;
	ystripid +=nInMO;
	ystripid<<=3;
	ystripid +=nInCH;
	ystripid<<=7;
	ystripid +=yptsall[iz][iy];
	ystripid<<=5;
	ystripid +=0;
	ystripid<<=3;
	ystripid +=TMath::Min(3,7);
	
	InoStrip*  Ystrip = new InoStrip(); //VALGRIND
	Ystrip->SetId(ystripid);
	Ystrip->SetpdgId(tmpdigipdgid);
	int itdc = yptsall[iz][iy]%8;
	if(inoTDCy_pointer->ytdctiming[iz][itdc].size()) {
	  Ystrip->SetSmrTime(inoTDCy_pointer->ytdctiming[iz][itdc][0]);
	} else {
	  Ystrip->SetSmrTime(-10000000);
	}
	Ystrip->SetTrueTime(iz);
	Ystrip->SetStripNumLoc(yptsall[iz][iy]);
	Ystrip->SetPulse(tmpdigienr);
	G4ThreeVector tmp3v(tmpdigipx, tmpdigipy, tmpdigipz);
	Ystrip->SetMomentum(tmp3v.mag());
	Ystrip->SetTheta(tmp3v.theta());
	Ystrip->SetPhi(tmp3v.phi());
	Ystrip->SetfNoise(tmpdiginoise);
	Ystrip->SetGenPosX(tmpdigivx);
	Ystrip->SetGenPosY(tmpdigivy);
	Ystrip->SetGenPosZ(tmpdigivz);
	inoStripY_pointer->InoStripY_list.push_back(Ystrip);
      } // for (int ix=0; ix<xptsall[iz].size(); ix++) {

    } // for(int iz=0; iz<nLayer; iz++) {

  } else {
    // cout<<"-->Processing Evt Nb. "<<ixt<<" "<<pAnalysis->ievt<<" has "<<pAnalysis->ndigiht<<" hits."<<endl;
    for(unsigned ij=0;ij<pAnalysis->ndigiht;ij++) {
      unsigned istrp = pAnalysis->stripid[ij];
      unsigned strpxx = istrp;
      strpxx>>=8;
      int nnInXX = strpxx%128;
      strpxx>>=7;
      int iRPCMod = strpxx;
      int nInCH = strpxx%8; //nInCH;
      strpxx>>=3;
      int nInMO = strpxx%8; //nInMO;
      strpxx>>=3;
      int nInLA = strpxx%256; //nInLA;
      strpxx >>=8;
      int nInDT = strpxx%4; //nInDT;
      // cout<<"Hi... "<<ij<<"nInLA "<<nInLA<<" "<<endl;

      InoStrip*  Xstrip = new InoStrip(); //VALGRIND
      Xstrip->SetId(istrp);
      Xstrip->SetpdgId(pAnalysis->digipdgid[ij]);
      Xstrip->SetSmrTime(pAnalysis->digitime[ij]);
      Xstrip->SetTrueTime(pAnalysis->digitruetime[ij]);
      Xstrip->SetPulse(pAnalysis->digienr[ij]);
      
      G4ThreeVector tmp3v(pAnalysis->digipx[ij], pAnalysis->digipy[ij], pAnalysis->digipz[ij]);
      Xstrip->SetMomentum(tmp3v.mag());
      Xstrip->SetTheta(tmp3v.theta());
      Xstrip->SetPhi(tmp3v.phi());
      
      Xstrip->SetfNoise(pAnalysis->diginoise[ij]);
      Xstrip->SetGenPosX(pAnalysis->digivx[ij]);
      Xstrip->SetGenPosY(pAnalysis->digivy[ij]);
      Xstrip->SetGenPosZ(pAnalysis->digivz[ij]);
      
      if ((istrp>>31)==0) { //Most significant bit is X/Y
	inoStripX_pointer->InoStripX_list.push_back(Xstrip);
	// cout<<"nInX "<<nInDT<<" "<<nInLA<<" "<<nInMO<<" "<<nInCH<<" "<<nnInXX<<endl;
      } else {
	inoStripY_pointer->InoStripY_list.push_back(Xstrip);
	// cout<<"nInY "<<nInDT<<" "<<nInLA<<" "<<nInMO<<" "<<nInCH<<" "<<nnInXX<<endl;
      }
    }
  }
  pAnalysis->hw_trig = -1;
  
  // if((pAnalysis->trigx==4) || (pAnalysis->trigy==4)) {pAnalysis->hw_trig = 1;}
  
  // cout<<"Digi Input Reading Complete ... Reorganising the strips."<<endl;
  // cout<<"rearranging X strips"<<inoStripX_pointer->InoStripX_list.size()<<endl;
  for (unsigned ij=0; ij<inoStripX_pointer->InoStripX_list.size(); ij++) {
    unsigned istrp = inoStripX_pointer->InoStripX_list[ij]->GetId();
    
    double energy = istrp%8;
    istrp >>=8;
    int nInX = istrp%128;
    istrp>>=7;
    int iRPCMod = istrp%65536; //  2**16
    int nInCH = istrp%8;
    istrp>>=3;
    int nInMO = istrp%8;
    istrp>>=3;
    int nInLA = istrp%256;
      
    istrp>>=8;
    int nInDT = istrp%4;      
    istrp>>=2;
    // cout<<"nInX 22  "<<nInDT<<" "<<nInLA<<" "<<nInMO<<" "<<nInCH<<" "<<nInX<<endl;
    inoStripX_pointer->InoStripX_list[ij]->SetPlaneView(istrp);
    inoStripX_pointer->InoStripX_list[ij]->SetPlane(nInLA);
    inoStripX_pointer->InoStripX_list[ij]->SetRPCmod(iRPCMod);
    inoStripX_pointer->InoStripX_list[ij]->SetStrip(numberInX*numberInMO*nInDT+numberInX*nInMO+nInX);

    double xpos, ypos, zpos;
    if(DetectorType==0) { // ICAL
      xpos = (1/m)*( (nInDT-1)*(2*parino[0]+gapino) - parlay[0]  + (2*nInMO+1)*parmod[0] -pargas[0] + Xstrwd*(nInX+0.5) + ShiftInX);
      ypos = (1/m)*(- parmod[1]  + (2*nInCH+1)*parchm[1] + ShiftInY);
      zpos = ZLayerPos[nInLA]/m;
    } else { // mICAL
      xpos = (1/m)*(-pargas[0] + Xstrwd*(nInX+0.5) + ShiftInX);
      ypos = (1/m)*(-parlay[1] + ShiftInY);
      zpos = ZLayerPos[nInLA]/m;
    }
    //(1/m)*(-(nLayer-1)*(parirlay[2]+parlay[2])+(nInLA)*2*(parirlay[2] + parlay[2]) + ShiftInZ); //AAR:** changes for Central Iron Layer **
      
    //if smearing //RANDOM
    inoStripX_pointer->InoStripX_list[ij]->SetXYPos(xpos);
    inoStripX_pointer->InoStripX_list[ij]->SetZPos(zpos);
    if(TrgLayer[0]==nInLA || TrgLayer[1]==nInLA || TrgLayer[2]==nInLA || TrgLayer[3]==nInLA) {
      nTriggerX[nInLA]++;
    }
    // cout<<"11 "<<ij<<" "<<xpos<<" "<<zpos<<" "<<inoStripX_pointer->InoStripX_list[ij]->GetXYPos() << " " << inoStripX_pointer->InoStripX_list[ij]->GetZPos() << " " << inoStripX_pointer->InoStripX_list[ij]->GetGenPosX() << " "<< inoStripX_pointer->InoStripX_list[ij]->GetGenPosZ()<<endl;
    // cout<<"pAnalysis "<<pAnalysis<<" "<<pAnalysis->pPosX<<endl;
    pAnalysis->pPosX->Fill(100*xpos - 0.1*inoStripX_pointer->InoStripX_list[ij]->GetGenPosX());
    pAnalysis->pPosZ->Fill(100*zpos - 0.1*inoStripX_pointer->InoStripX_list[ij]->GetGenPosZ());
    pAnalysis->pPosXX->Fill(100*xpos, 100*xpos - 0.1*inoStripX_pointer->InoStripX_list[ij]->GetGenPosX());
    pAnalysis->pPosZZ->Fill(100*zpos, 100*zpos - 0.1*inoStripX_pointer->InoStripX_list[ij]->GetGenPosZ());
    grecoi->xlayer_occu[nInLA]->Fill(nInX);
    layXmult[nInLA]++;
    StripBankX[nInLA].push_back(inoStripX_pointer->InoStripX_list[ij]);
    // cout<<"pAnalysis->(pPosX,pPosZ,pPosXX,pPosZZ)->Fill();"<<endl;
    if (energy >100000 || abs(ypos)>100000) cout <<"ypos "<<ypos<<" "<<energy<<endl;
  }
  
  // cout<<"rearranging Y strips"<<inoStripY_pointer->InoStripY_list.size()<<endl;
  for (unsigned ij=0; ij<inoStripY_pointer->InoStripY_list.size(); ij++) {
    unsigned istrp = inoStripY_pointer->InoStripY_list[ij]->GetId();
    double energy = istrp%8;
    istrp >>=8;
    int nInY = istrp%128;
    istrp>>=7;
    int iRPCMod = istrp%65536; //  2**16
    int nInCH = istrp%8;
    istrp>>=3;
    int nInMO = istrp%8;
    istrp>>=3;
    int nInLA = istrp%256;
    istrp>>=8;
    int nInDT = istrp%4;      
    istrp>>=2;
      
    inoStripY_pointer->InoStripY_list[ij]->SetPlaneView(istrp);
    inoStripY_pointer->InoStripY_list[ij]->SetPlane(nInLA);
    inoStripY_pointer->InoStripY_list[ij]->SetRPCmod(iRPCMod);
    inoStripY_pointer->InoStripY_list[ij]->SetStrip(numberInY*nInCH+nInY);
    double xpos, ypos, zpos;
    if(DetectorType==0) { // ICAL
      xpos = (1/m)*( (nInDT-1)*(2*parino[0]+gapino) - parlay[0]  + (2*nInMO+1)*parmod[0] + ShiftInX);
      ypos = (1/m)*(-parmod[1] +(2*nInCH+1)*parchm[1] -pargas[1] + Ystrwd*(nInY+0.5) + ShiftInY);
      zpos = ZLayerPos[nInLA]/m;
    } else { // mICAL
      xpos = (1/m)*(-parlay[0] + ShiftInX);
      ypos = (1/m)*(-pargas[1] + Ystrwd*(nInY+0.5) + ShiftInY);
      zpos = ZLayerPos[nInLA]/m;
    }


    // (1/m)*(-(nLayer-1)*(parirlay[2]+parlay[2])+(nInLA)*2*(parirlay[2] + parlay[2]) + ShiftInZ); //AAR:** changes for Central Iron Layer **
    // double zpos = inoStripY_pointer->InoStripY_list[ij]->GetGenPosZ()/m;
    inoStripY_pointer->InoStripY_list[ij]->SetXYPos(ypos);
    inoStripY_pointer->InoStripY_list[ij]->SetZPos(zpos);

    // cout<<"22 "<<ij<<" "<<ypos<<" "<<zpos<<" "<<inoStripY_pointer->InoStripY_list[ij]->GetXYPos() << " " << inoStripY_pointer->InoStripY_list[ij]->GetZPos() << " " << inoStripY_pointer->InoStripY_list[ij]->GetGenPosY() << " "<< inoStripY_pointer->InoStripY_list[ij]->GetGenPosZ()<<endl;

    if(TrgLayer[0]==nInLA || TrgLayer[1]==nInLA || TrgLayer[2]==nInLA || TrgLayer[3]==nInLA) {
      nTriggerY[nInLA]++;
    }
    
    pAnalysis->pPosY->Fill(100*ypos - 0.1*inoStripY_pointer->InoStripY_list[ij]->GetGenPosY());
    pAnalysis->pPosZ->Fill(100*zpos - 0.1*inoStripY_pointer->InoStripY_list[ij]->GetGenPosZ());
    pAnalysis->pPosYY->Fill(100*ypos, 100*ypos - 0.1*inoStripY_pointer->InoStripY_list[ij]->GetGenPosY());
    pAnalysis->pPosZZ->Fill(100*zpos, 100*zpos - 0.1*inoStripY_pointer->InoStripY_list[ij]->GetGenPosZ());

    grecoi->ylayer_occu[nInLA]->Fill(nInY);
    layYmult[nInLA]++;
    StripBankY[nInLA].push_back(inoStripY_pointer->InoStripY_list[ij]);
    // cout<<"pAnalysis->(pPosY,pPosZ,pPosYY,pPosZZ)->Fill();"<<endl;
    // cout<<"Ical0cal0SD::EndOfEvent(G4HCofThisEvent*) {....."<<endl;
    if (energy >100000 || abs(xpos)>100000) cout <<"xpos "<<xpos<<" "<<energy<<endl;
  }
  
  // cout<<"Strip rearrangement complete..."<<endl;
  // cout<<"hhh "<<pAnalysis->momin[0]<<" "<<cos(pAnalysis->thein[0])<<" "<<pAnalysis->phiin[0]*180/3.14<<endl;
  
  for(int il=0; il<nLayer; il++) {
    grecoi->xlayer_mult[il]->Fill(layXmult[il]);
    grecoi->ylayer_mult[il]->Fill(layYmult[il]);
  }
  
  // cout<<"Strip Rearrangement complete."<<endl;
}


void InoRecoAlg::PerformTrackReconstruction() {

  if(debug_hits) {
    cout<<"StripBank X"<<endl;
    for(int ij=0; ij<nLayer; ij++) {
      for(unsigned jk=0; jk<StripBankX[ij].size(); jk++) {
	cout<< "InoHits():" 
	    <<std::setw(4) <<jk <<" "
	    << " pln="   <<std::setw(4)<< StripBankX[ij][jk]->GetPlane() 
	    << " view="   <<std::setw(4)<< StripBankX[ij][jk]->GetPlaneView()
	    << " ID=" <<std::setw(14)<< StripBankX[ij][jk]->GetId()
	    << " Strp=" <<std::setw(5)<<StripBankX[ij][jk]->GetStripNumLoc()
	    << " time="<<std::setw(9)<<StripBankX[ij][jk]->GetSmrTime()
	    <<endl;
      }
    }
    cout<<"StripBank Y"<<endl;
    for(int ij=0; ij<nLayer; ij++) {
      for(unsigned jk=0; jk<StripBankY[ij].size(); jk++) {
	cout<< "InoHits():" 
	    <<std::setw(4) <<jk <<" "
	    << " pln="   <<std::setw(4)<< StripBankY[ij][jk]->GetPlane() 
	    << " view="   <<std::setw(4)<< StripBankY[ij][jk]->GetPlaneView()
	    << " ID=" <<std::setw(14)<< StripBankY[ij][jk]->GetId()
	    << " Strp=" <<std::setw(5)<<StripBankY[ij][jk]->GetStripNumLoc()
	    << " time="<<std::setw(9)<<StripBankY[ij][jk]->GetSmrTime()
	    <<endl;
      }
    }
  }

  pAnalysis->strtnhitsx=-1;
  pAnalysis->strtnhitsy=-1;
  pAnalysis->strtxexpec=-10;
  pAnalysis->strtyexpec=-10;
  pAnalysis->strtchisqx=-10.0;
  pAnalysis->strtchisqy=-10.0;
  pAnalysis->strtintercptx=-10.0;
  pAnalysis->strtintercpty=-10.0;
  pAnalysis->strtslopex=-10.0;
  pAnalysis->strtslopey=-10.0;

  pAnalysis->simpleradii=-100000;
  pAnalysis->simplecurv=-100000;
  pAnalysis->simplex0=-100000;
  pAnalysis->simplez0=-100000;
  pAnalysis->simplechisqpos=-100000;
  pAnalysis->simplechisqneg=-100000;
  pAnalysis->simplechisqcndn=-100000;
  pAnalysis->simpleavgxpos=-100000;
  pAnalysis->simpleavgxneg=-100000;
  pAnalysis->simpleavgxcndn=-100000;
  pAnalysis->simplenhits = -1;
  pAnalysis->simplexexpec = -10;
  pAnalysis->momdiff1 = -10.0;
  pAnalysis->radialdiff1 = -10.0;
  
  pAnalysis->ntdc1x=0; pAnalysis->ntstrp1x=0;
  pAnalysis->ntdc2x=0; pAnalysis->ntstrp2x=0;
  pAnalysis->ntdc1y=0; pAnalysis->ntstrp1y=0;
  pAnalysis->ntdc2y=0; pAnalysis->ntstrp2y=0;
  pAnalysis->ntrecord1x=0;
  pAnalysis->ntrecord1y=0;
  pAnalysis->ntrecord2x=0;
  pAnalysis->ntrecord2y=0;
  pAnalysis->nvisclst=0;
  pAnalysis->ntrkt=0;
  pAnalysis->nvisht = 0;
  pAnalysis->sw_trigx = 0;
  pAnalysis->sw_trigy = 0;
  for(int tr1=0; tr1<nLayer; tr1++) {
    for(int tr2=0; tr2<triglays; tr2++) {
      if(tr1==TrgLayer[tr2]) {
	if(nTriggerX[tr1]>0) {pAnalysis->sw_trigx++;}
	if(nTriggerY[tr1]>0) {pAnalysis->sw_trigy++;}
      }
    }
  }
  
  if(1) {//pAnalysis->sw_trigx>3 && pAnalysis->sw_trigy>3) {
    InoTrackFinder trackfinder;
    if(debug_hits) cout<<"Running the finder now ...."<<endl;
    trackfinder.RunTheFinder();
    if(debug_hits) cout<<"Track finding completed..."<<endl;
    
    InoTrack_Manager *pinotrack = InoTrack_Manager::APointer;
    if (pinotrack) {
      if(debug_hits)
	cout<<"Going to Track Fitter..."<<endl;
      InoTrackFitAlg trackfitter;
      trackfitter.RunAlg();
      InoTrackCand_Manager *pfitTrack = InoTrackCand_Manager::APointer;
      if(debug_hits)
	cout<<"fitted track manager "<<pfitTrack<<endl;
      if (pfitTrack) {
	if(debug_hits) cout<<"fitted track list size "<<" "<<pfitTrack->InoTrackCand_list.size()<<endl;
	unsigned ij=0;
	if(pfitTrack->InoTrackCand_list.size()) {
	  for (unsigned jk=0; jk<pfitTrack->InoTrackCand_list.size() ; jk++) {
	    if (ij <pAnalysis->ntrkmx) {
	      SaveRecoDataToRootFiles(ij,pfitTrack->InoTrackCand_list[jk]);
	      // SaveRecoDataToRootFiles;
	      // cout<<"input "<<pAnalysis->momin[ij]<<" "<<pAnalysis->thein[ij]<<" "<< pAnalysis->phiin[ij]<<endl;
	      // cout<<"momval "<<pAnalysis->trkmm[ij]<<" "<<pAnalysis->trkth[ij]<<" "<< pAnalysis->trkph[ij]<<" "<<pAnalysis->endzplane[ij]<<" "<<pAnalysis->nhits_last<<" "<<pAnalysis->nhits_last_m1<<endl;
	    }
	    ij++;
	  } //for (unsigned jk=0; jk<pfitTrack->InoTrackCand_list.size() ; jk++)
	  pAnalysis->ntrkt = ij;
	} else {
	  if(debug_hits) cout <<"XXXXXXXXXXXXXXXXXXXXXXX No tracks in InoTrackCand_Manager::APointer "<<endl;
	}
      } else { // else of if (pfitTrack) {
	if(debug_hits) cout <<"XXXXXXXXXXXXXXXXXXXXXXX InoTrackCand_Manager::APointer is not found "<<endl;
      } // else if (pfitTrack) {
    } else {
      if(debug_hits) cout <<"XXXXXXXXXXXXXXXXXXXXXXX InoTrack_Manager::APointer is not found "<<endl;
    }
  } else {
    if(debug_hits) cout<<"x-x-x-x-x-x-x-x-x- No hit in trigger layer -x-x-x-x-x-x-x-x"<<endl;
  }
  // cout<<"ntrkt = "<<pAnalysis->ntrkt<<" "<<pAnalysis->nvisclst<<" "<<pAnalysis->nvisht<<endl;
  pAnalysis->pEventTree->Fill();
}

void InoRecoAlg::PerformHadronReconstruction() {
  cout<<"No hadron analysis code...."<<endl;
}

void InoRecoAlg::SaveRecoDataToRootFiles(int iTrackNum, InoTrackCand* pfittedTrack) {
  pAnalysis->trkcmom[iTrackNum] =  pfittedTrack->GetCircleMom();
  pAnalysis->trkcs[iTrackNum] =  pfittedTrack->GetCircleChisq();
  pAnalysis->itype[iTrackNum] =  pfittedTrack->GetFitType();
  pAnalysis->nhits[iTrackNum] = (pfittedTrack->GetNDOF()+5)/2;
  pAnalysis->chisq[iTrackNum] =  pfittedTrack->GetChi2();
  pAnalysis->cvalue[iTrackNum] = pfittedTrack->Getcval();
  // pAnalysis->track_merge_flg[iTrackNum] = pfittedTrack->GetTrackMergeVar();
  pAnalysis->fc_or_pc[iTrackNum] = pfittedTrack->GetFCPC();
  // if(pAnalysis->fc_or_pc[iTrackNum]<0 || pAnalysis->fc_or_pc[iTrackNum]>255) {
  //   cout<<"Error Error Error Error Error....333 "<<ij<<" "<<pAnalysis->fc_or_pc[iTrackNum]<<" "<<pfittedTrack->GetFCPC()<<endl;
  // }
  // pAnalysis->recomomnew[iTrackNum] = pfittedTrack->GetNewMomentum();
  // pAnalysis->newrecoflag[iTrackNum] = pfittedTrack->GetNewFitOut();
  // cout<<"Track # --- "<<iTrackNum<<endl;
  for(unsigned int kpc=0; kpc<pfittedTrack->ClustsInTrack.size(); kpc++) {
    // cout<<"save "<<pfittedTrack->ClustsInTrack[kpc]->GetXPos()<<" "<<pfittedTrack->ClustsInTrack[kpc]->GetYPos()<<" "<<pfittedTrack->ClustsInTrack[kpc]->GetZPos()<<" "<<pfittedTrack->ClustsInTrack[kpc]->GetTime()<<endl;
  }

  vector<int> totlay;
  if(iTrackNum==0) {
    pAnalysis->nvisht = pfittedTrack->zfitpos1.size();
    for(unsigned int kpx=0; kpx<pAnalysis->nvisht; kpx++) {
      pAnalysis->fitposxx[kpx] = pfittedTrack->xfitpos1[kpx];
      pAnalysis->fitposyy[kpx] = pfittedTrack->yfitpos1[kpx];
      pAnalysis->fitposzz[kpx] = pfittedTrack->zfitpos1[kpx]; 
      pAnalysis->fitlayzz[kpx] = pfittedTrack->zfitlay1[kpx]; 
      pAnalysis->fitlayx2[kpx] = pfittedTrack->filteredx2[kpx]; 
      pAnalysis->fitlayx3[kpx] = pfittedTrack->filteredx3[kpx]; 
      pAnalysis->fitlayx4[kpx] = pfittedTrack->filteredx4[kpx]; 
      pAnalysis->fitlaymom[kpx] = pfittedTrack->filteredmom[kpx];
      pAnalysis->fitlaythe[kpx] = pfittedTrack->filteredthe[kpx];
      pAnalysis->fitlayphi[kpx] = pfittedTrack->filteredphi[kpx];
      pAnalysis->extrapolxx[kpx] = pfittedTrack->extrapolx0[kpx]; 
      pAnalysis->extrapolyy[kpx] = pfittedTrack->extrapolx1[kpx]; 
      pAnalysis->extrapolmom[kpx] = pfittedTrack->extrapolmom[kpx]; 
    }

    pAnalysis->momdiff1 = 0.0;
    for(unsigned kpe=0; kpe<pfittedTrack->momvecdiff1.size(); kpe++) {
      pAnalysis->momdiff1 += pfittedTrack->momvecdiff1[kpe];
    }
    
    pAnalysis->radialdiff1 = 0.0;
    for(unsigned kpe=0; kpe<pfittedTrack->radialdiff1.size(); kpe++) {
      pAnalysis->radialdiff1 += pfittedTrack->radialdiff1[kpe];
    }

    pAnalysis->nvisclst = pfittedTrack->ClustsInTrack.size();
    for(unsigned int kpc=0; kpc<pAnalysis->nvisclst; kpc++) {
      pAnalysis->clstposxx[kpc] = pfittedTrack->ClustsInTrack[kpc]->GetXPos();
      pAnalysis->clstposyy[kpc] = pfittedTrack->ClustsInTrack[kpc]->GetYPos();
      pAnalysis->clstposzz[kpc] = pfittedTrack->ClustsInTrack[kpc]->GetZPos();
      pAnalysis->clstposzpln[kpc] = pfittedTrack->ClustsInTrack[kpc]->GetZPlane();
      if(debug_hits) {
	cout<< "Final Track Cluster "<<std::setw(3)<<kpc
	    << " pln="    <<std::setw(3)<< pfittedTrack->ClustsInTrack[kpc]->GetZPlane()
	    << " beX=" <<std::setw(4)<< pfittedTrack->ClustsInTrack[kpc]->GetBegXStrip()
	    << "  -  " <<std::setw(4)<< pfittedTrack->ClustsInTrack[kpc]->GetEndXStrip()
	    << " beY=" <<std::setw(4)<< pfittedTrack->ClustsInTrack[kpc]->GetBegYStrip()
	    << "  -  " <<std::setw(4)<< pfittedTrack->ClustsInTrack[kpc]->GetEndYStrip()
	    << " Xpos="<<std::setw(8)<< pfittedTrack->ClustsInTrack[kpc]->GetXPos()
	    << " Ypos="<<std::setw(8)<< pfittedTrack->ClustsInTrack[kpc]->GetYPos()
	    << " chg=" <<std::setw(8)<< pfittedTrack->ClustsInTrack[kpc]->GetPulse()
	    << " time="  <<std::setw(8)<< pfittedTrack->ClustsInTrack[kpc]->GetTime()
	    << " timex="  <<std::setw(8)<< pfittedTrack->ClustsInTrack[kpc]->GetBegTime()
	    << " timey="  <<std::setw(8)<< pfittedTrack->ClustsInTrack[kpc]->GetEndTime()
	    << endl;
      }
      bool fild = true;
      int plnnum = pfittedTrack->ClustsInTrack[kpc]->GetZPlane();
      grecoi->fit_nLayer->Fill(plnnum);
      for(unsigned sz=0; sz<totlay.size(); sz++) {
	if(plnnum == totlay[sz]) {
	  fild=false; break;
	}
      }
      if(fild) totlay.push_back(plnnum);
    }
    
  }
  grecoi->nlayer_fit->Fill(totlay.size());
  
  pAnalysis->trkmm[iTrackNum] = pfittedTrack->GetMomentum();
  pAnalysis->trkth[iTrackNum] = pfittedTrack->GetTheta();
  pAnalysis->trkph[iTrackNum] = pfittedTrack->GetPhi();
  // cout<<"input "<<pAnalysis->momin[iTrackNum]<<" "<<pAnalysis->thein[iTrackNum]<<" "<< pAnalysis->phiin[iTrackNum]<<endl;
  // cout<<"momval "<<pAnalysis->trkmm[iTrackNum]<<" "<<pAnalysis->trkth[iTrackNum]<<" "<< pAnalysis->trkph[iTrackNum]<<endl;
  // pAnalysis->therr[iTrackNum] = pfittedTrack->GetThErr();
  // pAnalysis->pherr[iTrackNum] = pfittedTrack->GetPhErr();
  
  pAnalysis->momvx[iTrackNum] = pfittedTrack->GetMomentumCurve();
  pAnalysis->thevx[iTrackNum] = acos(pfittedTrack->GetDirCosZ());
  pAnalysis->phivx[iTrackNum] = atan2(pfittedTrack->GetDirCosV(),
				      pfittedTrack->GetDirCosU());
  pAnalysis->posxvx[iTrackNum] = pfittedTrack->GetVtxU();
  pAnalysis->posyvx[iTrackNum] = pfittedTrack->GetVtxV();
  pAnalysis->poszvx[iTrackNum] = pfittedTrack->GetVtxZ();

  // pAnalysis->strpxvx[iTrackNum] = 
  // pAnalysis->strpyvx[iTrackNum] = pfittedTrack->GetVtxV();

  
  pAnalysis->momend[iTrackNum] = pfittedTrack->GetEndMomentumCurve();
  pAnalysis->theend[iTrackNum] = acos(pfittedTrack->GetEndDirCosZ());
  pAnalysis->phiend[iTrackNum] = atan2(pfittedTrack->GetEndDirCosV(),
				       pfittedTrack->GetEndDirCosU());
  pAnalysis->tx_end[iTrackNum] = pfittedTrack->GetEndDirCosU();
  pAnalysis->ty_end[iTrackNum] = pfittedTrack->GetEndDirCosV();
  
  pAnalysis->posxend[iTrackNum] = pfittedTrack->GetEndU();
  pAnalysis->posyend[iTrackNum] = pfittedTrack->GetEndV();
  pAnalysis->poszend[iTrackNum] = pfittedTrack->GetEndZ();

  pAnalysis->strtslopex = pfittedTrack->GetStraightLineSlopeX();
  pAnalysis->strtslopey = pfittedTrack->GetStraightLineSlopeY();
  pAnalysis->strtintercptx = pfittedTrack->GetStraightLineInterceptX();
  pAnalysis->strtintercpty = pfittedTrack->GetStraightLineInterceptY();
  pAnalysis->strtchisqx = pfittedTrack->GetStraightLineChi2X();
  pAnalysis->strtchisqy = pfittedTrack->GetStraightLineChi2Y();
  pAnalysis->strtnhitsx = pfittedTrack->GetStraightLineNhitsX();
  pAnalysis->strtnhitsy = pfittedTrack->GetStraightLineNhitsY();
  double p2x1 = (1000*pfittedTrack->GetStraightLineXExpec())-ShiftInX+pargas[0];
  int x2s1 = int(p2x1/Xstrwd);
  pAnalysis->strtxexpec = x2s1;
  double p2y1 = (1000*pfittedTrack->GetStraightLineYExpec())-ShiftInY+pargas[0];
  int y2s1 = int(p2y1/Xstrwd);
  pAnalysis->strtyexpec = y2s1;

  pAnalysis->simpleradii = pfittedTrack->GetSimpleRadii();
  pAnalysis->simplecurv = pfittedTrack->GetSimpleCurv();
  pAnalysis->simplex0 = pfittedTrack->GetSimpleX0();
  pAnalysis->simplez0 = pfittedTrack->GetSimpleZ0();
  pAnalysis->simplechisqpos = pfittedTrack->GetSimpleChi2Pos();
  pAnalysis->simplechisqneg = pfittedTrack->GetSimpleChi2Neg();
  pAnalysis->simplechisqcndn = pfittedTrack->GetSimpleChi2Cndn();
  pAnalysis->simpleavgxpos = pfittedTrack->GetSimpleAvgXPos();
  pAnalysis->simpleavgxneg = pfittedTrack->GetSimpleAvgXNeg();
  pAnalysis->simpleavgxcndn = pfittedTrack->GetSimpleAvgXCndn();
  pAnalysis->simpleavgxmeas = pfittedTrack->GetSimpleAvgXMeas();
  pAnalysis->simplenhits = pfittedTrack->GetSimpleNhits();
  double p2x2 = (1000*pfittedTrack->GetSimpleXExpec())-ShiftInX+pargas[0];
  int x2s2 = int(p2x2/Xstrwd);
  pAnalysis->simplexexpec = x2s2;

  double pos1x = (1000*pfittedTrack->GetEndU())-ShiftInX+pargas[0];
  double pos2x = (1000*pfittedTrack->GetEndV())-ShiftInY+pargas[1];
  // cout<<"planes "<<pfittedTrack->GetVtxPlane()<<" "<<pfittedTrack->GetEndPlane()<<endl;
  int xtmp1 = (pos1x/Xstrwd);
  int ytmp1 = (pos2x/Xstrwd);
  
  // cout<<"pos1x "<<pos1x<<" "<<xtmp1<<" "<<pos2x<<" "<<ytmp1<<" "<<1000*pfittedTrack->GetEndU()<<" "<<1000*pfittedTrack->GetEndV()<<" "<<ShiftInX<<" "<<ShiftInY<<" "<<pargas[0]<<" "<<pargas[1]<<endl;
  pAnalysis->strpxend[iTrackNum] = xtmp1;
  pAnalysis->strpyend[iTrackNum] = ytmp1;
    
  pAnalysis->momds[iTrackNum] = pfittedTrack->GetMomentumdS();
  pAnalysis->momrg[iTrackNum] = pfittedTrack->GetMomentumRange();
  
  pAnalysis->vtxzplane[iTrackNum] = pfittedTrack->GetVtxPlane();
  pAnalysis->endzplane[iTrackNum] = pfittedTrack->GetEndPlane();
  
  if(iTrackNum==0) {
    pAnalysis->ntdc1x=0; pAnalysis->ntstrp1x=0;
    pAnalysis->ntdc2x=0; pAnalysis->ntstrp2x=0;
    pAnalysis->ntdc1y=0; pAnalysis->ntstrp1y=0;
    pAnalysis->ntdc2y=0; pAnalysis->ntstrp2y=0;
    pAnalysis->ntrecord1x=0;
    pAnalysis->ntrecord1y=0;
    pAnalysis->ntrecord2x=0;
    pAnalysis->ntrecord2y=0;
    pAnalysis->nhits_last_m1 = pfittedTrack->GetNhitsEndPlaneM1();
    pAnalysis->nhits_last = pfittedTrack->GetNhitsEndPlane();
    
    
    int epln = pAnalysis->endzplane[iTrackNum];
    if(epln>0) {

      pAnalysis->nhits_below = pfittedTrack->HitsNotInTrack.size();
      pAnalysis->ftime_last = pfittedTrack->GetNewTimeEndPlaneExp();

      int yij=0;
      for(unsigned jh1=0; jh1<StripBankX[epln].size(); jh1++) {
	int tmpid1 = StripBankX[epln][jh1]->GetId();
	tmpid1 >>= 8;
	int strpn = tmpid1&0x07F;//tmpid1%128;
	if(abs(pAnalysis->strpxend[0]-strpn)<3) {
	  for(unsigned kl1=0; kl1<inoTDCx_pointer->xtdctiming[epln][strpn/8].size(); kl1++) {
	    if(abs(TimeToDigiConv*inoTDCx_pointer->xtdctiming[epln][strpn/8][kl1]-pAnalysis->ftime_last)>10.0) {
	      pAnalysis->striprec1x[yij] = strpn;
	      pAnalysis->tdcrec1x[yij] = TimeToDigiConv*inoTDCx_pointer->xtdctiming[epln][strpn/8][kl1];
	      yij++;
	    }
	  }
	}
      }
      pAnalysis->ntrecord1x=yij;
      yij=0;
      for(unsigned jh1=0; jh1<StripBankY[epln].size(); jh1++) {
	int tmpid1 = StripBankY[epln][jh1]->GetId();
	tmpid1 >>= 8;
	int strpn = tmpid1&0x07F;//128;
	if(abs(pAnalysis->strpyend[0]-strpn)<3) {
	  for(unsigned kl1=0; kl1<inoTDCy_pointer->ytdctiming[epln][strpn/8].size(); kl1++) {
	    if(abs(TimeToDigiConv*inoTDCy_pointer->ytdctiming[epln][strpn/8][kl1]-pAnalysis->ftime_last)>10.0) {
	      pAnalysis->striprec1y[yij] = strpn;
	      pAnalysis->tdcrec1y[yij] = TimeToDigiConv*inoTDCy_pointer->ytdctiming[epln][strpn/8][kl1];
	      yij++;
	    }
	  }
	}
      }
      pAnalysis->ntrecord1y=yij;
      yij=0;

      for(unsigned jh1=0; jh1<StripBankX[epln-1].size(); jh1++) {
	int tmpid1 = StripBankX[epln-1][jh1]->GetId();
	tmpid1 >>= 8;
	int strpn = tmpid1&0x07F;//tmpid1%128;
	if(abs(pAnalysis->strpxend[0]-strpn)<3) {
	  for(unsigned kl1=0; kl1<inoTDCx_pointer->xtdctiming[epln-1][strpn/8].size(); kl1++) {
	    if(abs(TimeToDigiConv*inoTDCx_pointer->xtdctiming[epln-1][strpn/8][kl1]-pAnalysis->ftime_last)>10.0) {
	      pAnalysis->striprec2x[yij] = strpn;
	      pAnalysis->tdcrec2x[yij] = TimeToDigiConv*inoTDCx_pointer->xtdctiming[epln-1][strpn/8][kl1];
	      yij++;
	    }
	  }
	}
      }
      pAnalysis->ntrecord2x=yij;
      yij=0;

      for(unsigned jh1=0; jh1<StripBankY[epln-1].size(); jh1++) {
	int tmpid1 = StripBankY[epln-1][jh1]->GetId();
	tmpid1 >>= 8;
	int strpn = tmpid1&0x07F;//tmpid1%128;
	if(abs(pAnalysis->strpyend[0]-strpn)<3) {
	  for(unsigned kl1=0; kl1<inoTDCy_pointer->ytdctiming[epln-1][strpn/8].size(); kl1++) {
	    if(abs(TimeToDigiConv*inoTDCy_pointer->ytdctiming[epln-1][strpn/8][kl1]-pAnalysis->ftime_last)>10.0) {
	      pAnalysis->striprec2y[yij] = strpn;
	      pAnalysis->tdcrec2y[yij] = TimeToDigiConv*inoTDCy_pointer->ytdctiming[epln-1][strpn/8][kl1];
	      yij++;
	    }
	  }
	}
      }
      pAnalysis->ntrecord2y=yij;
      yij=0;

      
      if(debug_hits) {cout<<"epln = "<<epln<<endl;
	cout<<"ss1 "<<StripBankX[epln].size()<<endl;}
      for(unsigned jh1=0; jh1<StripBankX[epln].size(); jh1++) {
	int tmpid1 = StripBankX[epln][jh1]->GetId();
	tmpid1 >>= 8;
	pAnalysis->StrpID1x[jh1] = tmpid1%128;
      }
      pAnalysis->ntstrp1x = StripBankX[epln].size();
      if(debug_hits) cout<<"ss2 "<<StripBankY[epln].size()<<endl;
      for(unsigned jh1=0; jh1<StripBankY[epln].size(); jh1++) {
	int tmpid1 = StripBankY[epln][jh1]->GetId();
	tmpid1 >>= 8;
	pAnalysis->StrpID1y[jh1] = tmpid1%128;
      }
      pAnalysis->ntstrp1y = StripBankY[epln].size();

      if(debug_hits) cout<<"ss3 "<<StripBankX[epln-1].size()<<endl;
      for(unsigned jh1=0; jh1<StripBankX[epln-1].size(); jh1++) {
	int tmpid1 = StripBankX[epln-1][jh1]->GetId();
	tmpid1 >>= 8;
	pAnalysis->StrpID2x[jh1] = tmpid1%128;
      }
      pAnalysis->ntstrp2x = StripBankX[epln-1].size();
      if(debug_hits) cout<<"ss4 "<<StripBankY[epln-1].size()<<endl;
      for(unsigned jh1=0; jh1<StripBankY[epln-1].size(); jh1++) {
	int tmpid1 = StripBankY[epln-1][jh1]->GetId();
	tmpid1 >>= 8;
	pAnalysis->StrpID2y[jh1] = tmpid1%128;
      }
      pAnalysis->ntstrp2y = StripBankY[epln-1].size();
    
      int lpct=0;
      for(int kg1=0; kg1<8; kg1++) {
	for(unsigned jh1=0; jh1<inoTDCx_pointer->xtdctiming[epln][kg1].size(); jh1++) {
	  pAnalysis->tdcID1x[lpct] = kg1;
	  pAnalysis->TDCval1x[lpct] = TimeToDigiConv*inoTDCx_pointer->xtdctiming[epln][kg1][jh1];
	  if(debug_hits) cout<<"tdc1x "<<kg1<<" "<<jh1<<" "<<TimeToDigiConv<<" "<<TimeToDigiConv*inoTDCx_pointer->xtdctiming[epln][kg1][jh1]<<" "<<pAnalysis->TDCval1x[lpct]<<endl;
	  lpct++;
	}
      }
      pAnalysis->ntdc1x = lpct;
      lpct=0;
      for(int kg1=0; kg1<8; kg1++) {
	for(unsigned jh1=0; jh1<inoTDCy_pointer->ytdctiming[epln][kg1].size(); jh1++) {
	  pAnalysis->tdcID1y[lpct] = kg1;
	  pAnalysis->TDCval1y[lpct] = TimeToDigiConv*inoTDCy_pointer->ytdctiming[epln][kg1][jh1];
	  if(debug_hits) cout<<"tdc1y "<<kg1<<" "<<jh1<<" "<<TimeToDigiConv<<" "<<TimeToDigiConv*inoTDCy_pointer->ytdctiming[epln][kg1][jh1]<<" "<<pAnalysis->TDCval1y[lpct]<<endl;
	  lpct++;
	}
      }
      pAnalysis->ntdc1y = lpct;

      lpct=0;
      for(int kg1=0; kg1<8; kg1++) {
	for(unsigned jh1=0; jh1<inoTDCx_pointer->xtdctiming[epln-1][kg1].size(); jh1++) {
	  pAnalysis->tdcID2x[lpct] = kg1;
	  pAnalysis->TDCval2x[lpct] = TimeToDigiConv*inoTDCx_pointer->xtdctiming[epln-1][kg1][jh1];
	  lpct++;
	}
      }
      pAnalysis->ntdc2x = lpct;
      lpct=0;
      for(int kg1=0; kg1<8; kg1++) {
	for(unsigned jh1=0; jh1<inoTDCy_pointer->ytdctiming[epln-1][kg1].size(); jh1++) {
	  pAnalysis->tdcID2y[lpct] = kg1;
	  pAnalysis->TDCval2y[lpct] = TimeToDigiConv*inoTDCy_pointer->ytdctiming[epln-1][kg1][jh1];
	  lpct++;
	}
      }
      pAnalysis->ntdc2y = lpct;

    }
  }
  
  pAnalysis->xxerr[iTrackNum] = pfittedTrack->GetVtxUError(); //GMA14 all these 15 tesrms
  pAnalysis->yyerr[iTrackNum] = pfittedTrack->GetVtxVError();
  pAnalysis->txerr[iTrackNum] = pfittedTrack->GetVtxdUError();
  pAnalysis->tyerr[iTrackNum] = pfittedTrack->GetVtxdVError();
  pAnalysis->qperr[iTrackNum] = pfittedTrack->GetVtxQPError();
	    
  pAnalysis->xxenderr[iTrackNum] = pfittedTrack->GetEndUError();
  pAnalysis->yyenderr[iTrackNum] = pfittedTrack->GetEndVError();
  pAnalysis->txenderr[iTrackNum] = pfittedTrack->GetEnddUError();
  pAnalysis->tyenderr[iTrackNum] = pfittedTrack->GetEnddVError();
  pAnalysis->qpenderr[iTrackNum] = pfittedTrack->GetEndQPError();
	    
  // pAnalysis->xxin[iTrackNum] = pfittedTrack->GetVtxXX();
  // pAnalysis->yyin[iTrackNum] = pfittedTrack->GetVtxYY();
  // pAnalysis->txin[iTrackNum] = pfittedTrack->GetVtxTX();
  // pAnalysis->tyin[iTrackNum] = pfittedTrack->GetVtxTY();
  // pAnalysis->tx[iTrackNum] = pfittedTrack->GetVtxdU();
  // pAnalysis->ty[iTrackNum] = pfittedTrack->GetVtxdV();
  pAnalysis->nhits_finder[iTrackNum] = pfittedTrack->GetNFinderHits();
  //Initialise with this mainly for noise track
	      
  // cout<<"---------------------------------------------------------------"<<endl;
  // cout<<"p = "<<pAnalysis->trkmm[iTrackNum]<<", momds = "<<pAnalysis->momds[iTrackNum]<<", E_mu = "<<pAnalysis->momin[0]<<endl;
  // cout<<"---------------------------------------------------------------"<<endl;

  pAnalysis->momgnvx[iTrackNum]  = pAnalysis->momgnend[iTrackNum] = 0.0;
  pAnalysis->thegnvx[iTrackNum] = pAnalysis->thegnend[iTrackNum] = 10.;
  pAnalysis->phignvx[iTrackNum] = pAnalysis->phignend[iTrackNum] = 10.;
	    
  vector<InoCluster*> tmpclusts = pfittedTrack->ClustsInTrack; // fTrack->ClustsInTrack;
  
  int plane = pfittedTrack->GetVtxPlane();
	    
  // cout <<"genplane "<<plane<<endl;
  for (unsigned kl = 0; kl <tmpclusts.size(); kl++)	{
    if (tmpclusts[kl]->GetZPlane()==plane) {
      pAnalysis->momgnvx[iTrackNum] = 0.001*tmpclusts[kl]->HitsInCluster[0]->GetMomentum();
      pAnalysis->thegnvx[iTrackNum] = tmpclusts[kl]->HitsInCluster[0]->GetTheta();
      pAnalysis->phignvx[iTrackNum] = tmpclusts[kl]->HitsInCluster[0]->GetPhi();
      break;
    } // if (tmpclusts[kl]->GetZPlane()==plane)
  } //  for (unsigned kl = 0; kl <tmpclusts.size(); kl++)
	    
  plane = pfittedTrack->GetEndPlane();
	    
  for (unsigned kl = 0; kl <tmpclusts.size(); kl++) {
    if (tmpclusts[kl]->GetZPlane()==plane) {
      pAnalysis->momgnend[iTrackNum] = 0.001*tmpclusts[kl]->HitsInCluster[0]->GetMomentum();
      pAnalysis->thegnend[iTrackNum] = tmpclusts[kl]->HitsInCluster[0]->GetTheta();
      pAnalysis->phignend[iTrackNum] = tmpclusts[kl]->HitsInCluster[0]->GetPhi();
      break;
    }
  }
  //tmpclusts.clear();
}
