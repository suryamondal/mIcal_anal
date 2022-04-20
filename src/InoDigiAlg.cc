#include "InoDigiAlg.hh"


InoDigiAlg::InoDigiAlg() {
  pAnalysis = MultiSimAnalysisDigi::AnPointer;
  paradef = DetectorParameterDef::AnPointer;
  DetectorType = paradef->GetDetectorType();
  nLayer = paradef->GetnLayer();
  nIRLayer = paradef->GetnIRLayer();
  Xstrwd = paradef->GetXStrwd();
  Ystrwd = paradef->GetYStrwd();
  gapino = paradef->GetGapino();
  ShiftInX = paradef->GetRPCShift(0); // + paradef->GetStackShift(0);
  ShiftInY = paradef->GetRPCShift(1); // + paradef->GetStackShift(1);
  ShiftInZ = paradef->GetRPCShift(2); // + paradef->GetStackShift(2);
  for(int ij=0; ij<3;ij++) {
    parino[ij] = paradef->GetParino(ij);
    parlay[ij] = paradef->GetParlay(ij);
    parmod[ij] = paradef->GetParmod(ij);
    pargas[ij] = paradef->GetPargas(ij);
    parchm[ij] = paradef->GetParchm(ij);
  }
  for(int jk=0; jk<nLayer; jk++) {
    ZLayerPos[jk] = paradef->GetRPCLayerPosZ(jk);
  }
  numberInMO = paradef->GetnModule();
  numberInCH = paradef->GetnChamber();
  numberInX = paradef->GetnXStrip();
  numberInY = paradef->GetnYStrip();
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

  TrgLayer[0] = 6;
  TrgLayer[1] = 7;
  TrgLayer[2] = 8;
  TrgLayer[3] = 9;
  ntriglay = 4;
  for(int jk=0; jk<nLayer; jk++) {
    TrgDataX[jk] = 0; TrgDataY[jk] = 0;
  }
  trigStoreX = 0;
  trigStoreY = 0;

  SetCorrTimeSmear(0.0);
  SetUnCorrTimeSmear(0.0);
  SetCorrInefficiency(0.00);
  SetUnCorrXInefficiency(0.00);
  SetUnCorrYInefficiency(0.00);
  SetTimeToDigiConv(0.1);
  SetSignalSpeed(0.15);

  NewMultiplicity = 1;
  
  inocal0hit_pointer = new InoCal0Hit_Manager();
  inoStripX_pointer = new InoStripX_Manager();
  inoStripY_pointer = new InoStripY_Manager();
}

InoDigiAlg::~InoDigiAlg() {
  if(inoStripX_pointer) {inoStripX_pointer->InoStripX_list.clear(); delete inoStripX_pointer; inoStripX_pointer=0;}
  if(inoStripY_pointer) {inoStripY_pointer->InoStripY_list.clear(); delete inoStripY_pointer; inoStripY_pointer=0;}
  if(inocal0hit_pointer) {inocal0hit_pointer->InoCal0Hit_list.clear(); delete inocal0hit_pointer; inocal0hit_pointer=0;}
  // cout << "Deleting InoDigiAlg ..." << endl;
}

void InoDigiAlg::ReadEvent(int ixt) {
  pAnalysis->inputRootFile->cd();
  pAnalysis->inputEventTree->GetEntry(ixt);
  // cout <<"siminput "<< pAnalysis->nsimht<<endl;
  for(unsigned ij=0;ij<pAnalysis->nsimht;ij++) {
    InoCal0Hit* newHit = new InoCal0Hit();
    G4ThreeVector mom(pAnalysis->simpx[ij],pAnalysis->simpy[ij],pAnalysis->simpz[ij]);
    G4ThreeVector pos(pAnalysis->simvx[ij],pAnalysis->simvy[ij],pAnalysis->simvz[ij]);
    newHit->SetHitId(pAnalysis->detid[ij]);
    newHit->SetpdgId(pAnalysis->simpdgid[ij]);
    newHit->SetEdep( pAnalysis->simenr[ij] );
    newHit->SetTime( pAnalysis->simtime[ij] );
    newHit->SetPos(pos);
    newHit->SetMom(mom);

    newHit->SetLocalXPos(pAnalysis->simlocvx[ij]);
    newHit->SetLocalYPos(pAnalysis->simlocvy[ij]);
    // cout<<"ij "<<ij<<" "<<pos<<" "<<pAnalysis->simlocvx[ij]<<" "<<pAnalysis->simlocvy[ij]<<endl;
    inocal0hit_pointer->InoCal0Hit_list.push_back(newHit);
  }
  // cout<<"Event reading complete. "<<inocal0hit_pointer->InoCal0Hit_list.size()<<" Total sim hit stored."<<endl;
}

void InoDigiAlg::DigitiseSimData() {
  // cout<<"Digitising data..."<<endl;

  int iMnT = 10000; //Should we use these at all ? GMA151001
  int iMxT = -10000;
  double eMx = 100;
  int nHits = 0;
  const int MxStrip=3;
  // cout<<"inocal0hit_pointer-> "<<inocal0hit_pointer->InoCal0Hit_list.size()<<endl;
  for (unsigned ij=0; ij<inocal0hit_pointer->InoCal0Hit_list.size(); ij++) {
    double edep = inocal0hit_pointer->InoCal0Hit_list[ij]->GetEdep();
    eMx +=edep;
    nHits++;
    G4ThreeVector posvec2 = inocal0hit_pointer->InoCal0Hit_list[ij]->GetPos();

    int nInX[MxStrip]={-1, -1, -1};
    int nInY[MxStrip]={-1, -1, -1}; 

    unsigned long detid = inocal0hit_pointer->InoCal0Hit_list[ij]->GetHitId();
    nInY[0] = detid%128;
    detid>>=7;
    nInX[0] = detid%128;
    if (nInX[0] <0 && nInY[0] <0) { continue;}
    detid>>=7;
    int iRPCMod = detid;
    int nInCH = detid%8; //nInCH;
    detid>>=3;
    int nInMO = detid%8; //nInMO;
    detid>>=3;
    int nInLA = detid%256; //nInLA;
    detid >>=8;
    int nInDT = detid%4;
    
    // if(pAnalysis->collatedIn) {
    //   CorrIneffiPar = pAnalysis->inefficiency_corx[nInLA]->GetBinContent(nInX[0]+1,nInY[0]+1);
    // }
    // if(gRandom->Rndm(0) < CorrIneffiPar) continue;

    int pdgid = inocal0hit_pointer->InoCal0Hit_list[ij]->GetpdgId();
    double atime = inocal0hit_pointer->InoCal0Hit_list[ij]->GetTime(); 

    double CorrTimeSmr = gRandom->Gaus(0,TimeCorrSmr);
    double tmpatimeX = atime + SignalSpeed*(nInY[0]+0.5) + CorrTimeSmr; // + gRandom->Gaus(0,TimeUnCorrSmr);
    double tmpatimeY = atime + SignalSpeed*(nInX[0]+0.5) + CorrTimeSmr; // + gRandom->Gaus(0,TimeUnCorrSmr);

    int nInT = int(atime/TimeToDigiConv); // Assuming Minimum scale of timing ~100 ps = 0.1 ns
    if (nInT < iMnT) { iMnT = nInT;} 
    if (nInT > iMxT) { iMxT = nInT;}
    
    double gapX = (pargas[0] + inocal0hit_pointer->InoCal0Hit_list[ij]->GetLocalXPos() - nInX[0]*Xstrwd)/Xstrwd  - 0.5;
    double gapY = (pargas[1] + inocal0hit_pointer->InoCal0Hit_list[ij]->GetLocalYPos() - nInY[0]*Ystrwd)/Ystrwd  - 0.5;

    int nxmul=1;
    int nymul=1;
    if(pAnalysis->collatedIn) {
      double CorrIneffiParX = pAnalysis->inefficiency_uncx[nInLA]->GetBinContent(nInX[0]+1,nInY[0]+1);
      double CorrIneffiParY = pAnalysis->inefficiency_uncy[nInLA]->GetBinContent(nInX[0]+1,nInY[0]+1);
      if(gRandom->Rndm(0) < CorrIneffiParX) nInX[0] = -1;
      if(gRandom->Rndm(0) < CorrIneffiParY) nInY[0] = -1;
    }
    if(nInX[0] <0 && nInY[0] <0) { continue;}
    // cout<<"nInX[0] "<<nInX[0] <<" "<<nInY[0]<<" "<<nInLA<<endl;
    
    if(nInX[0] >=0 && NewMultiplicity) {
      if(pAnalysis->collatedIn) {
	// nxmul = GetRandomXY(gapX,pAnalysis->strp_xmulsim_cor[nInLA]);
	nxmul = GetRandomXY(gapX,pAnalysis->block_xmulsim[nInLA][int(nInX[0]/4.)][int(nInY[0]/4.)]);
	// cout << " something x " << nxmul << endl;
      } else {
	double arand=gRandom->Rndm();
	if (arand<0.1) { //10% case three strip hits
	  nxmul=3;
	} else { 
	  if (gRandom->Rndm() < 3.2*gapX*gapX) {
	    nxmul=2;
	  } else {
	    nxmul=1;
	  }
	}
      }
      // nxmul = 1;
      if (nxmul==3) { //10% case three strip hits
	nInX[1] = nInX[0] + 1; //int(gapX/abs(max(1.e-12,gapX)));
	nInX[2] = nInX[0] - 1; //int(gapX/abs(max(1.e-12,gapX)));
      } else if (nxmul==2) { 
	// f(x) = ax**2, => a=3.2
	nInX[1] = nInX[0] + int(gapX/(max(1.e-12,abs(gapX))));
      } 
    } // if (nInX[0] >=0 && NewMultiplicity) {

    if (nInY[0] >=0 && NewMultiplicity) {
      if(pAnalysis->collatedIn) {
	// nymul = GetRandomXY(gapY,pAnalysis->strp_ymulsim_cor[nInLA]);
	nymul = GetRandomXY(gapY,pAnalysis->block_ymulsim[nInLA][int(nInX[0]/4.)][int(nInY[0]/4.)]);
	// cout << " something y " << nymul << endl;
      } else {
	double arand=gRandom->Rndm();
	if (arand<0.1) {nymul = 3;}
	else { 
	  if (gRandom->Rndm() < 3.2*gapY*gapY) {nymul=2;}
	  else {nymul=1;}
	}
      }
      // nymul = 1;
      if (nymul==3) { //10% case three strip hits
	nInY[1] = nInY[0] + 1; //int(gapY/abs(max(1.e-12,gapY)));
	nInY[2] = nInY[0] - 1; //int(gapY/abs(max(1.e-12,gapY)));
      } else if (nymul==2) { 
	nInY[1] = nInY[0] + int(gapY/(max(1.e-12,abs(gapY))));
      }
    }

    // cout<<"nInX[0] "<<nInX[0] <<" "<<nInY[0]<<" "<<nInLA<<endl;
    
    for (int ix=0; ix<MxStrip; ix++) { 
      if(!NewMultiplicity && ix>0) continue;
      if (nInX[ix] <0 || nInX[ix]>=numberInX) continue;

      // if(pAnalysis->collatedIn && nInLA!=5) {
      // 	UnCorrXIneffiPar = pAnalysis->inefficiency_uncx[nInLA]->GetBinContent(nInX[ix]+1,nInY[0]+1);
      // }
      // if(gRandom->Rndm(0) < UnCorrXIneffiPar) continue;
      
      // double trigeffiX = 0.0;
      // if(pAnalysis->collatedIn) {
      // 	trigeffiX = pAnalysis->triggereffi_xevt[nInLA]->GetBinContent(nInX[ix]+1,nInY[0]+1);
      // 	for(int trglx=0; trglx<ntriglay; trglx++) {
      // 	  if((nInLA == TrgLayer[trglx]) && (gRandom->Rndm()<(trigeffiX))) {
      // 	    TrgDataX[TrgLayer[trglx]]++;
      // 	  }
      // 	}
      // } // if(pAnalysis->collatedIn) {
      
      G4double atimeX = tmpatimeX +  gRandom->Gaus(0,TimeUnCorrSmr);
      int nInTX = int(atimeX/TimeToDigiConv);
      int iold = 0;
	    
      for (unsigned jk=0; jk<inoStripX_pointer->InoStripX_list.size(); jk++) {
	InoStrip* Xstrip =inoStripX_pointer->InoStripX_list[jk]; 
	if (Xstrip->GetRPCmod()==iRPCMod && 
	    Xstrip->GetStrip()%numberInX==nInX[ix]) { 
	  inoStripX_pointer->InoStripX_list[jk]->AddPulse(edep); //GMA151001 for large multiplicty share this energy
	  if (inoStripX_pointer->InoStripX_list[jk]->GetSmrTime() >nInTX) {
	    inoStripX_pointer->InoStripX_list[jk]->SetSmrTime(nInTX);
	  }
	  iold = 1; break;
	}
      }
	    
      //GMA Take precaution of these one strip hits
      // 1. Segement direction, for X/Y Z-value might be different 
      //    consequently direction
      // 2. 
      //
      //
      if (iold==0) {
	      
	InoStrip *  Xstrip = new InoStrip(); //VALGRIND
	Xstrip->SetStrip(numberInX*numberInMO*nInDT+numberInX*nInMO+nInX[ix]);
	      
	Xstrip->SetpdgId( pdgid);
	Xstrip->SetTrueTime(nInT);
	Xstrip->SetSmrTime(nInTX); //nInT
	      
	// cout <<"xtime "<< nInT<<" "<<nInTX<<" "<<Xstrip->GetTrueTime()<<" "<<Xstrip->GetSmrTime()<<" "<<pdgid<<" "<<Xstrip->GetpdgId()<<endl;
	      
	Xstrip->SetPulse(edep);
	Xstrip->SetRPCmod(iRPCMod);
	    
	G4ThreeVector trkmom = inocal0hit_pointer->InoCal0Hit_list[ij]->GetMom();
	    
	Xstrip->SetMomentum(trkmom.mag());
	Xstrip->SetTheta(trkmom.theta());
	Xstrip->SetPhi(trkmom.phi());
	      
	G4ThreeVector posvec3 = inocal0hit_pointer->InoCal0Hit_list[ij]->GetPos();
	Xstrip->SetGenPosX(posvec3.x());
	Xstrip->SetGenPosY(posvec3.y());
	Xstrip->SetGenPosZ(posvec3.z());
	      
	G4int xstripid = 0;
	xstripid<<=2;
	xstripid +=nInDT;
	      
	xstripid<<=8;
	xstripid +=nInLA;
	      
	xstripid<<=3;
	xstripid +=nInMO;
	      
	xstripid<<=3;
	xstripid +=nInCH;
	      
	xstripid<<=7;
	xstripid +=nInX[ix];
	      
	xstripid<<=5;
	xstripid +=0; //nInT;   
	      
	xstripid<<=3;
	xstripid +=TMath::Min(int(edep/16),7);
	      
	Xstrip->SetId(xstripid);
	    
	//cout <<"iold = "<<iold<<", nInLA = "<<nInLA<<", nInX = "<<nInX<<endl;
	// pAnalysis->DeadStripX->Fill(nInX[ix]);
	// cout<<"pAnalysis->DeadStripX->Fill(nInX);="<<nInX<<endl;
	inoStripX_pointer->InoStripX_list.push_back(Xstrip);
      } // if (iold==0)
    } // for (int ix=0; ix<MxStrip; ix++)

    for (int jy=0; jy<MxStrip; jy++) { 
      if(!NewMultiplicity && jy>0) continue;
      if (nInY[jy] <0 || nInY[jy]>=numberInY) continue;

      // if(pAnalysis->collatedIn && nInLA!=5) {
      // 	UnCorrYIneffiPar = pAnalysis->inefficiency_uncy[nInLA]->GetBinContent(nInX[0]+1,nInY[jy]+1);
      // }
      // if(gRandom->Rndm(0) < UnCorrYIneffiPar) continue;

      // double trigeffiY = 0.0;
      // if(pAnalysis->collatedIn) {
      // 	trigeffiY = pAnalysis->triggereffi_yevt[nInLA]->GetBinContent(nInX[0]+1,nInY[jy]+1);
      // 	for(int trgly=0; trgly<ntriglay; trgly++) {
      // 	  if((nInLA == TrgLayer[trgly]) && (gRandom->Rndm()<(trigeffiY))) {
      // 	    TrgDataY[TrgLayer[trgly]]++;
      // 	  }
      // 	}
      // } // if(pAnalysis->collatedIn) {
      
      G4double atimeY = tmpatimeY +  gRandom->Gaus(0,TimeUnCorrSmr);
      int nInTY = int(atimeY/TimeToDigiConv);
	    
      int iold = 0;
      for (unsigned jk=0; jk<inoStripY_pointer->InoStripY_list.size(); jk++) {
	InoStrip* Ystrip =inoStripY_pointer->InoStripY_list[jk]; 
	if (Ystrip->GetRPCmod()==iRPCMod && 
	    Ystrip->GetStrip()%numberInY==nInY[jy]) { 
	  inoStripY_pointer->InoStripY_list[jk]->AddPulse(edep); 
	  if (inoStripY_pointer->InoStripY_list[jk]->GetSmrTime() >nInTY) {
	    inoStripY_pointer->InoStripY_list[jk]->SetSmrTime(nInTY);
	  }
	  iold = 1; break;
	}
      }
	    
      if (iold==0) {
	InoStrip* Ystrip = new InoStrip(); //VALGRIND
	Ystrip->SetStrip(numberInY*nInCH+nInY[jy]);
	      
	Ystrip->SetpdgId(pdgid);
	Ystrip->SetTrueTime(nInT);
	Ystrip->SetSmrTime(nInTY);
	// cout <<"ytime "<< nInT<<" "<<nInTY<<" "<<Ystrip->GetTrueTime()<<" "<<Ystrip->GetSmrTime()<<" "<<pdgid<<" "<<Ystrip->GetpdgId()<<endl;

	Ystrip->SetPulse(edep);
	Ystrip->SetRPCmod(iRPCMod);
	G4ThreeVector trkmom = inocal0hit_pointer->InoCal0Hit_list[ij]->GetMom();
	    
	Ystrip->SetMomentum(trkmom.mag());
	Ystrip->SetTheta(trkmom.theta());
	Ystrip->SetPhi(trkmom.phi());
	      
	G4ThreeVector posvec = inocal0hit_pointer->InoCal0Hit_list[ij]->GetPos();
	Ystrip->SetGenPosX(posvec.x());
	Ystrip->SetGenPosY(posvec.y());
	Ystrip->SetGenPosZ(posvec.z());
	      
	G4int ystripid = 1;
	      
	ystripid<<=2;
	ystripid +=nInDT;
	      
	ystripid<<=8;
	ystripid +=nInLA;
	      
	ystripid<<=3;
	ystripid +=nInMO;
	      
	ystripid<<=3;
	ystripid +=nInCH;
	    
	ystripid<<=7;
	ystripid +=nInY[jy];
	      
	ystripid<<=5;
	ystripid +=0; //05/01/2009 nInT;   
	      
	ystripid<<=3;
	ystripid +=TMath::Min(int(edep/16),7);
	      
	Ystrip->SetId(ystripid);

	//cout <<"iold = "<<iold<<", nInLA = "<<nInLA<<", nInY = "<<nInY[jy]<<endl;
	// pAnalysis->DeadStripY->Fill(nInY[jy]);
	// cout<<"pAnalysis->DeadStripY->Fill(nInY[jy]);="<<nInY<<endl;
	inoStripY_pointer->InoStripY_list.push_back(Ystrip);
      } // if (iold==0)
    } // for (int jy=0; jy<MxStrip; jy++)
  } // for (unsigned ij=0; ij<inocal0hit_pointer->InoCal0Hit_list.size(); ij++) {

  // cout<<"Digitising of data complete."<<endl;
  // cout<<"StripList sizes X: "<<inoStripX_pointer->InoStripX_list.size()<<" , Y: "<<inoStripY_pointer->InoStripY_list.size()<<endl;

}

void InoDigiAlg::CalculateTrigger() {
  for(int trgi=0; trgi<nLayer; trgi++) {
    if(TrgDataX[trgi]>0) {
      trigStoreX++;
    }
    if(TrgDataY[trgi]>0) {
      trigStoreY++;
    }
  }
}

void InoDigiAlg::SaveDigiData() {
  pAnalysis->pRootFile->cd();
  pAnalysis->trigx = trigStoreX;
  pAnalysis->trigy = trigStoreY;
  pAnalysis->triggeracceptance = (trigStoreX>3 || trigStoreY>3) ? 1 : 0;

  pAnalysis->ndigiht = inoStripX_pointer->InoStripX_list.size() 
    + inoStripY_pointer->InoStripY_list.size();
  if (pAnalysis->ndigiht >pAnalysis->ndigihtmx) pAnalysis->ndigiht =pAnalysis->ndigihtmx;
  for (unsigned ij=0; ij<inoStripX_pointer->InoStripX_list.size() && ij<pAnalysis->ndigiht ; ij++) {
    pAnalysis->stripid[ij] =inoStripX_pointer->InoStripX_list[ij]->GetId();
    pAnalysis->digipdgid[ij] =inoStripX_pointer->InoStripX_list[ij]->GetpdgId();
    pAnalysis->digitime[ij] = inoStripX_pointer->InoStripX_list[ij]->GetSmrTime();
    pAnalysis->digitruetime[ij] = inoStripX_pointer->InoStripX_list[ij]->GetTrueTime();
    pAnalysis->digienr[ij] =inoStripX_pointer->InoStripX_list[ij]->GetPulse();
    pAnalysis->digivx[ij] =inoStripX_pointer->InoStripX_list[ij]->GetGenPosX();
    pAnalysis->digivy[ij] =inoStripX_pointer->InoStripX_list[ij]->GetGenPosY();
    pAnalysis->digivz[ij] =inoStripX_pointer->InoStripX_list[ij]->GetGenPosZ();
	
    G4ThreeVector trkmom(1000,1,1);
    trkmom.setMag(inoStripX_pointer->InoStripX_list[ij]->GetMomentum());
    trkmom.setTheta(inoStripX_pointer->InoStripX_list[ij]->GetTheta());
    trkmom.setPhi(inoStripX_pointer->InoStripX_list[ij]->GetPhi());

    pAnalysis->digipx[ij] = trkmom.x();
    pAnalysis->digipy[ij] = trkmom.y();
    pAnalysis->digipz[ij] = trkmom.z();
    if (ij >=pAnalysis->ndigihtmx) break; //redundant
    // cout<<"ij "<<ij<<" "<<pAnalysis->digivx[ij]<<" "<<pAnalysis->digivy[ij]<<" "<<pAnalysis->digivz[ij]<<" "<<(pAnalysis->stripid[ij]>>8)<<endl;
  }
  unsigned jk = inoStripX_pointer->InoStripX_list.size();
  for (unsigned ij=0; ij<inoStripY_pointer->InoStripY_list.size() && jk< pAnalysis->ndigiht; ij++, jk++) {
    pAnalysis->stripid[jk] =inoStripY_pointer->InoStripY_list[ij]->GetId();
    pAnalysis->digipdgid[jk] =inoStripY_pointer->InoStripY_list[ij]->GetpdgId();
    pAnalysis->digitime[jk] =inoStripY_pointer->InoStripY_list[ij]->GetSmrTime();
    pAnalysis->digitruetime[jk] = inoStripY_pointer->InoStripY_list[ij]->GetTrueTime();
    pAnalysis->digienr[jk] =inoStripY_pointer->InoStripY_list[ij]->GetPulse();
    pAnalysis->digivx[jk] =inoStripY_pointer->InoStripY_list[ij]->GetGenPosX();
    pAnalysis->digivy[jk] =inoStripY_pointer->InoStripY_list[ij]->GetGenPosY();
    pAnalysis->digivz[jk] =inoStripY_pointer->InoStripY_list[ij]->GetGenPosZ();
	
    G4ThreeVector trkmom(1000,0,0);
    trkmom.setMag(inoStripY_pointer->InoStripY_list[ij]->GetMomentum());
    trkmom.setTheta(inoStripY_pointer->InoStripY_list[ij]->GetTheta());
    trkmom.setPhi(inoStripY_pointer->InoStripY_list[ij]->GetPhi());
	
    pAnalysis->digipx[jk] = trkmom.x();
    pAnalysis->digipy[jk] = trkmom.y();
    pAnalysis->digipz[jk] = trkmom.z();
    if (jk >=pAnalysis->ndigihtmx) break; //redundant
    // cout<<"jk "<<jk<<" "<<pAnalysis->digivx[jk]<<" "<<pAnalysis->digivy[jk]<<" "<<pAnalysis->digivz[jk]<<" "<<(pAnalysis->stripid[jk]>>8)<<endl;
  }
  // cout<<"digioutput "<<pAnalysis->ndigiht<<" "<<pAnalysis->trigx<<" "<<pAnalysis->trigy<<endl;
  pAnalysis->pEventTree->Fill();
}

int InoDigiAlg::GetRandomXY(double& GapX, TH2D* tmphistx) {
  double sumX = gRandom->Rndm();
  int xbinf = tmphistx->GetXaxis()->FindBin(GapX);
  int nmult = -1;
  int iiter = 0;
  int nmxusedhits = 3;
  while(nmult<=0) {
    sumX = gRandom->Rndm();
    double valY = 0.0;
    for (int ijf=0; ijf<=nmxusedhits+1; ijf++) {
      valY += tmphistx->GetBinContent(xbinf, ijf+1);
      if (valY > sumX) {
	nmult = ijf; 
	break;
      }
    } // for (int ijf=0; ijf<=nmxusedhits; ijf++) {
    if (iiter++==100) {
      nmult = 1;
    }
  } // while(nmult<=0) { 

  return nmult;
}








  
