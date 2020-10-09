#include "InoCluster.h"
#include "InoHit.h"
#include "TMath.h"
#include <iostream>
InoCluster::InoCluster(InoHit* hit) :
  fZPlane(-1),fRPCmod(-1),fBegXStrip(-1), fEndXStrip(-1),fBegYStrip(-1), fEndYStrip(-1), 
  fBegTime(0.), fEndTime(0.),
  fBegXPos(-999.), fEndXPos(999.),
  fBegYPos(-999.), fEndYPos(999.),
  fXPos(-999.), fYPos(-999.), fZPos(0.), fXPulse(0.), fYPulse(0.), 
  fTrkFlag(0), fShwFlag(0),
  fTrkPlnFlag(0), fShwPlnFlag(0),
  //  fPlaneView(-1),
  fDigits(0),fNDFlag(1),       //asmQ what is the difference between fDigit and fView
  fXPosErr(999.),
  fYPosErr(999.),
  //  fXPosErr(0.02/sqrt(12.)),
  //  fYPosErr(0.02/sqrt(12.)),
  fView (2),
  InTrack(false), isStraight(true),
  StripXWidth(0.0196),
  StripYWidth(0.0196)
{
  paradef = DetectorParameterDef::AnPointer; //AAR: 
  StripXWidth = paradef->GetXStrwd()/1000;
  StripYWidth = paradef->GetYStrwd()/1000;
  //LayerThickness = (1/m)*2*(paradef->GetParlay(2)+paradef->GetParirlay(2)); //(1/m)*2*paradef->GetParlay(2);
  debug_clust = false;
  this->AddHit(hit);
}


InoCluster::~InoCluster()
{
  HitsInCluster.clear();
}


void InoCluster::AddHit(InoHit* hit) {
  if(HitsInCluster.size()==0) {
    HitsInCluster.push_back(hit);
    fZPlane=hit->GetZPlane();
    fRPCmod = hit->GetRPCmod();
    if (hit->GetXPosErr()<100) {
      fBegXStrip=hit->GetXStrip()->GetStrip();
      fEndXStrip=hit->GetXStrip()->GetStrip();
      fBegXPos=hit->GetXPos();
      fEndXPos=hit->GetXPos();
    }
    if (hit->GetYPosErr()<100) {
      fBegYStrip=hit->GetYStrip()->GetStrip();
      fEndYStrip=hit->GetYStrip()->GetStrip();
      fBegYPos=hit->GetYPos();
      fEndYPos=hit->GetYPos();
    }

    fBegTime=hit->GetTime();
    fEndTime=hit->GetTime();
    fZPos=hit->GetZPos();
    fView = hit->GetView();
    fXPos = hit->GetXPos();
    fYPos = hit->GetYPos();
    fXPosErr = hit->GetXPosErr();
    fYPosErr = hit->GetYPosErr();
    //    fPlaneView=hit->GetPlaneView();            //asmQ what is fView for why is it set to 2
  } else {
    if(this->ContainsHit(hit)==true) {return;}
    HitsInCluster.push_back(hit);   
    
    if (hit->GetXPosErr()<100) {
      if(hit->GetXStrip()->GetStrip()<fBegXStrip) fBegXStrip=hit->GetXStrip()->GetStrip();
      if(hit->GetXStrip()->GetStrip()>fEndXStrip) fEndXStrip=hit->GetXStrip()->GetStrip();
      if(hit->GetXPos()<fBegXPos) fBegXPos=hit->GetXPos();
      if(hit->GetXPos()>fEndXPos) fEndXPos=hit->GetXPos();
      if (fView ==1) fView = 2;
    }
    
    if (hit->GetYPosErr()<100) {
      if(hit->GetYStrip()->GetStrip()<fBegYStrip) fBegYStrip=hit->GetYStrip()->GetStrip();
      if(hit->GetYStrip()->GetStrip()>fEndYStrip) fEndYStrip=hit->GetYStrip()->GetStrip();
      if(hit->GetYPos()<fBegYPos) fBegYPos=hit->GetYPos();
      if(hit->GetYPos()>fEndYPos) fEndYPos=hit->GetYPos();
      if (fView ==0) fView = 2;
    }
    if(hit->GetTime()<fBegTime) fBegTime=hit->GetTime();
    if(hit->GetTime()>fEndTime) fEndTime=hit->GetTime();
  }
  
  //fDigit gives total # of X+Y  strips  //nXstrip -> total # of x strip in this cluster
  if (hit->GetXPosErr()<100) { 
    fDigits += 1; // hit->GetCandStripHandle()->GetNDaughters();
    unsigned int nxstrip = GetXEntries();
    fXPos = (fXPos*(nxstrip-1)+hit->GetXPos())/(1.0*nxstrip);
    if (nxstrip==1) {  
      fXPosErr = hit->GetXPosErr();
    } else {
      fXPosErr = hit->GetXPosErr()/2; //GMA 09/02/09 Put these separately from data   //asm why so?
    } 
    fXPulse += hit->GetXPulse();
  }
  
  if (hit->GetYPosErr()<100) { 
    fDigits += 1;       
    unsigned int nystrip = GetYEntries();
    fYPos = (fYPos*(nystrip-1)+hit->GetYPos())/(1.0*nystrip);
    if (nystrip==1) {
      fYPosErr = hit->GetYPosErr();
    } else {
      fYPosErr = hit->GetYPosErr()/2; //GMA 09/02/09 Put these separately from data
    } 
    fYPulse += hit->GetYPulse();
  }
  return;
}

bool InoCluster::ContainsHit(InoHit* hit) {
  for(unsigned int ij=0; ij<HitsInCluster.size(); ++ij) {
    if(hit==HitsInCluster[ij]) {return true;}
  }
  return false;
}

unsigned int InoCluster::GetXEntries() {
  unsigned int nxhit = 0;
  for(unsigned int ij=0; ij<HitsInCluster.size(); ++ij) {
    if (HitsInCluster[ij]->GetXPosErr() < 100) nxhit++;
  }
  return nxhit;
}

unsigned int InoCluster::GetYEntries() {
  unsigned int nyhit = 0;
  for(unsigned int ij=0; ij<HitsInCluster.size(); ++ij) {
    if (HitsInCluster[ij]->GetYPosErr() < 100) nyhit++;
  }
  return nyhit;
}

unsigned int InoCluster::GetNXStripsInClust() {
  // cout<<"InoCluster::GetNXStripsInClust()..."<<endl;
  vector <int> allrecord;
  allrecord.clear();
  for(unsigned int ij=0; ij<HitsInCluster.size(); ij++) {
    if(HitsInCluster[ij]->GetXPosErr() < 100) {
      if(allrecord.size()) {
	int ncnt=0;
	for(unsigned jk=0; jk<allrecord.size(); jk++) {
	  if(allrecord[jk]==HitsInCluster[ij]->GetXStripNum()) ncnt++;
	}
	if(ncnt==0) {allrecord.push_back(HitsInCluster[ij]->GetXStripNum());}
      } else {
	allrecord.push_back(HitsInCluster[ij]->GetXStripNum());
      }
    }
  }
  // cout<<".... InoCluster::GetNXStripsInClust()..."<<endl;
  return allrecord.size();
}

unsigned int InoCluster::GetNYStripsInClust() {
  // cout<<"InoCluster::GetNYStripsInClust()..."<<endl;
  vector <int> allrecord;
  allrecord.clear();
  for(unsigned int ij=0; ij<HitsInCluster.size(); ij++) {
    if(HitsInCluster[ij]->GetYPosErr() < 100) {
      if(allrecord.size()) {
	int ncnt=0;
	for(unsigned jk=0; jk<allrecord.size(); jk++) {
	  if(allrecord[jk]==HitsInCluster[ij]->GetYStripNum()) ncnt++;
	}
	if(ncnt==0) {allrecord.push_back(HitsInCluster[ij]->GetYStripNum());}
      } else {
	allrecord.push_back(HitsInCluster[ij]->GetYStripNum());
      }
    }
  }
  // cout<<".... InoCluster::GetNYStripsInClust()..."<<endl;
  return allrecord.size();
}

unsigned int InoCluster::GetXProjEntries() {
  unsigned int nxhit = 0; 
   
  for(unsigned int ij=0; ij<HitsInCluster.size(); ++ij) {
    bool is_identicl=0;
    if (HitsInCluster[ij]->GetXPosErr() < 100){ 
      for ( unsigned int jk=0; jk<ij;jk++){
	if( HitsInCluster[ij]->GetXStripNum() == HitsInCluster[ij]->GetXStripNum()){ 
	  is_identicl=1; break;
	}
      }
      if(!is_identicl)nxhit++; 
    }
  }
  return nxhit;
}
unsigned int InoCluster::GetYProjEntries() {
  unsigned int nyhit = 0;
  
  for(unsigned int ij=0; ij<HitsInCluster.size(); ++ij) {
    bool is_identicl=0;
    if (HitsInCluster[ij]->GetYPosErr() < 100){
      for ( unsigned int jk=0; jk<ij;jk++){
	if( HitsInCluster[ij]->GetYStripNum() == HitsInCluster[ij]->GetYStripNum()) { 
	  is_identicl=1; break;
	}
      }
      if(!is_identicl)nyhit++;
    }
  }
  return nyhit;
}

int InoCluster::IsHitAssoc(InoHit* hit) const {
  double TimeWindow = 5.0;//9999.9; //GMA 290709 Optimise this
  
  if ((hit->GetZPlane()==fZPlane) &&
      (hit->GetXPosErr()>100 || (hit->GetXPos()>(fBegXPos-2.1*StripXWidth) && hit->GetXPos()<(fEndXPos+2.1*StripXWidth))) &&
      (hit->GetYPosErr()>100 || (hit->GetYPos()>(fBegYPos-2.1*StripYWidth) && hit->GetYPos()<(fEndYPos+2.1*StripYWidth))) &&
      (hit->GetTime()>(fBegTime-TimeWindow) && hit->GetTime()<(fEndTime+TimeWindow))) {
    return 1; 
  }  else {
    return 0;
  }
}

int InoCluster::IsShwAssoc(InoCluster* clust) const {
  //Need to check for signle plane hits
  double TimeWindow = 99.9; int ShwAssocNum = 0;
  
  if( (fEndTime-clust->GetBegTime())>-TimeWindow || (clust->GetEndTime()-fBegTime)>-TimeWindow ) {
    
    if( abs(clust->GetZPlane()-fZPlane)<5 &&
	(clust->GetXPos()<-100 ||
	 ((clust->GetEndXPos()-fBegXPos)>-6*StripXWidth && (fEndXPos-clust->GetBegXPos())>-6*StripXWidth)) &&
	(clust->GetYPos()<-100 ||
	 ((clust->GetEndYPos()-fBegYPos)>-6*StripYWidth && (fEndYPos-clust->GetBegYPos())>-6*StripYWidth))) {
      
      if( ( abs(clust->GetZPlane()-fZPlane)<3  &&
	    (clust->GetXPos()<-100 || 
	     ((clust->GetEndXPos()-fBegXPos)>-StripXWidth && (fEndXPos-clust->GetBegXPos())>-StripXWidth)) &&
	    (clust->GetYPos()<-100 || 
	     ((clust->GetEndYPos()-fBegYPos)>-StripYWidth && (fEndYPos-clust->GetBegYPos())>-StripYWidth)))
	  
          || ( clust->GetZPlane()==fZPlane &&
	       (clust->GetXPos()<-100 || 
	        ((clust->GetEndXPos()-fBegXPos)>-3*StripXWidth && (fEndXPos-clust->GetBegXPos())>-3*StripXWidth)) &&
	       (clust->GetYPos()<-100 || 
               ((clust->GetEndYPos()-fBegYPos)>-3*StripYWidth && (fEndYPos-clust->GetBegYPos())>-3*StripYWidth ))) ) {
	ShwAssocNum=2;
      } else {
	ShwAssocNum=1;
      }
    }
  }
  
  if(ShwAssocNum==2 && this->GetHitEntries()<3 && clust->GetHitEntries()<2) {ShwAssocNum=1;}
  
  return ShwAssocNum;
}


int InoCluster::IsTrkAssoc(InoCluster* clustm, InoCluster* clustp) const {
  double TimeWindow = 15.0; int TrkAssocNum = 0;
  double NDScale=2; //GMAA 1;
  //  double fact; //asm
  //  double min=0.2; double max=0.8;
  
  // Configure for correct detector instrumentation
  int PlaneGap = 1;
  
  // Check timing proximity
  if(debug_clust) {
    cout<<"-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-b-"<<endl;
    cout<<"ij-1 Plane =  z "<<clustm->GetZPlane() << " x " <<clustm->GetXPos() <<" y "<<clustm->GetYPos()<<endl;
    cout<<"ij Plane =  z "<< fZPlane << " x " << fXPos <<" y "<< fYPos <<endl;
    cout<<"ij+1 Plane =  z "<<clustp->GetZPlane() << " x " <<clustp->GetXPos() <<" y "<<clustp->GetYPos()<<endl;
  }
  if (fXPos <-100 && clustm->GetXPos()<-100 &&  clustp->GetXPos()<-100) {
    if(debug_clust) cout<<"Abnormal X return 0;"<<endl;
    return TrkAssocNum;
  }
  if (fYPos <-100 && clustm->GetYPos()<-100 &&  clustp->GetYPos()<-100) {
    if(debug_clust) cout<<"Abnormal Y return 0;"<<endl;
    return TrkAssocNum;
  }
  
  if(debug_clust) cout<<"ij - 1 fBegTime "<<clustm->GetBegTime()<<" fEndTime "<<clustm->GetEndTime()<<endl;
  if(debug_clust) cout<<"ij fBegTime "<<fBegTime<<" fEndTime "<<fEndTime<<endl;
  if(debug_clust) cout<<"ij + 1 fBegTime "<<clustp->GetBegTime()<<" fEndTime "<<clustp->GetEndTime()<<endl;
  if(debug_clust) cout<<"Condn Time : "<< ((fEndTime-clustm->GetBegTime())>-TimeWindow) << ((clustm->GetEndTime()-fBegTime)>-TimeWindow) << ((fEndTime-clustp->GetBegTime())>-TimeWindow) << ((clustp->GetEndTime()-fBegTime)>-TimeWindow) << " " <<((fEndTime-clustm->GetBegTime())>-TimeWindow && (clustm->GetEndTime()-fBegTime)>-TimeWindow) << ((fEndTime-clustp->GetBegTime())>-TimeWindow && (clustp->GetEndTime()-fBegTime)>-TimeWindow) << " " << (((fEndTime-clustm->GetBegTime())>-TimeWindow && (clustm->GetEndTime()-fBegTime)>-TimeWindow) && ((fEndTime-clustp->GetBegTime())>-TimeWindow && (clustp->GetEndTime()-fBegTime)>-TimeWindow)) << endl;
  if(( (fEndTime-clustm->GetBegTime())>-TimeWindow && (clustm->GetEndTime()-fBegTime)>-TimeWindow) &&
     ( (fEndTime-clustp->GetBegTime())>-TimeWindow && (clustp->GetEndTime()-fBegTime)>-TimeWindow) ) {
    
    // If more than two planes away, scale back width of cluster
    // and then treat as if only two planes away
    
    double mXPos=clustm->GetXPos();
    double mYPos=clustm->GetYPos();
    
    if((fZPlane-clustm->GetZPlane())>=PlaneGap) {
      double mScale = double(PlaneGap)/double(fZPlane-clustm->GetZPlane());
      
      mXPos=fXPos+mScale*(mXPos-fXPos);
      mYPos=fYPos+mScale*(mYPos-fYPos);
    }
    
    double pXPos=clustp->GetXPos();
    double pYPos=clustp->GetYPos();
    
    if((clustp->GetZPlane()-fZPlane)>=PlaneGap) {
      double pScale = double(PlaneGap)/double(clustp->GetZPlane()-fZPlane);
      
      pXPos=fXPos+pScale*(pXPos-fXPos);
      pYPos=fYPos+pScale*(pYPos-fYPos);
    }

    //asm[][][][][][][][][][][][][]This function is changed to get tracks with gracing angle[][][][]
    if(debug_clust) cout<<"Before Condn. TrkAssocNum = "<<TrkAssocNum<<endl;
    if(debug_clust) cout<<"Condn. Pos X "<< (fXPos >-100) << (clustm->GetXPos()>-100) << (clustp->GetXPos()>-100) << (TMath::Abs(mXPos+pXPos-2*fXPos) < NDScale*1.1*StripXWidth) << (TMath::Abs(mXPos-fXPos) < 30*StripXWidth) << (TMath::Abs(pXPos-fXPos) < 30*StripXWidth) <<" "<< (fXPos >-100 && clustm->GetXPos()>-100 &&  clustp->GetXPos()>-100 &&	TMath::Abs(mXPos+pXPos-2*fXPos) < NDScale*1.1*StripXWidth && TMath::Abs(mXPos-fXPos) < 30*StripXWidth && TMath::Abs(pXPos-fXPos) < 30*StripXWidth) << endl;
    if(debug_clust) cout<<"(TMath::Abs(mXPos+pXPos-2*fXPos) < NDScale*1.1*StripXWidth) "<< TMath::Abs(mXPos+pXPos-2*fXPos) << " " << NDScale << " " << StripXWidth << " " << NDScale*1.1*StripXWidth << endl;
    
    if (fXPos >-100 && clustm->GetXPos()>-100 &&  clustp->GetXPos()>-100) {
      if (  TMath::Abs(mXPos+pXPos-2*fXPos) < NDScale*1.1*StripXWidth &&
	    TMath::Abs(mXPos-fXPos) < 10*StripXWidth &&
	    TMath::Abs(pXPos-fXPos) < 10*StripXWidth) {
	TrkAssocNum = 1;
	if(debug_clust) cout<<"Condn. X 1"<<endl;
      } else if (TMath::Abs(mXPos-fXPos) < 5*StripXWidth &&
		 TMath::Abs(pXPos-fXPos) < 5*StripXWidth) {
	TrkAssocNum = 1;
	if(debug_clust) cout<<"Condn. X 2"<<endl;
      }
    } else {
      if(debug_clust) cout<<"Condn. X not satisfied."<<endl;
    }
   // if (fXPos >-100 && clustm->GetXPos()>-100 &&  clustp->GetXPos()>-100 &&
   // 	// TMath::Abs(mXPos+pXPos-2*fXPos) < NDScale*fact*StripXWidth &&
   // 	TMath::Abs(mXPos+pXPos-2*fXPos) < NDScale*1.1*StripXWidth &&
   // 	TMath::Abs(mXPos-fXPos) < 30*StripXWidth &&
   // 	TMath::Abs(pXPos-fXPos) < 30*StripXWidth) {
   //    TrkAssocNum = 1;
   //    if(debug_clust) cout<<"CriX1"<<endl;
   //  } else if (fXPos >-100 && clustm->GetXPos()>-100) {
   //    if(debug_clust) cout<<"Condn. X 2"<< endl;
   //    if(debug_clust) cout<<"TMath::Abs(mXPos-fXPos) "<< TMath::Abs(mXPos-fXPos) << " " << StripXWidth << " " << 5*StripXWidth << endl;
   //    if (TMath::Abs(mXPos-fXPos) < 5*StripXWidth) {
   // 	TrkAssocNum = 1;
   // 	if(debug_clust) cout<<"CriX2"<<endl;
   //    }
   //  } else if (fXPos >-100 && clustp->GetXPos()>-100) {
   //    if(debug_clust) cout<<"Condn. X 3"<< endl;
   //    if (TMath::Abs(pXPos-fXPos) < 5*StripXWidth) {
   // 	TrkAssocNum = 1;
   // 	if(debug_clust) cout<<"CriX3"<<endl;
   //    }
   //  } else if (clustm->GetXPos() >-100 && clustp->GetXPos()>-100) {
   //    if(debug_clust) cout<<"Condn. X 4"<< endl;
   //    if (TMath::Abs(pXPos-mXPos) < 5*StripXWidth) {
   // 	TrkAssocNum = 1;
   // 	if(debug_clust) cout<<"CriX4"<<endl;
   //    }
   //  } else {
   //    TrkAssocNum = 1;
   //    if(debug_clust) cout<<"else X"<<endl;
   //  }
    if(debug_clust) cout<<"After Condn. Pos X TrkAssocNum = "<< TrkAssocNum <<endl;
    
    if(debug_clust) cout<<"Condn. Pos Y "<< (fYPos >-100) << (clustm->GetYPos()>-100) << ( clustp->GetYPos()>-100) << (TMath::Abs(mYPos+pYPos-2*fYPos) < NDScale*1.1*StripYWidth ) << (TMath::Abs(mYPos-fYPos) < 30*StripYWidth) << (TMath::Abs(pYPos-fYPos) < 30*StripYWidth) << " " << (fYPos >-100 && clustm->GetYPos()>-100 &&  clustp->GetYPos()>-100 && TMath::Abs(mYPos+pYPos-2*fYPos) < NDScale*1.1*StripYWidth  && TMath::Abs(mYPos-fYPos) < 30*StripYWidth && TMath::Abs(pYPos-fYPos) < 30*StripYWidth) << endl;
    if(debug_clust) cout<<"(TMath::Abs(mYPos+pYPos-2*fYPos) < NDScale*1.1*StripYWidth) "<< TMath::Abs(mYPos+pYPos-2*fYPos) << " " << NDScale << " " << StripYWidth << " " << NDScale*1.1*StripYWidth << endl;

    if (fYPos >-100 && clustm->GetYPos()>-100 &&  clustp->GetYPos()>-100) {
      if (  TMath::Abs(mYPos+pYPos-2*fYPos) < NDScale*1.1*StripYWidth &&
	    TMath::Abs(mYPos-fYPos) < 10*StripYWidth &&
	    TMath::Abs(pYPos-fYPos) < 10*StripYWidth) {
	if(debug_clust) cout<<"Condn. Y 1"<<endl;
	TrkAssocNum += 2;
      } else if (TMath::Abs(mYPos-fYPos) < 5*StripYWidth &&
		 TMath::Abs(pYPos-fYPos) < 5*StripYWidth) {
	TrkAssocNum += 2;
	if(debug_clust) cout<<"Condn. Y 2"<<endl;
      }
    } else {
      if(debug_clust) cout<<"Condn. Y not satisfied."<<endl;
    }
    // if (fYPos >-100 && clustm->GetYPos()>-100 &&  clustp->GetYPos()>-100 &&
    // 	// TMath::Abs(mYPos+pYPos-2*fYPos) < NDScale*fact*StripYWidth &&
    // 	TMath::Abs(mYPos+pYPos-2*fYPos) < NDScale*1.1*StripYWidth  &&
    // 	TMath::Abs(mYPos-fYPos) < 30*StripYWidth &&
    // 	TMath::Abs(pYPos-fYPos) < 30*StripYWidth)  {
    //   TrkAssocNum += 2;
    //   if(debug_clust) cout<<"CriY1"<<endl;
    // } else if (fYPos >-100 && clustm->GetYPos()>-100) {
    //   if(debug_clust) cout<<"Condn. Y 2"<< endl;
    //   if(debug_clust) cout<<"TMath::Abs(mYPos-fYPos) "<< TMath::Abs(mYPos-fYPos) <<endl;
    //   if (TMath::Abs(mYPos-fYPos) < 5*StripYWidth) {
    // 	TrkAssocNum += 2;
    // 	if(debug_clust) cout<<"CriY2"<<endl;
    //   }
    // } else if (fYPos >-100 && clustp->GetYPos()>-100) {
    //   if(debug_clust) cout<<"Condn. Y 3"<< endl;
    //   if (TMath::Abs(pYPos-fYPos) < 5*StripYWidth) {
    // 	TrkAssocNum += 2;
    // 	if(debug_clust) cout<<"CriY3"<<endl;
    //   }
    // } else if (clustm->GetYPos() >-100 && clustp->GetYPos()>-100) {
    //   if(debug_clust) cout<<"Condn. Y 4"<< endl;
    //   if (TMath::Abs(pYPos-mYPos) < 5*StripYWidth) {
    // 	TrkAssocNum += 2;
    // 	if(debug_clust) cout<<"CriY4"<<endl;
    //   }
    // } else {
    //   TrkAssocNum += 2;
    //   if(debug_clust) cout<<"else Y"<<endl;
    // }
    if(debug_clust) cout<<"After Condn. Pos Y TrkAssocNum = "<< TrkAssocNum <<endl;
    
    // if (TrkAssocNum!=3)// {
    //   NDScale= 3;    
    //   if(debug_clust) cout<<"(TrkAssocNum!=3) "<<endl;
    //   if(TrkAssocNum==1 || TrkAssocNum==0){
    // 	if (TrkAssocNum==0) NDScale= 2;  
    // 	if(debug_clust) cout<<"{{ NDScale = "<<NDScale<<endl;
    // 	if(debug_clust) cout << "NDScale*3.1*StripYWidth "<< NDScale*3.1*StripYWidth << endl;
    // 	if(debug_clust) cout << "Condition "<< (TMath::Abs(mYPos+pYPos-2*fYPos) < NDScale*3.1*StripYWidth) << ( fYPos >-100 && clustm->GetYPos()>-100 &&  clustp->GetYPos()>-100 && TMath::Abs(mYPos+pYPos-2*fYPos) < NDScale*3.1*StripYWidth && ((TMath::Abs(mYPos-fYPos) >= -25*StripYWidth && TMath::Abs(mYPos-fYPos) < 50*StripYWidth  && TMath::Abs(pYPos-fYPos) >= -25*StripYWidth && TMath::Abs(pYPos-fYPos) < 50*StripYWidth)|| (TrkAssocNum==0&&TMath::Abs(mYPos-fYPos) < 50*StripYWidth && TMath::Abs(pYPos-fYPos)<50*StripYWidth)))<< endl;
    // 	if( fYPos >-100 && clustm->GetYPos()>-100 &&  clustp->GetYPos()>-100 &&
    //         TMath::Abs(mYPos+pYPos-2*fYPos) < NDScale*3.1*StripYWidth &&
    // 	    ((TMath::Abs(mYPos-fYPos) >= -25*StripYWidth && TMath::Abs(mYPos-fYPos) < 50*StripYWidth  &&
    // 	      TMath::Abs(pYPos-fYPos) >= -25*StripYWidth && TMath::Abs(pYPos-fYPos) < 50*StripYWidth)||
    // 	     (TrkAssocNum==0&&TMath::Abs(mYPos-fYPos) < 50*StripYWidth && TMath::Abs(pYPos-fYPos)<50*StripYWidth)))
    // 	  { 
    // 	    TrkAssocNum += 2; 
    // 	  } 
    //   }
    //   if(debug_clust) cout<<"After Condn. Pos Y 11TrkAssocNum = "<< TrkAssocNum <<endl;
    //   if (TrkAssocNum==2 ){
    // 	if(debug_clust) cout<<"{{ NDScale = "<<NDScale<<endl;
    // 	if(debug_clust) cout << "Condition " << (TMath::Abs(mXPos+pXPos-2*fXPos) < NDScale*3.1*StripXWidth) << (fXPos >-100 && clustm->GetXPos()>-100 &&  clustp->GetXPos()>-100 && TMath::Abs(mXPos+pXPos-2*fXPos) < NDScale*3.1*StripXWidth && ((TMath::Abs(mXPos-fXPos) >=-25*StripXWidth && TMath::Abs(mXPos-fXPos) < 50*StripXWidth && TMath::Abs(pXPos-fXPos) >=-25*StripXWidth && TMath::Abs(pXPos-fXPos) < 50*StripYWidth )|| (NDScale==2&&TMath::Abs(mYPos-fYPos) < 50*StripYWidth&&TMath::Abs(pYPos-fYPos)<50*StripYWidth))) << endl;
    // 	if (fXPos >-100 && clustm->GetXPos()>-100 &&  clustp->GetXPos()>-100 &&
    // 	    TMath::Abs(mXPos+pXPos-2*fXPos) < NDScale*3.1*StripXWidth &&
    // 	    ((TMath::Abs(mXPos-fXPos) >=-25*StripXWidth && TMath::Abs(mXPos-fXPos) < 50*StripXWidth &&
    // 	      TMath::Abs(pXPos-fXPos) >=-25*StripXWidth && TMath::Abs(pXPos-fXPos) < 50*StripYWidth )||
    // 	     (NDScale==2&&TMath::Abs(mYPos-fYPos) < 50*StripYWidth&&TMath::Abs(pYPos-fYPos)<50*StripYWidth))) 
    // 	  {
    // 	    TrkAssocNum=3; 
    // 	  } 
    //   }
    //   NDScale=2;  
    // }
    
    // if (TrkAssocNum<3)// {
    //   TrkAssocNum = 0;
    //   if(debug_clust) cout<<"Condition 333 "<< NDScale <<endl;
      
    //   if((TMath::Abs(mYPos-fYPos) >=-50*StripYWidth && TMath::Abs(mYPos-fYPos) < 100*StripYWidth  &&
    // 	  TMath::Abs(pYPos-fYPos) >=-50*StripYWidth && TMath::Abs(pYPos-fYPos) < 100*StripYWidth  &&
    // 	  TMath::Abs(pXPos-fXPos) < 100*StripYWidth  ) ||
    // 	 (TMath::Abs(mXPos-fXPos) >=-50*StripXWidth && TMath::Abs(mXPos-fXPos)< 100*StripXWidth &&
    // 	  TMath::Abs(pXPos-fXPos) >=-50*StripXWidth && TMath::Abs(pXPos-fXPos)< 100*StripXWidth &&
    // 	  TMath::Abs(pYPos-fYPos) < 100*StripYWidth )){
	
    // 	if (fYPos >-100 && clustm->GetYPos()>-100 &&  clustp->GetYPos()>-100 &&
    // 	    TMath::Abs(mYPos+pYPos-2*fYPos) > NDScale*3.1*StripYWidth &&
    // 	    TMath::Abs(mYPos+pYPos-2*fYPos) < NDScale*8.1*StripYWidth) {
    // 	  TrkAssocNum = 1; 
    // 	}
    //     if (fXPos >-100 && clustm->GetXPos()>-100 &&  clustp->GetXPos()>-100 &&
    // 	    TMath::Abs(mXPos+pXPos-2*fXPos) > NDScale*3.1*StripXWidth && 
    // 	    TMath::Abs(mXPos+pXPos-2*fXPos) < NDScale*8.1*StripXWidth) {
    // 	  TrkAssocNum += 2; 
    // 	}
    //   }
    // }
    //asm[][][][][][][][][][][][][][][][][][][][]Over[][][][][][][][][][][]
  }

  //  cout <<"xxpostyyyy "<< mXPos<<" "<< fXPos<<" "<<pXPos<<" "<<mYPos<<" "<< mYPos<<" "<<pYPos<<endl;
  if(debug_clust) cout<<"In the end : IsTrkAssoc = "<< TrkAssocNum <<endl;
  return TrkAssocNum;
}

InoHit* InoCluster::GetHit(unsigned int ij) const {
  if(ij<HitsInCluster.size()) {return HitsInCluster[ij];}
  else {return 0;}
}

int InoCluster::IsDiffuseShwAssoc(InoCluster* clr) const {
  double win = 99.9;
  int assoc = 0;
  if(fEndTime-clr->GetBegTime()>-win || clr->GetEndTime()-fBegTime>-win){
    if( clr->GetZPlane()-fZPlane<9 && clr->GetZPlane()-fZPlane>-9 && 
	clr->GetEndXStrip()-fBegXStrip>-21 && fEndXStrip-clr->GetBegXStrip()>-21 &&
	clr->GetEndYStrip()-fBegYStrip>-21 && fEndYStrip-clr->GetBegYStrip()>-21 ) {
      if( ( clr->GetZPlane()-fZPlane<5 && clr->GetZPlane()-fZPlane>-5 && 
	    clr->GetEndXStrip()-fBegXStrip>-11 && fEndXStrip-clr->GetBegXStrip()>-11 && 
	    clr->GetEndYStrip()-fBegYStrip>-11 && fEndYStrip-clr->GetBegYStrip()>-11) ) assoc=2;
      else assoc=1;
    }
  }
  return assoc;
}

bool InoCluster::isIdentical(InoCluster* icls) {
  if ((GetHitEntries()) == (icls->GetHitEntries())) {
    for (unsigned ij=0; ij<GetHitEntries(); ij++) {
      if ( !GetHit(ij)->isIdentical( icls->GetHit(ij))) return false;
    }
  } else { return false ;}
  return true;
}
