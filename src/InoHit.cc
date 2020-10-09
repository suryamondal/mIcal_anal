//12/07/09 Put position of hit in this class from the strip information
// Where should we use mapping ? Electronic channel to geometrical channels 

#include "InoHit.h"
//#include "Validity/VldContext.h"
#include <cmath>
#include "TMath.h"
#include <iostream>

InoHit::InoHit() :
   fXStrip(0), fYStrip(0), fUid(0), 
   fXStripNum(-1), fXTime(-999.), fXTrueTime(-999.),
   fYStripNum(-1), fYTime(-999.), fYTrueTime(-999.), 
   fXpdgId(-25), fYpdgId(-25),
    fTrackFlag(0), fShowerFlag(0), fXPulse(0.), fXPlane(-1),
   fYPulse(0.), fYPlane(-1), fXPos(-999.), fXPosErr(999.),
   fYPos(-999.), fYPosErr(999.),
   genposx(-999.),genposy(-999.),
   fZPlane(-1), fZPos(-999.), fView(-1),
   fMomentum(0), fTheta(0), fPhi(0)
{}

InoHit::InoHit(InoStrip* fx, InoStrip* fy) :
   fXStrip(0), fYStrip(0), fUid(0), 
   fXStripNum(-1), fXTime(-999.), fXTrueTime(-999.),
   fYStripNum(-1), fYTime(-999.), fYTrueTime(-999.),
   fXpdgId(-25), fYpdgId(-25),
   fTrackFlag(0), fShowerFlag(0), fXPulse(0.), fXPlane(-1),
   fYPulse(0.), fYPlane(-1), fXPos(-999.), fXPosErr(999.),
   fYPos(-999.), fYPosErr(999.),
   genposx(-999.),genposy(-999.),
   fZPlane(-1), fZPos(-999.), fView(-1),
   fMomentum(0), fTheta(0), fPhi(0)
   {

  paradef = DetectorParameterDef::AnPointer; //AAR:
  pAnalysis = MultiSimAnalysisDigi::AnPointer;
  
  StripXWidth = paradef->GetXStrwd()/1000;
  StripYWidth = paradef->GetYStrwd()/1000;

  // StripXWidth = 0.01*sqrt(12)/1000;
  // StripYWidth = 0.01*sqrt(12)/1000;

  DigiToTimeConv = pAnalysis->GetTimeToDigiConvVal();
  SignalSpeed = pAnalysis->GetSignalSpeedVal();
  if (fx->GetPlane() != fy->GetPlane()) {
    std::cout <<"Strips are not in the same plane"<<std::endl;
  }
  if (fx->GetPlaneView()==fy->GetPlaneView()) {
    std::cout <<"Strips of parallel strips, donot fill up the variables"<<std::endl;
  } else if (fx->GetPlaneView()==0) {
    
    fXStripNum = fx->GetStrip(); //duplicate, but let it be
    // fXstriploc = fx->GetStripNumLoc();
    fXstriploc = (fx->GetId()>>8)&0x7F;
    fXTime = fx->GetSmrTime();
    fXTrueTime = fx->GetTrueTime();
    fXpdgId = fx->GetpdgId();
    fXPulse = fx->GetPulse();
    fXPlane = fx->GetPlaneView();
    fXPos = fx->GetXYPos();
    fXPosErr = StripXWidth/sqrt(12.); // GMA Need to put through database 
    fXStrip = fx;
    fZPos = fx->GetZPos();
    fZPlane = fx->GetPlane(); 
    genposx = fx->GetGenPosX();
    fYStripNum = fy->GetStrip(); //duplicate, but let it be
    // fYstriploc = fy->GetStripNumLoc();
    fYstriploc = (fy->GetId()>>8)&0x7F;
    fYTime = fy->GetSmrTime();
    fYTrueTime = fy->GetTrueTime();
    fYpdgId = fy->GetpdgId();
    fYPulse = fy->GetPulse();
    fYPlane = fy->GetPlaneView();
    fYPos = fy->GetXYPos();
    fYPosErr = StripYWidth/sqrt(12.); // GMA Need to put through database 
    fYStrip = fy;
    genposy = fy->GetGenPosY();
  } else {

    fXStripNum = fy->GetStrip(); //duplicate, but let it be
    // fXstriploc = fy->GetStripNumLoc();
    fXstriploc = (fy->GetId()>>8)&0x7F;
    fXTime = fy->GetSmrTime();
    fXTrueTime = fy->GetTrueTime();
    fXpdgId = fy->GetpdgId();
    fXPulse = fy->GetPulse();
    fXPlane = fy->GetPlaneView();
    fXPos = fy->GetXYPos();
    fXPosErr = StripXWidth/sqrt(12.); // GMA Need to put through database 
    fXStrip = fy;
    fZPos = fy->GetZPos();
    fZPlane = fy->GetPlane(); 
    genposx = fy->GetGenPosX();
    fYStripNum = fx->GetStrip(); //duplicate, but let it be
    fYstriploc = (fx->GetId()>>8)&0x7F;
    fYTime = fx->GetSmrTime();
    fYTrueTime = fx->GetTrueTime();
    fYpdgId = fx->GetpdgId();
    fYPulse = fx->GetPulse();
    fYPlane = fx->GetPlaneView();
    fYPos = fx->GetXYPos();
    fYPosErr = StripXWidth/sqrt(12.); // GMA Need to put through database 
    fYStrip = fx;
    genposy = fx->GetGenPosY();
  }
  fMomentum = fx->GetMomentum();
  fTheta = fx->GetTheta();
  fPhi = fx->GetPhi();
  fView = 2;
}

InoHit::InoHit(InoStrip* fx) :
   fXStrip(0), fYStrip(0), fUid(0), 
   fXStripNum(-1), fXTime(-999.), fXTrueTime(-999.),
   fYStripNum(-1), fYTime(-999.), fYTrueTime(-999.),
   fXpdgId(-25), fYpdgId(-25),
   fTrackFlag(0), fShowerFlag(0), fXPulse(0.), fXPlane(-1),
   fYPulse(0.), fYPlane(-1), fXPos(-999.), fXPosErr(999.),
   fYPos(-999.), fYPosErr(999.),
   fZPlane(-1), fZPos(-999.), fView(-1),
   fMomentum(0), fTheta(0), fPhi(0)
{

paradef = DetectorParameterDef::AnPointer; //AAR:
pAnalysis = MultiSimAnalysisDigi::AnPointer;

StripXWidth = paradef->GetXStrwd()/1000;
StripYWidth = paradef->GetYStrwd()/1000;
DigiToTimeConv = pAnalysis->GetTimeToDigiConvVal();
SignalSpeed = pAnalysis->GetSignalSpeedVal();

  if (fx->GetPlaneView()==0) {
    fXStripNum = fx->GetStrip(); //duplicate, but let it be
    fXstriploc = (fx->GetId()>>8)&0x7F;
    // fXstriploc = fx->GetStripNumLoc();
    fXTime = fx->GetSmrTime();
    fXTrueTime = fx->GetTrueTime();
    fXpdgId = fx->GetpdgId();
    fXPulse = fx->GetPulse();
    fXPlane = fx->GetPlaneView();
    fXPos = fx->GetXYPos();
    fXPosErr = 0.02/sqrt(12.); // Need to put through database   //AAR why is this 0.02 and not 0.0196
    fXStrip = fx;
    fZPos = fx->GetZPos();
    fZPlane = fx->GetPlane(); 
    fMomentum = fx->GetMomentum();
    fTheta = fx->GetTheta();
    fPhi = fx->GetPhi();
    genposx = fx->GetGenPosX();
    fView = 0;
  } else {
    fYStripNum = fx->GetStrip(); //duplicate, but let it be
    fYstriploc = (fx->GetId()>>8)&0x7F;
    fYTime = fx->GetSmrTime();
    fYTrueTime = fx->GetTrueTime();
    fYpdgId = fx->GetpdgId();
    fYPulse = fx->GetPulse();
    fYPlane = fx->GetPlaneView();
    fYPos = fx->GetXYPos();
    fYPosErr = 0.02/sqrt(12.); // Need to put through database 
    fYStrip = fx;
    fZPos = fx->GetZPos();
    fZPlane = fx->GetPlane(); 
    fMomentum = fx->GetMomentum();
    fTheta = fx->GetTheta();
    fPhi = fx->GetPhi();
    genposy = fx->GetGenPosY();
    fView = 1;
  }
}
// fView= 0=> X only, 1=>Y only , 2=> both X and Y
InoHit::~InoHit()
{
}


int InoHit::IsShwAssoc(InoHit* hit) const {
  
  cout<<" stripwdt:IsShwAssoc " << StripXWidth << " "<< StripYWidth<<endl;
  double win=99.9;
  if(fXStrip!=0 && fYStrip!=0) {
    if( TMath::Abs(hit->GetTime()-this->GetTime())<win ) {
      if( TMath::Abs(hit->GetZPlane()-this->GetZPlane())<3 && 
	  TMath::Abs(hit->GetXStrip()->GetStrip()-this->GetXStrip()->GetStrip())<4&&
	  TMath::Abs(hit->GetYStrip()->GetStrip()-this->GetYStrip()->GetStrip())<4) {
	return 2;
      } else if( TMath::Abs(hit->GetZPlane()-this->GetZPlane())<5 && 
		 TMath::Abs(hit->GetXStrip()->GetStrip()-this->GetXStrip()->GetStrip())<6 &&
		 TMath::Abs(hit->GetYStrip()->GetStrip()-this->GetYStrip()->GetStrip())<6) {
	return 1;
      }
    }
  } else if (fXStrip!=0 ) {
    if( TMath::Abs(hit->GetTime()-this->GetTime())<win ) {
      if( TMath::Abs(hit->GetZPlane()-this->GetZPlane())<3 && 
	  TMath::Abs(hit->GetXStrip()->GetStrip()-this->GetXStrip()->GetStrip())<4) {
	return 2;
      } else if( TMath::Abs(hit->GetZPlane()-this->GetZPlane())<5 && 
		 TMath::Abs(hit->GetXStrip()->GetStrip()-this->GetXStrip()->GetStrip())<6 ) {
	return 1;
      }
    }
  } else if (fYStrip!=0 ) {
    if( TMath::Abs(hit->GetTime()-this->GetTime())<win ) {
      if( TMath::Abs(hit->GetZPlane()-this->GetZPlane())<3 && 
	  TMath::Abs(hit->GetYStrip()->GetStrip()-this->GetYStrip()->GetStrip())<4) {
	return 2;
      } else if( TMath::Abs(hit->GetZPlane()-this->GetZPlane())<5 && 
		 TMath::Abs(hit->GetYStrip()->GetStrip()-this->GetYStrip()->GetStrip())<6) {
	return 1;
      }
    }
  }
  return 0;
}
// asmS:2=> closely packed hits than 1 , 0=> no association

int InoHit::IsDiffuseShwAssoc(InoHit* hit) const {
  cout<<" stripwdt:IsDiffShwAssoc " << StripXWidth << " "<< StripYWidth<<endl;
  double win=99.9;
  
  if (fXStrip!=0 && fYStrip!=0 ) {
    if( TMath::Abs(hit->GetTime()-this->GetTime())<win ) {
      if( TMath::Abs(hit->GetZPlane()-this->GetZPlane())<5 && 
	  TMath::Abs(hit->GetXStrip()->GetStrip()-this->GetXStrip()->GetStrip())<11&&
	  TMath::Abs(hit->GetYStrip()->GetStrip()-this->GetYStrip()->GetStrip())<11 ) {
	return 2;
      } else if( TMath::Abs(hit->GetZPlane()-this->GetZPlane())<9 && 
		 TMath::Abs(hit->GetXStrip()->GetStrip()-this->GetXStrip()->GetStrip())<21 &&
		 TMath::Abs(hit->GetYStrip()->GetStrip()-this->GetYStrip()->GetStrip())<21 ) {
	return 1;
      }
    }
  } else if(fXStrip!=0){
    if( TMath::Abs(hit->GetTime()-this->GetTime())<win ) {
      if( TMath::Abs(hit->GetZPlane()-this->GetZPlane())<5 && 
	  TMath::Abs(hit->GetXStrip()->GetStrip()-this->GetXStrip()->GetStrip())<11 ) {
	return 2;
      } else if( TMath::Abs(hit->GetZPlane()-this->GetZPlane())<9 && 
		 TMath::Abs(hit->GetXStrip()->GetStrip()-this->GetXStrip()->GetStrip())<21 ) {
	return 1;
      }
    }
  } else if (fYStrip!=0){
    if( TMath::Abs(hit->GetTime()-this->GetTime())<win ) {
      if( TMath::Abs(hit->GetZPlane()-this->GetZPlane())<5 && 
	  TMath::Abs(hit->GetYStrip()->GetStrip()-this->GetYStrip()->GetStrip())<11 ) {
	return 2;
      } else if( TMath::Abs(hit->GetZPlane()-this->GetZPlane())<9 && 
		 TMath::Abs(hit->GetYStrip()->GetStrip()-this->GetYStrip()->GetStrip())<21) {
	return 1;
      }
    }
  }
  return 0;
}
//asmS: looks for diffuse hit return value mean same as above.

double InoHit::GetXTimeCorr() const {
  // int nInY = (fYStrip->GetId()>>8)&0x7F;
  return DigiToTimeConv*fXTime - (fYstriploc + 0.5)*SignalSpeed - fXtOffset;
}

double InoHit::GetYTimeCorr() const {
  // int nInX = (fXStrip->GetId()>>8)&0x7F;
  return DigiToTimeConv*fYTime - (fXstriploc + 0.5)*SignalSpeed - fYtOffset;
}

int InoHit::isNoiseHit() const {
  if (fView==2) {
    if(fXStrip->isNoiseStrip()==0 && fYStrip->isNoiseStrip()==0) {return 0;}
    else {return 1;}
  } else  if (fView==1) {
    return fYStrip->isNoiseStrip();
  } else if (fView==0) {
    return fXStrip->isNoiseStrip();
  }
}

double InoHit::GetTime() const {
  if (fView==2) {
    // int nInY = (fYStrip->GetId()>>8)&0x7F;
    // int nInX = (fXStrip->GetId()>>8)&0x7F;
    double fXTimeReturn = DigiToTimeConv*fXTime - (fYstriploc + 0.5)*SignalSpeed - fXtOffset;
    double fYTimeReturn = DigiToTimeConv*fYTime - (fXstriploc + 0.5)*SignalSpeed - fYtOffset;
    return 0.5*(fXTimeReturn+fYTimeReturn);
  } else if (fView==1) {
    return fYTime - 4.8; 
  } else if (fView==0) {
    return fXTime - 4.8; 
  }
  return -999.;
}

bool InoHit::isIdentical(InoHit* hit) {

  if ((hit->GetXStripNum() == GetXStripNum()) &&
      (hit->GetYStripNum() == GetYStripNum())) {
    return true;
  }
  return false;
}

int InoHit::GetRPCmod() const {
  if (fView%2==0) return fXStrip->GetRPCmod();
  if (fView>0) return fYStrip->GetRPCmod();
  return 0;
}

