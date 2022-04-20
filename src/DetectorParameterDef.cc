#include "DetectorParameterDef.hh"
//#include "G4SIunits.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
#include "ParameterMessenger.hh"
using namespace std;
DetectorParameterDef *DetectorParameterDef::AnPointer;

DetectorParameterDef::DetectorParameterDef() {
  
  AnPointer = this;

  ParameterMessenger* detectorConfig = ParameterMessenger::AnPointer;

  DetectorType = detectorConfig->GetDetectorType();
  cout<<" DetectorType "<<DetectorType<<endl;
  if(DetectorType==0) {
    icalparadef = new Ical0DetectorParameterDef();
  } else {
    miniparadef = new micalDetectorParameterDef();
  }
  
  UpdateDetectorPars();
  PrintDetectorParameters();
}

DetectorParameterDef::~DetectorParameterDef() {
  if(icalparadef) {delete icalparadef; icalparadef=0;}
  if(miniparadef) {delete miniparadef; miniparadef=0;}
}

void DetectorParameterDef::UpdateDetectorPars() {
  if(DetectorType==0) {
    for(int ixi=0;ixi<3;ixi++) {
      parworld[ixi] = icalparadef->GetParworld(ixi);
      parino[ixi] = icalparadef->GetParino(ixi);
      parirlay[ixi] = icalparadef->GetParirlay(ixi);
      parlay[ixi] = icalparadef->GetParlay(ixi);
      pargas[ixi] = icalparadef->GetPargas(ixi);
      parchm[ixi] = icalparadef->GetParchm(ixi);
      parmod[ixi] = icalparadef->GetParmod(ixi);
      parhcoil[ixi] = icalparadef->GetParhcoil(ixi);	
      parcoilspacerpc[ixi] = icalparadef->GetParcoilspacerpc(ixi);
      RPCShift[ixi] = icalparadef->GetRPCShift(ixi);
      StackShift[ixi] = 0.0;
    }
    nLayer = icalparadef->GetnLayer();
    nIRLayer = icalparadef->GetnIRLayer();
    nModule = icalparadef->GetnModule();
    nChamber = icalparadef->GetnChamber();
    for(int iyi=0; iyi<nIRLayer; iyi++) {
      IRONLayerPosZ[iyi] = icalparadef->GetIRONLayerPosZ(iyi);
      if(iyi<nIRLayer-1) {
	RPCLayerPosZ[iyi] = icalparadef->GetRPCLayerPosZ(iyi);
	ShiftInZ[iyi] = icalparadef->GetRPCShift(iyi);
	LayerZdim[iyi] = icalparadef->GetParlay(2)*2;
	IronLayerZdim[iyi] = icalparadef->GetParirlay(2)*2;
      }
    }
    gapino = icalparadef->GetGapino();
    Xstrwd = icalparadef->GetXStrwd();
    Ystrwd = icalparadef->GetYStrwd();
    numino = icalparadef->GetNumino();
    nXStrip = icalparadef->GetnXStrip();
    nYStrip = icalparadef->GetnXStrip();
  } else {
    for(int ixi=0;ixi<3;ixi++) {
      parworld[ixi] = miniparadef->GetParworld(ixi);
      parino[ixi] = miniparadef->GetParino(ixi);
      parirlay[ixi] = miniparadef->GetParirlay(ixi);
      parlay[ixi] = miniparadef->GetParlay(ixi);
      pargas[ixi] = miniparadef->GetPargas(ixi);
      parhcoil[ixi] = miniparadef->GetParhcoil(ixi);
      parcoilspacerpc[ixi] = miniparadef->GetParcoilspacerpc(ixi);
      RPCShift[ixi] = miniparadef->GetRPCShift(ixi);
      StackShift[ixi] = miniparadef->GetStackPosInRoom(ixi) + miniparadef->GetINOroomPos(ixi);
      // StackPosInRoom[ixi] = miniparadef->GetStackPosInRoom(ixi);
      // INOroomPos[ixi] = miniparadef->GetINOroomPos(ixi);
    }
    nLayer = miniparadef->GetnLayer();
    nIRLayer = miniparadef->GetnIRLayer();
    nModule = 1;
    nChamber = 1;
    for(int iyi=0; iyi<nIRLayer; iyi++) {
      IRONLayerPosZ[iyi] = miniparadef->GetIRONLayerPosZ(iyi);
      IronLayerZdim[iyi] = miniparadef->GetIronLayerZdim(iyi);
      if(iyi<nIRLayer-1) {
	RPCLayerPosZ[iyi] = miniparadef->GetRPCLayerPosZ(iyi);
	ShiftInZ[iyi] = miniparadef->GetShiftInZ(iyi);
	LayerZdim[iyi] = miniparadef->GetLayerZdim(iyi);
      }
    }
    gapino = 0;
    Xstrwd = miniparadef->GetXStrwd();
    Ystrwd = miniparadef->GetYStrwd();
    numino = 1;
    nXStrip = miniparadef->GetnXStrip();
    nYStrip = miniparadef->GetnYStrip();
    
  }
  
}

void DetectorParameterDef::PrintDetectorParameters() {
  cout << "------------------------------------------------------------"<<endl;
  cout <<"Detector Parameters"<<endl; 
  cout <<"parworld "<<parworld[0]<<"*mm, "<<parworld[1]<<"*mm, "<<parworld[2]<<"*mm"<<endl;
  cout <<"parino "<<parino[0]<<"*mm, "<<parino[1]<<"*mm, "<<parino[2]<<"*mm"<<endl;
  cout <<"parlay "<<parlay[0]<<"*mm, "<<parlay[1]<<"*mm, "<<parlay[2]<<"*mm"<<endl;
  cout <<"parirlay "<<parirlay[0]<<"*mm, "<<parirlay[1]<<"*mm, "<<parirlay[2]<<"*mm"<<endl;
  cout <<"pargas "<<pargas[0]<<"*mm, "<<pargas[1]<<"*mm, "<<pargas[2]<<"*mm"<<endl;
  cout <<"parmod "<<parmod[0]<<"*mm, "<<parmod[1]<<"*mm, "<<parmod[2]<<"*mm"<<endl;
  cout <<"parchm "<<parchm[0]<<"*mm, "<<parchm[1]<<"*mm, "<<parchm[2]<<"*mm"<<endl;
  cout <<"parhcoil "<<parhcoil[0]<<"*mm, "<<parhcoil[1]<<"*mm, "<<parhcoil[2]<<"*mm"<<endl;
  cout <<"parcoilsupport "<<parcoilsupport[0]<<"*mm, "<<parcoilsupport[1]<<"*mm, "<<parcoilsupport[2]<<"*mm"<<endl;
  cout <<"parcoilspacerpc "<<parcoilspacerpc[0]<<"*mm, "<<parcoilspacerpc[1]<<"*mm, "<<parcoilspacerpc[2]<<"*mm"<<endl;
  cout<<"RPC shift "<<RPCShift[0]<<" "<<RPCShift[1]<<" "<<RPCShift[2]<<endl;
  cout<<"Stack shift "<<StackShift[0]<<" "<<StackShift[1]<<" "<<StackShift[2]<<endl;
  cout<<"No. of Strips (X:Y) ("<<nXStrip<<","<<nYStrip<<")"<<endl;
  cout << "------------------------------------------------------------"<<endl;
 
}
