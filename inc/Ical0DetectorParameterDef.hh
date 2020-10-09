#ifndef Ical0DetectorParameterDef_h
#define Ical0DetectorParameterDef_h 1
//#include "G4SIunits.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include "G4ios.hh"
#include <G4String.hh>
// #include "TFile.h"
// #include "TNtuple.h"
// #include <libpq-fe.h>
#include <strings.h>

// *** Box geometry of all volume.
//    Half lenght of each rectangle
//   
//    PARINO : mother volume                           IRON    
//    GAPINO : Gap between two INO module              AIR
//    NUMINO : Number of INO modules                   AIR
//    PARLAY : Each layer    division along z cor    AIR
//    PARMOD : each module   division along x cor    AIR
//    PARCHM : each chamber  division along y cor    IRON
//    PARAIR : inside PARCHM (air and RPC volume)    AIR
//    PARIRN :
//    PARCUP : inside PARAIR (Cu + G10 + RPC + GAS)  
//    PARG10 : inside PARCUP (G10 + RPC + GAS)
//    PARQRZ : inside PARG10 (RPC + GAS)
//    PARGAS : inside PARQRZ (GAS)
//    PARIRNLAY : Top and bottom iron layer
//    STRWD  : Strip width
class Ical0DetectorMessenger;

class Ical0DetectorParameterDef
{

  public:
  
  Ical0DetectorParameterDef();
  ~Ical0DetectorParameterDef(){};

  public:
  static Ical0DetectorParameterDef* AnPointer;
  void UpdateDetectorParameterDef();
  void UpdateDetectorParameterDef1();
  bool SetParameterDB();
  void ReadStripInfoFile();

  void SetParameterLocation(G4String value){parameterLocation = value;}
  G4String GetParameterLocation() {return parameterLocation;}
  
  double GetParworld(int i) { return parworld[i];}
  double GetParino(int i) { return parino[i];}
  double GetGapino()       { return gapino;} 
  int    GetNumino()       { return numino;} 
  double GetParcoilspacerpc(int i) { return parcoilspacerpc[i];}
  double GetParcoilspaceiron(int i) { return parcoilspaceiron[i];}
  double GetParairgap1(int i) { return parairgap1[i];}
  double GetParairgap2(int i) { return parairgap2[i];}
  double GetParairgap3(int i) { return parairgap3[i];}
  double GetParairgap4(int i) { return parairgap4[i];}
  double GetParlay(int i) { return parlay[i];}
  double GetParmod(int i) { return parmod[i];}
  double GetParchm(int i) { return parchm[i];}
  double GetPariron(int i) { return pariron[i];}
  double GetParirmod(int i) { return parirmod[i];}
  double GetParirlay(int i) { return parirlay[i];}
  double GetParspacer1(int i) { return parspacer1[i];}
  double GetParspacer2(int i) { return parspacer2[i];}
  double GetParspacer3(int i) { return parspacer3[i];}
  double GetParspacer4(int i) { return parspacer4[i];}
  double GetParspacer5(int i) { return parspacer5[i];}
  double GetParspacer6(int i) { return parspacer6[i];}
  double GetParFrpBox(int i) { return parfrpbox[i];}
  double GetParAirBox(int i) { return parairbox[i];}
  double GetParG10Trap1(int i) { return parg10Trap1[i];}
  double GetParG10Trap2(int i) { return parg10Trap2[i];}
  // double GetParg10(int i) { return parg10[i];}
  double GetParal(int i) { return paral[i];}
  double GetParALCutBig(int i) { return paralCutBig[i];}
  double GetParALCutSmall(int i) { return paralCutSmall[i];}
  double GetParhoneycomb(int i) { return parhoneycomb[i];}
  double GetParHoneyCombCutBig(int i) { return parhoneycombCutBig[i];}
  double GetParHoneyCombCutSmall(int i) { return parhoneycombCutSmall[i];}
  double GetParcup(int i) { return parcup[i];}
  double GetParCupCutBig(int i) { return parcupCutBig[i];}
  double GetParCupCutSmall(int i) { return parcupCutSmall[i];}
  double GetParmylar(int i) { return parmylar[i];}
  double GetParMylarCutBig(int i) { return parmylarCutBig[i];}
  double GetParMylarCutSmall(int i) { return parmylarCutSmall[i];}
  double GetParcoat(int i) { return parcoat[i];}
  double GetParCoatCutBig(int i) { return parcoatCutBig[i];}
  double GetParCoatCutSmall(int i) { return parcoatCutSmall[i];}
  double GetParqurz(int i) { return parqurz[i];}
  double GetParQurzCutBig(int i) { return parqurzCutBig[i];}
  double GetParQurzCutSmall(int i) { return parqurzCutSmall[i];}
  double GetPargas(int i) { return pargas[i];}
  double GetParGasCutBig(int i) { return pargasCutBig[i];}
  double GetParGasCutSmall(int i) { return pargasCutSmall[i];}
  double GetParvcoil(int i) { return parvcoil[i];}
  double GetParhcoil(int i) { return parhcoil[i];}
  double GetParcurvedcoil(int i) { return parcurvedcoil[i];}
  double GetParcoilsupport(int i) { return parcoilsupport[i];}
  double GetXStrwd() { return Xstrwd;}
  double GetYStrwd() { return Ystrwd;}
  //double GetFeThickness(){ return nFeThickness;}
  //double GetAirGap_btwn_FeLayers(){ return nAirGap;}
  int GetnXStrip() { return nXStrip;} //AAR:
  int GetnYStrip() { return nYStrip;} //AAR:
  int GetnLayer() { return nLayer;}
  int GetnModule() { return nModule;}
  int GetnIron(){return  nIron;}
  int GetnChamber() { return nChamber;} 
  int GetnIRLayer(){return  nIRLayer;}
  int GetnIRModule(){return  nIRModule;}
  int GetnSpacer1(){return nSpacer1;}
  int GetnSpacer2(){return nSpacer2;}
  int GetnSpacer3(){return nSpacer3;}
  int GetnSpacer4(){return nSpacer4;}
  int GetnSpacer5(){return nSpacer5;}
  int GetnSpacer6(){return nSpacer6;}
  int GetnCoil(){return nCoil;}
  int GetnCoilSupport(){return nCoilSupport;}
  unsigned long int StripInfo[3][150][8][8][2][64];
  double GetRPCLayerPosZ(int j){return RPCLayerPosZ[j];}
  double GetIRONLayerPosZ(int j){return IRONLayerPosZ[j];}
  double GetRPCShift(int i) {return RPCShift[i];}

  private : 

  double RPCShift[3];

// *** Defines dimension of USER'S VOLUMES  

//       define MOTHER volume
  double nFeThickness; //5.6 cm Halflength//thickness of iron layer.
  float nAirGap ;     //4 cm Halflength//thickness of Air volume between to Iron layers. >2.9cm

  //# of strips in a chamber :AAR
  int nXStrip;
  int nYStrip;
  // # of layer
  int nLayer; //150 //will be calculated from nFeThickness and nAirGap values set
  // # of chamber
  int nChamber; //8
  // # ofcolumns of  module 
  int nModule; //8
  //# of iron plates in a module
  int nIron; //4
  // # of iron layer
  int nIRLayer; //151 //will be calculated from nLayer
  //# of iron modules in a layer
  int nIRModule; //8
  //# of spacer of type 1 per layer
  int nSpacer1; //4
  //# of spacer of type 2 per layer
  int nSpacer2; //8
  //# of spacer of type 3 per layer
  int nSpacer3; //6
  //# of spacer of type 4 per layer
  int nSpacer4; //21
  //# of spacer of type 5 per layer
  int nSpacer5; //28
  //# of spacer of type 6 per layer
  int nSpacer6; //14
  //# of coils per module
  int nCoil; //4
  //# of supports per coil
  int nCoilSupport; //3
  
  
  // INO detector size   
  float parworld[3];
  float parino[3]; //={1606.01*cm, 706.01*cm, 598.0*cm};
  // Gap between two INO module (full lenght not halflenght like other paameters)
  float gapino; // = 20*cm;
  // NUmber of INO modules
  int   numino; // = 3;
  //space for coil in RPC layer
  float parcoilspacerpc[3];
  //space for coil in IRON layer
  float parcoilspaceiron[3];
  //space for horizontal air gap in IRON layer
  float parairgap1[3];
  //space for vertical air gap in IRON layer
  float parairgap2[3];
  //space for air gap between coils in IRON layer
  float parairgap3[3];
  //space for air gap between coils in IRON layer
  float parairgap4[3];
  // layer
  float parlay[3]; //={1600.0*cm, 700.02*cm, 4.26*cm};
  // module
  float parmod[3]; //={100.01*cm, 700.01*cm, 4.25*cm};
  // chamber size
  float parchm[3]; //={100.0*cm,  100.0*cm, 4.24*cm};
  //iron plate
  float pariron[3];
  //iron module
  float parirmod[3];
  //iron layer
  float parirlay[3];
  //spacer of type 1
  float parspacer1[3];
  //spacer of type 2
  float parspacer2[3];
  //spacer of type 3
  float parspacer3[3];
  //spacer of type 4
  float parspacer4[3];
  //spacer of type 5
  float parspacer5[3];
  //spacer of type 6
  float parspacer6[3];
  

  // //  g10 which includes g10, qrz and gas
  // float parg10[3];//={98.003*cm,  98.003*cm, 0.40*cm};
  //  FRP Tray which includes g10 electronics and RPC
  float parfrpbox[3];
  // Inner Walls of FRP Tray.
  float parairbox[3];
  // Triangle Electronics Board for DFE
  float parg10Trap1[8];
  // Triangle Electronics Board for H.V. supply
  float parg10Trap2[8];

  //Dimensions in (mm)

  //aluminium containing honeycomb
  float paral[3]; //={870.06, 917.56, 9.38}; 
  float paralCutBig[4];//={95, 136.35, 136.35, 1};
  float paralCutSmall[4];//={22.3, 33.537, 33.537, 1};
  
  //honeycomb containing air volume
  float parhoneycomb[3];//={870.05, 917.55, 9.23};
  float parhoneycombCutBig[4];//={95, 136.35, 136.35, 9.23};
  float parhoneycombCutSmall[4];//={22.3, 33.537, 33.537, 9.23};
  
  // cu plate which includes cu, g10, qrz and gas
  float parcup[3];//={870.04, 917.54, 4.23};
  float parcupCutBig[4];//={95, 136.35, 136.35, 4.23};
  float parcupCutSmall[4];//={22.3, 33.537, 33.537, 4.23};
  
  //mylar on graphite
  float parmylar[3];//={870.03, 917.53, 4.13};
  float parmylarCutBig[4];//={95, 136.35, 136.35, 4.13};
  float parmylarCutSmall[4];//={22.3, 33.537, 33.537, 4.13};
  
  // graphite coating on glass
  float parcoat[3];//={870.02, 917.52, 4.03};
  float parcoatCutBig[4];//={95, 136.35, 136.35, 4.03};
  float parcoatCutSmall[4];//={22.3, 33.537, 33.537, 4.03};
  
  // rpc glass which includes qrz and gas
  float parqurz[3];//={870.01, 917.51, 4};
  float parqurzCutBig[4];//={95, 136.35, 136.35, 4};
  float parqurzCutSmall[4];//={22.3, 33.537, 33.537, 4};
  
  // rpc gas module
  float pargas[3];//={870, 917.5, 1};
  float pargasCutBig[4];//={95, 136.35, 136.35, 1};
  float pargasCutSmall[4];//={22.3, 33.537, 33.537, 1};
  
  //coil in vertical direction
  float parvcoil[3];
  
  //coil in horizontal direction
  float parhcoil[3];
  
  //curved portion of coil
  double parcurvedcoil[5];
  
  //coil support
  float parcoilsupport[3];
  
  //  strip width
  float Xstrwd;// = 2.0*cm;
  float Ystrwd;// = 2.0*cm;
  float XYstrwd;//
  

  unsigned long long int strip;
  unsigned long int iRPCmod;
  unsigned int TimeDelay, Ampli, Threshold, Health;

  unsigned long int nInDT, nInLA, nInCH, nInMO, nXY, nStrip;
  
  unsigned long int st1;

  float delta;

  float SpacerWidthEdge; // = 20*mm;
  float SpacerWidthAll; // = 40*mm;
  float Spacer1Length; // = 125*mm;
  float Spacer2Length; // = 250*mm;
  float Spacer3Length; // = 252.5*mm;

  float ChamberSize; // = 1000*mm;

  float CutSmall; // = 919.99*mm;
  float CutBig; // = 919.99*mm;
  float ChamberLength; // = 1000*mm;

  float CoilHLength; // = 3782*mm;
  float CoilVLength; // = 7250*mm;
  float CoilThickness; // = 40*mm;
  float CoilWidth; // = 312.5*mm;
  float CurvedCoilInRadii; // = 178*mm;
  float CurvedCoilOutRadii; // = 258*mm;
  float CurvedCoilPhiMin; // = 0*rad;
  float CurvedCoilPhiMax; // = M_PI/2*rad;
  float CoilSupportLength; // = 100*mm;

  float AirGapIronLayer; // = 2.5*mm;

  //RPC
  float GasChamberSizeX; // = 870*mm; 
  float GasChamberSizeY; // = 917.5*mm;

  float GasGapThickness; // = 1*mm;
  float QurzThickness; // = 3*mm;
  float CoatThickness; // = 0.03*mm;
  float MylarThickness; // = 0.1*mm;
  float CopperThickness; // = 0.1*mm;
  float HoneyCombThickness; // = 5*mm;
  float AluminiumThickness; // = 0.150*mm;

  float BigCutRPC; // = 95*mm;
  float SmallCutRPC; // = 22.3*mm;
  float DeltaChm; // = 0.01*mm
  float DeltaCut; // = 2*mm;
  float FRPBoxXDim; // = 955*mm;
  float FRPBoxYDim; // = 970*mm;
  float FRPBoxZDim; // = 17*mm;
  float FRPThickness; // = 5*mm;
  float FRPZThickness;

  float G10Thickness; // = 2.5*mm;
  float G10TrapDX1; // = 0.01*mm;
  float G10Triangle1Length; // = 125.0*mm;
  float G10Triangle2Length; // = 100.0*mm;
  float G10TrapDY; // = 2.5*mm;

  //To calculate Shift of RPC from FRP Box Center
  float RPCShiftX; // = 119*mm;
  float RPCShiftY; // = 86*mm;
  float RPCShiftZ; // = 9*mm;

  //To calculate pos of G10Trap1 in FRPBox
  float G10ShiftX1; // = 92*mm;
  float G10ShiftY1; // = 60*mm;

  //To calculate pos of G10Trap2 in FRPBox
  float G10ShiftX2; // = 13*mm;
  float G10ShiftY2; // = 1*mm;

  //Ellipsoid World SemiAxis;
  float WorldXX;
  float WorldYY;
  float WorldZZ;

  G4String parameterLocation;
  std::string UserName;
  double RPCLayerPosZ[300];
  double IRONLayerPosZ[300];

};
#endif
