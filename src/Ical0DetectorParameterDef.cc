#include "Ical0DetectorParameterDef.hh"
//#include "G4SIunits.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
#include "ParameterMessenger.hh"

using namespace std;

Ical0DetectorParameterDef *Ical0DetectorParameterDef::AnPointer;


Ical0DetectorParameterDef::Ical0DetectorParameterDef(){
  //       define MOTHER volume
  // INO detector size
  
  AnPointer = this;

  ParameterMessenger* detectorConfig = ParameterMessenger::AnPointer;
  //  cout<<" DatabaseManager* DBmanager = new DatabaseManager();___1"<<endl;
  // DatabaseManager* DBmanager = new DatabaseManager();
  //  cout<<" DatabaseManager* DBmanager = new DatabaseManager();___2"<<endl;
  SetParameterLocation(detectorConfig->GetParameterLocation());


  //GMA
  nFeThickness= 56*mm;
  nFeThickness= detectorConfig->GetnFeThickness();
  nAirGap = 40*mm;
  nAirGap = detectorConfig->GetnAirGap();
  XYstrwd = 30;
  XYstrwd = detectorConfig->GetXYstrwd();
  if(XYstrwd < 15) {
    cout<<"-------- WWWW ------- inoical Exception-START -------- WWWW -------"<<endl;
    cout<<"*** inoical Exception : strip width < 15 mm"<<endl;
    cout<<"issued by : Ical0DetectorParameterDef::Ical0DetectorParameterDef()"<<endl;
    cout<<"*** inoical Solution : strip width = 15"<<endl;
    cout<<"*** This is just a warning message. ***"<<endl;
    cout<<"-------- WWWW -------- inoical Exception-END --------- WWWW -------"<<endl;
    XYstrwd = 15;
  }

  nLayer = 150; //GMA250808 140; // 182; // 140; // 182  //Z-axis
  nLayer = detectorConfig->GetnLayer();
  if(nLayer>256) {
    cout<<"-------- WWWW ------- inoical Exception-START -------- WWWW -------"<<endl;
    cout<<"*** inoical Exception : No. of Layers > 256"<<endl;
    cout<<"issued by : Ical0DetectorParameterDef::Ical0DetectorParameterDef()"<<endl;
    cout<<"*** inoical Solution : No. of Layers = 150"<<endl;
    cout<<"*** This is just a warning message. ***"<<endl;
    cout<<"-------- WWWW -------- inoical Exception-END --------- WWWW -------"<<endl;
    nLayer = 150;
  }
  nChamber = 8; //7; 08/07/08           //no.of RPCs per chamber
  nModule =  8; //16; 08/07/08          //no. of columns of modules
  nIron = 4;                          //no.of iron plates in a module
  nIRModule = 8;                       //no.of iron plate modules
  nIRLayer = nLayer+1;                     //no. of iron layers
  nSpacer1 = 4;                    //no. of spacers of type 1 in a layer
  nSpacer2 = 8;                    //no. of spacers of type 2 in a layer
  nSpacer3 = 6;                    //no. of spacers of type 3 in a layer
  nSpacer4 = 21;                    //no. of spacers of type 4 in a layer
  nSpacer5 = 28;                    //no. of spacers of type 5 in a layer
  nSpacer6 = 14;                    //no. of spacers of type 6 in a layer
  nCoil = 4;                          //no. of coils in a module
  nCoilSupport = 3;                    //no. of supports per coil
  
  // Gap between two INO module (full lenght not halflenght like other paameters)
  gapino = 20*cm;
  // NUmber of INO modules
  numino = 3; //sd=1

  CoilThickness = 40*mm;
  CoilWidth = 312.5*mm;
  CoilHLength = 3782*mm;
  CoilVLength = 7250*mm;
  CurvedCoilInRadii = 178*mm;
  CurvedCoilOutRadii = 258*mm;
  CurvedCoilPhiMin = 0*rad;
  CurvedCoilPhiMax = M_PI/2*rad;
  CoilSupportLength = 100*mm;

  AirGapIronLayer = 2.5*mm;

  SpacerWidthEdge = 20*mm;
  SpacerWidthAll = 40*mm;
  Spacer1Length = 125*mm;
  Spacer2Length = 250*mm;
  Spacer3Length = 252.5*mm;

  ChamberSize = 1000*mm;
  GasChamberSizeX = 870*mm; 
  GasChamberSizeY = 917.5*mm;

  GasGapThickness = 1*mm;
  QurzThickness = 3*mm;
  CoatThickness = 0.03*mm;
  MylarThickness = 0.1*mm;
  CopperThickness = 0.1*mm;
  HoneyCombThickness = 5*mm;
  AluminiumThickness = 0.150*mm;

  BigCutRPC = 95*mm;
  SmallCutRPC = 22.3*mm;
  DeltaChm = 0.01*mm;
  DeltaCut = 2*mm;

  FRPBoxXDim = 955*mm;
  FRPBoxYDim = 970*mm;
  FRPBoxZDim = 17*mm;
  FRPThickness = 5*mm;

  G10Thickness = 2.5*mm;
  G10TrapDX1 = 0.01*mm;
  G10Triangle1Length = 125.0*mm;
  G10Triangle2Length = 100.0*mm;
  G10TrapDY = 2.5*mm;

  //To calculate Shift of RPC from FRP Box Center
  RPCShiftX = 119*mm;
  RPCShiftY = 86*mm;
  RPCShiftZ = 9*mm;

  //To calculate pos of G10Trap1 in FRPBox
  G10ShiftX1 = 92*mm;
  G10ShiftY1 = 60*mm;

  //To calculate pos of G10Trap2 in FRPBox
  G10ShiftX2 = 13*mm;
  G10ShiftY2 = 1*mm;

  WorldXX = 1*km;
  WorldYY = 1*km;
  WorldZZ = 1*km;

  // Add the I/O from the database inoical_db
  // if(GetParameterLocation() == "database") {
  //   if(SetParameterDB()) {
  //     SetParameterLocation("database");
  //   } else {
  //     SetParameterLocation("local");
  //   }
  // }
  

  UpdateDetectorParameterDef();
}

void Ical0DetectorParameterDef::UpdateDetectorParameterDef(){
  // chamber size
  parchm[0] = ChamberSize; //100.0*cm;
  parchm[1] = ChamberSize; //100.0*cm;
  parchm[2] =nAirGap*0.5; //3.24*cm;
  
  //Ellipsoid World Size
  parworld[0] = WorldXX;
  parworld[1] = WorldYY;
  parworld[2] = WorldZZ;
 
  //INO
  parino[0]=parchm[0]*nModule; // 1606.01*cm;
  parino[1]=parchm[1]*nChamber; //706.01*cm;
  parino[2]=(nLayer*nAirGap*0.5)+(nIRLayer*nFeThickness*0.5) + 280*mm; //750.8*cm; //GMA14
   
  // RPC layer
  parlay[0]=parchm[0]*nModule; //800.0*cm; //1600.0*cm;
  parlay[1]=parchm[1]*nChamber ; // 700.02*cm;
  parlay[2]=nAirGap*0.5;// 3.26*cm;
  
  //space for coil in RPC layer
  parcoilspacerpc[0]=CoilThickness; //4*cm;
  parcoilspacerpc[1]=CoilWidth; //31.25*cm;
  parcoilspacerpc[2]= parlay[2];
  
  //RPC module
  parmod[0]=parchm[0];
  parmod[1]=parchm[1]*nChamber; //700.01*cm;
  parmod[2]=nAirGap*0.5; //3.25*cm
  
  // iron plate
  pariron[0]=parchm[0];
  pariron[1]=2*parchm[1]; //700.01*cm;
  pariron[2]=nFeThickness*0.5; //3.25*cm
  
  // iron module
  parirmod[0]=pariron[0];
  parirmod[1]=pariron[1]*nIron; //700.01*cm;
  parirmod[2]=nFeThickness*0.5;
  
  // iron layer
  parirlay[0]=pariron[0]*nIRModule;
  parirlay[1]=pariron[1]*nIron; //700.01*cm;
  parirlay[2]=nFeThickness*0.5;
  
  //space for coil in iron layer
  parcoilspaceiron[0]=CoilThickness; //4*cm;
  parcoilspaceiron[1]=CoilWidth; //31.25*cm;
  parcoilspaceiron[2]= parirlay[2];
  
  //space for horizontal air gap in iron layer
  parairgap1[0] = parirlay[0] - 5*mm;
  parairgap1[1] = AirGapIronLayer; //2.5*mm;
  parairgap1[2] = parirlay[2];
  
  //space for vertical air gap in iron layer
  parairgap2[0] = AirGapIronLayer; //2.5*mm;
  parairgap2[1] = pariron[1] - 2*parairgap1[1];
  parairgap2[2] = parirlay[2];
  
  //space for air gap between coil in iron layer
  parairgap3[0] = AirGapIronLayer; //2.5*mm;
  parairgap3[1] = (0.5*pariron[1] - parcoilspaceiron[1]- parairgap1[1])*0.5;
  parairgap3[2] = parirlay[2];
  
  //space for air gap between coil in iron layer
  parairgap4[0] = AirGapIronLayer; //2.5*mm;
  parairgap4[1] = 2*parairgap3[1];
  parairgap4[2] = parirlay[2];
  
  //spacer1
  parspacer1[0]=SpacerWidthEdge; //2*cm;
  parspacer1[1]=Spacer1Length; //12.5*cm;
  parspacer1[2]=nAirGap*0.5;
  
  //spacer2
  parspacer2[0]=SpacerWidthEdge; //2*cm;
  parspacer2[1]=Spacer2Length; //25*cm;
  parspacer2[2]=nAirGap*0.5;
  
  //spacer3
  parspacer3[0]=SpacerWidthEdge; //2*cm;
  parspacer3[1]=Spacer3Length; //25.25*cm;
  parspacer3[2]=nAirGap*0.5;
  
  //spacer4
  parspacer4[0]=SpacerWidthAll; //4*cm;
  parspacer4[1]=Spacer3Length; //25.25*cm;
  parspacer4[2]=nAirGap*0.5;
  
  //spacer5
  parspacer5[0]=SpacerWidthAll; //4*cm;
  parspacer5[1]=Spacer2Length; //25*cm;
  parspacer5[2]=nAirGap*0.5;
  
  //spacer6
  parspacer6[0]=SpacerWidthAll; //4*cm;
  parspacer6[1]=Spacer1Length; //12.5*cm;
  parspacer6[2]=nAirGap*0.5;
  
  // rpc gas rectangle(2mm)
  pargas[0] = GasChamberSizeX; //870*mm; //919.99*mm;
  pargas[1] = GasChamberSizeY; //917.5*mm;
  pargas[2] = GasGapThickness;

  // rpc gas Big Cut (2mm)
  pargasCutBig[0] = BigCutRPC;
  pargasCutBig[1] = sqrt(2)*pargasCutBig[0] + DeltaCut; //HalfLength Along X
  pargasCutBig[2] = sqrt(2)*pargasCutBig[0] + DeltaCut; //HalfLength Along Y
  pargasCutBig[3] = GasGapThickness; ////HalfLength Along Z

  // rpc gas Small Cut (2mm)
  pargasCutSmall[0] = SmallCutRPC;
  pargasCutSmall[1] = sqrt(2)*pargasCutSmall[0] + DeltaCut; //HalfLength Along X
  pargasCutSmall[2] = sqrt(2)*pargasCutSmall[0] + DeltaCut; //HalfLength Along Y
  pargasCutSmall[3] = GasGapThickness; //HalfLength Along Z

  // rpc glass which includes qrz and gas (3mm)
  parqurz[0] =pargas[0] + DeltaChm;
  parqurz[1] =pargas[1] + DeltaChm;
  parqurz[2] =pargas[2] + QurzThickness; //0.4*cm;

  // qurz Big Cut 1 (3mm)
  parqurzCutBig[0] = BigCutRPC;
  parqurzCutBig[1] = sqrt(2)*parqurzCutBig[0] + DeltaCut; //HalfLength Along X
  parqurzCutBig[2] = sqrt(2)*parqurzCutBig[0] + DeltaCut; //HalfLength Along Y
  parqurzCutBig[3] = pargasCutBig[3] + QurzThickness; //HalfLength Along Z

  // qurz Small Cut 2 (3mm)
  parqurzCutSmall[0] = SmallCutRPC;
  parqurzCutSmall[1] = sqrt(2)*parqurzCutSmall[0] + DeltaCut; //HalfLength Along X
  parqurzCutSmall[2] = sqrt(2)*parqurzCutSmall[0] + DeltaCut; //HalfLength Along Y
  parqurzCutSmall[3] = pargasCutSmall[3] + QurzThickness; //HalfLength Along Z

  //graphite coat on glass (30um)
  parcoat[0] =pargas[0] + DeltaChm;
  parcoat[1] =pargas[1] + DeltaChm;
  parcoat[2] =parqurz[2] + CoatThickness; // 0.4030*cm;

  // graphite Big Cut 1 (30um)
  parcoatCutBig[0] = BigCutRPC;
  parcoatCutBig[1] = sqrt(2)*parcoatCutBig[0] + DeltaCut; //HalfLength Along X
  parcoatCutBig[2] = sqrt(2)*parcoatCutBig[0] + DeltaCut; //HalfLength Along Y
  parcoatCutBig[3] = parqurzCutBig[3] + CoatThickness; //HalfLength Along Z

  // graphite Small Cut 2 (30um)
  parcoatCutSmall[0] = SmallCutRPC;
  parcoatCutSmall[1] = sqrt(2)*parcoatCutSmall[0] + DeltaCut; //HalfLength Along X
  parcoatCutSmall[2] = sqrt(2)*parcoatCutSmall[0] + DeltaCut; //HalfLength Along Y
  parcoatCutSmall[3] = parqurzCutSmall[3] + CoatThickness; //HalfLength Along Z

  //mylar on graphite (100um)
  parmylar[0] =parcoat[0] + DeltaChm;
  parmylar[1] =parcoat[1] + DeltaChm;
  parmylar[2] =parcoat[2] + MylarThickness; // 0.4130*cm;

  // mylar Big Cut 1 (100um)
  parmylarCutBig[0] = BigCutRPC;
  parmylarCutBig[1] = sqrt(2)*parmylarCutBig[0] + DeltaCut; //HalfLength Along X
  parmylarCutBig[2] = sqrt(2)*parmylarCutBig[0] + DeltaCut; //HalfLength Along Y
  parmylarCutBig[3] = parcoatCutBig[3] + MylarThickness; //HalfLength Along Z

  // mylar Small Cut 2 (100um)
  parmylarCutSmall[0] = SmallCutRPC;
  parmylarCutSmall[1] = sqrt(2)*parmylarCutSmall[0] + DeltaCut; //HalfLength Along X
  parmylarCutSmall[2] = sqrt(2)*parmylarCutSmall[0] + DeltaCut; //HalfLength Along Y
  parmylarCutSmall[3] = parcoatCutSmall[3] + MylarThickness; //HalfLength Along Z

  // copper strip (100um)
  parcup[0]=parmylar[0] + DeltaChm;
  parcup[1]=parmylar[1] + DeltaChm;
  parcup[2]=parmylar[2] + CopperThickness; //0.4230*cm;
  
  // copper Big Cut 1 (100um)
  parcupCutBig[0] = BigCutRPC;
  parcupCutBig[1] = sqrt(2)*parcupCutBig[0] + DeltaCut; //HalfLength Along X
  parcupCutBig[2] = sqrt(2)*parcupCutBig[0] + DeltaCut; //HalfLength Along Y
  parcupCutBig[3] = pargasCutBig[3] + CopperThickness; //HalfLength Along Z

  // copper Small Cut 2 (100um)
  parcupCutSmall[0] = SmallCutRPC;
  parcupCutSmall[1] = sqrt(2)*parcupCutSmall[0] + DeltaCut; //HalfLength Along X
  parcupCutSmall[2] = sqrt(2)*parcupCutSmall[0] + DeltaCut; //HalfLength Along Y
  parcupCutSmall[3] = parmylarCutSmall[3] + CopperThickness; //HalfLength Along Z

  // honeycomb containing air volume (5mm)
  parhoneycomb[0]=parcup[0] + DeltaChm;
  parhoneycomb[1]=parcup[1] + DeltaChm;
  parhoneycomb[2]=parcup[2] + HoneyCombThickness; //0.9230*cm;
  
  // honeycomb Big Cut 1 (5mm)
  parhoneycombCutBig[0] = BigCutRPC;
  parhoneycombCutBig[1] = sqrt(2)*parhoneycombCutBig[0] + DeltaCut; //HalfLength Along X
  parhoneycombCutBig[2] = sqrt(2)*parhoneycombCutBig[0] + DeltaCut; //HalfLength Along Y
  parhoneycombCutBig[3] = parcupCutBig[3] + HoneyCombThickness; //HalfLength Along Z

  // honeycomb Small Cut 2 (5mm)
  parhoneycombCutSmall[0] = SmallCutRPC;
  parhoneycombCutSmall[1] = sqrt(2)*parhoneycombCutSmall[0] + DeltaCut; //HalfLength Along X
  parhoneycombCutSmall[2] = sqrt(2)*parhoneycombCutSmall[0] + DeltaCut; //HalfLength Along Y
  parhoneycombCutSmall[3] = parcupCutSmall[3] + HoneyCombThickness; //HalfLength Along Z

  // aluminium containing honeycomb (150um)
  paral[0]=parhoneycomb[0] + DeltaChm;
  paral[1]=parhoneycomb[1] + DeltaChm;
  paral[2]=parhoneycomb[2] + AluminiumThickness; //0.9380*cm;
  
  // aluminium Big Cut 1 (150um)
  paralCutBig[0] = BigCutRPC;
  paralCutBig[1] = sqrt(2)*paralCutBig[0] + DeltaCut; //HalfLength Along X
  paralCutBig[2] = sqrt(2)*paralCutBig[0] + DeltaCut; //HalfLength Along Y
  paralCutBig[3] = parhoneycombCutBig[3] + AluminiumThickness; //HalfLength Along Z

  // aluminium Small Cut 2 (150um)
  paralCutSmall[0] = SmallCutRPC;
  paralCutSmall[1] = sqrt(2)*paralCutSmall[0] + DeltaCut; //HalfLength Along X
  paralCutSmall[2] = sqrt(2)*paralCutSmall[0] + DeltaCut; //HalfLength Along Y
  paralCutSmall[3] = parhoneycombCutSmall[3] + AluminiumThickness; //HalfLength Along Z

  // FRP Tray Dimensions
  parfrpbox[0] = FRPBoxXDim; //955*mm; 
  parfrpbox[1] = FRPBoxYDim; //970*mm; 
  parfrpbox[2] = FRPBoxZDim; //17*mm;

  //  Domensions of InnerWall of FRP Tray. Airbox in which G10 Trapezoid and Aluminium are placed
  parairbox[0] = parfrpbox[0] - FRPThickness;
  parairbox[1] = parfrpbox[1] - FRPThickness;
  parairbox[2] = parfrpbox[2] - FRPThickness + 0.5*mm;

  // The position of RPC w.r.t. the inner walls of FRP Tray. Shift to be added in the data while reconstucting the position of hit.
  // ShiftInX = /*119*mm*/ RPCShiftX - parairbox[0] + paral[0];
  // ShiftInY = /*86*mm*/ RPCShiftY - parairbox[1] + paral[1];
  // ShiftInZ = /*9*mm*/ RPCShiftZ - parfrpbox[2] + parairbox[2];
  RPCShift[0] = /*119*mm*/ RPCShiftX - parairbox[0] + paral[0];
  RPCShift[1] = /*86*mm*/ RPCShiftY - parairbox[1] + paral[1];
  RPCShift[2] = /*9*mm*/ RPCShiftZ - parfrpbox[2] + parairbox[2];

  // //  g10 which includes g10, qrz and gas (5mm)
  // parg10[0] =pargas[0] + 0.001*cm;
  // parg10[1] =pargas[1] + 0.001*cm;
  // parg10[2] =1.4380*cm;

  //g10 trapzoid 1 (Triangle board for DFE)
  parg10Trap1[0] = G10Triangle1Length; //125*mm; //Side Length of Triangle
  parg10Trap1[1] = G10TrapDX1; //0.01*mm; //(Half-length along x at the surface positioned at -dz)
  parg10Trap1[2] = sqrt(2)*parg10Trap1[0]/2; //Hypotenuse of triangle (Half-length along x at the surface positioned at +dz)
  parg10Trap1[3] = G10Thickness; //2.5*mm; // Thickness of the board (Half-length along y)
  parg10Trap1[4] = (sqrt((parg10Trap1[0]*parg10Trap1[0])-(parg10Trap1[0]*parg10Trap1[0]/2.))/2.); //(Half-length along z)
  parg10Trap1[5] = /*92*mm*/ G10ShiftX1 - parairbox[0] + parg10Trap1[0]/2; // xpos wrt center of airbox
  parg10Trap1[6] = /*60*mm*/ G10ShiftY1 - parairbox[1] + parg10Trap1[0]/4; // ypos wrt center of airbox
  parg10Trap1[7] = -parairbox[2] + parg10Trap1[3]; // zpos wrt center of airbox
  
  //g10 trapzoid 2 (Triangle board for HV power supply)
  parg10Trap2[0] = G10Triangle2Length; //100*mm; //Side Length of Triangle
  parg10Trap2[1] = G10TrapDX1; //0.01*mm; //(Half-length along x at the surface positioned at -dz)
  parg10Trap2[2] = sqrt(2)*parg10Trap2[0]/2; //Hypotenuse of triangle (Half-length along x at the surface positioned at +dz)
  parg10Trap2[3] = G10Thickness; //2.5*mm; // Thickness of the board (Half-length along y)
  parg10Trap2[4] = (sqrt((parg10Trap2[0]*parg10Trap2[0])-(parg10Trap2[0]*parg10Trap2[0]/2.))/2.); //(Half-length along z)
  parg10Trap2[5] = parairbox[0] - parg10Trap2[0]/2 - G10ShiftX2; //13*mm; // xpos wrt center of airbox
  parg10Trap2[6] = parairbox[0] - parg10Trap2[0]/4 - G10ShiftY2; //1*mm; // ypos wrt center of airbox
  parg10Trap2[7] = -parairbox[2] + parg10Trap2[3]; // zpos wrt center of airbox
  
  //coil in vertical direction
  parvcoil[0]=CoilThickness; //4*cm;
  parvcoil[1]=CoilWidth; //31.25*cm;
  parvcoil[2]=CoilVLength; //725*cm ;
  
  //coil in horizontal direction
  parhcoil[0]=CoilHLength; //378.2*cm;
  parhcoil[1]=CoilWidth; //31.25*cm;
  parhcoil[2]=CoilThickness; //4*cm ;
  
  //curved portion of coil
  parcurvedcoil[0]= CurvedCoilInRadii; //17.8*cm;
  parcurvedcoil[1]= CurvedCoilOutRadii; //25.8*cm;
  parcurvedcoil[2]= CoilWidth; //31.25*cm ;
  parcurvedcoil[3]= CurvedCoilPhiMin; //0*rad ;
  parcurvedcoil[4]= CurvedCoilPhiMax; //M_PI/2*rad ;
  
  //coil support
  parcoilsupport[0]=CoilSupportLength; //10*cm;
  parcoilsupport[1]=CoilWidth; //31.25*cm;
  parcoilsupport[2]=CoilSupportLength; //10*cm ;
  
  //  strip width
  Xstrwd = XYstrwd*mm; //100 strips
  Ystrwd = XYstrwd*mm; // 100 strips
  nXStrip =  int(1.999*pargas[0]/Xstrwd)+1; //parchm[0]*2/Xstrwd; //GMA14
  nYStrip =  int(1.999*pargas[1]/Ystrwd)+1; //parchm[1]*2/Xstrwd;
  // cout<<"nXStrip = "<<nXStrip<<", nYStrip = "<<nYStrip<<endl;
  // ReadStripInfoFile();

  for(int ij=0; ij<nLayer; ij++) {
    RPCLayerPosZ[ij] = -((int(0.5*nIRLayer)-ij)*2*(parirlay[2]+ parlay[2]))+(parirlay[2]+ parlay[2]);
  }

  for(int ij=0; ij<nIRLayer; ij++) {
    IRONLayerPosZ[ij] = -( (int(0.5*nIRLayer)-ij)*2*(parirlay[2]+ parlay[2]) );
  }

   for(int ij=0; ij<10; ij++) {
    cout<<"RPC Layer["<<ij<<"] = "<<RPCLayerPosZ[ij]<<endl;
  }
  for(int ij=0; ij<11; ij++) {
    cout<<"IRON Layer["<<ij<<"] = "<<IRONLayerPosZ[ij]<<endl;
  }
  
}



