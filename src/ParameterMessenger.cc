#include "ParameterMessenger.hh"
#include <libconfig.h++>

using namespace libconfig;
using namespace std;

ParameterMessenger *ParameterMessenger::AnPointer;

ParameterMessenger::ParameterMessenger() {
  AnPointer = this;
  SetDetectorType(1);
  SetInputOutput(0);
  SetCollatedIn(1);
  SetVisualOutput(0);
  SetXTermOutput(0);
  SetParameterLocation("local");
  SetGeometryLocation("local");
  SetStripInfoLocation("local");
  SetnFeThickness(56.0);
  SetXYstrwd(30.0);
  if(DetectorType==0) {
    SetnLayer(150);
    SetnAirGap(40.0);
  } else {
    SetnLayer(10);
    SetnAirGap(45.0);
  }
  SetCorrTimeSmr(0.0);
  SetUnCorrTimeSmr(0.0);
  SetSignalSpeed(0.15);
  SetTimeToDigiConv(0.1);
}  
  
ParameterMessenger::ParameterMessenger(char* infile)
{
  AnPointer = this;
  //ctor
  Config cfg;
  // cout<<"===================================================="<<endl;
  // Read the config file. If there is an error,
  // report it and exit.
  try {
    // cfg.readFile("/products/GEANT4.10/ICALDOCS/detector_config.cfg");
    cfg.readFile(infile);
    // cout <<"In file detector_config.cfg =>"<<endl;
  } catch
    (const FileIOException& fioex) {
    G4cerr << "I/O error while reading file." << G4endl;
    G4cerr << "Setting the parameters to default values ...."
           << G4endl;
    SetDetectorType(1);
    SetInputOutput(0);
    SetVisualOutput(0);
    SetXTermOutput(0);
    SetParameterLocation("local");
    SetGeometryLocation("local");
    SetStripInfoLocation("local");
    SetXYstrwd(30.0);
    SetnFeThickness(56.0);
    SetnAirGap(45.0);
    SetnLayer(10);
    SetCorrTimeSmr(0.7);
    SetUnCorrTimeSmr(0.7);
    SetSignalSpeed(0.15);
    SetTimeToDigiConv(0.1);
    SetCollatedIn(0);
  }

  //Setup reading of the detector parameters
  int DetectorTypePar;
  bool cfgDetectorType = cfg.lookupValue("detector_type",DetectorTypePar);
  if(cfgDetectorType) {
    cout<<"cfgDetectorType "<<cfgDetectorType<<endl;
    SetDetectorType(DetectorTypePar);
    cout<<"Setting detector type = "<<DetectorTypePar<<" ... "<<endl;
  } else {
    SetDetectorType(0);
    cout << "No 'detector_type' setting in configuration file."
           << endl;
    cout << "Setting the detector to default value 0 ...." << endl;
  }

  int IsInOut;
  bool cfgIsInOut = cfg.lookupValue("input_output",IsInOut);
  if(cfgIsInOut) {
    cout<<"cfgIsInOut "<<cfgIsInOut<<endl;
    SetVisualOutput(IsInOut);
    cout<<"Setting visual output = "<<IsInOut<<" ... "<<endl;
  } else {
    SetInputOutput(0);
    cout << "No 'input_output' setting in configuration file."
           << endl;
    cout << "Setting the input/output to default value 0 ...." << endl;
  }
  
  int IsVisOut;
  bool cfgIsVisOut = cfg.lookupValue("visual_output",IsVisOut);
  if(cfgIsVisOut) {
    cout<<"cfgIsVisOut "<<cfgIsVisOut<<endl;
    SetVisualOutput(IsVisOut);
    cout<<"Setting visual output = "<<IsVisOut<<" ... "<<endl;
  } else {
    SetVisualOutput(0);
    cout << "No 'visual_output' setting in configuration file."
           << endl;
    cout << "Setting the visual output to default value 0 ...." << endl;
  }

  int IsCollatedIn;
  bool cfgIsCollatedIn = cfg.lookupValue("collated_input",IsCollatedIn);
  if(cfgIsCollatedIn) {
    cout<<"cfgIsCollatedIn "<<cfgIsCollatedIn<<endl;
    SetCollatedIn(IsCollatedIn);
    cout<<"Setting collated input = "<<IsCollatedIn<<" ... "<<endl;
  } else {
    SetCollatedIn(0);
    cout << "No 'collated_input' setting in configuration file."
           << endl;
    cout << "Setting the collated input to default value 0 ...." << endl;
  }
    
  int IsXTermOut;
  bool cfgIsXTermOut = cfg.lookupValue("xterm_output",IsXTermOut);
  if(cfgIsXTermOut) {
    cout<<"cfgIsXTermOut "<<cfgIsXTermOut<<endl;
    SetXTermOutput(IsXTermOut);
    cout<<"Setting visual output = "<<IsXTermOut<<" ... "<<endl;
  } else {
    SetXTermOutput(0);
    cout << "No 'xterm_output' setting in configuration file."
           << endl;
    cout << "Setting the xterm output to default value 0 ...." << endl;
  }
  
  // Get the parameter location.
  try {
    string location = cfg.lookup("parameter_location");
    cout <<"parameter location = "<< location<<endl;
    SetParameterLocation(location);
  } catch
    (const SettingNotFoundException& nfex) {
    SetParameterLocation("local");
    G4cerr << "No 'parameter_location' setting in configuration file."
           << G4endl;
    G4cerr << "Setting the parameters to default values ...."
           << G4endl;
  }

  // Get the geometry location
  try {
    string location = cfg.lookup("geometry_location");
    //    cout <<"geometry location = "<< location<<endl;
    SetGeometryLocation(location);
    cout<<"Setting GeomoetryLocation = "<<location<<" ... "<<endl;
  } catch
    (const SettingNotFoundException& nfex) {
    SetGeometryLocation("local");
    cout<<"Setting GeometryLocation = local ... "<<endl;
    G4cerr << "No 'geometry_location' setting in configuration file."
           << G4endl;
    G4cerr << "Setting the geometry to default value ...."
           << G4endl;
  }

  // Get the StripInfo location
  try {
    string location = cfg.lookup("stripInfo_location");
    //    cout <<"stripInfo location = "<< location<<endl;
    SetStripInfoLocation(location);
    cout<<"Setting StripInfoLocation = "<<location<<" ... "<<endl;
  } catch
    (const SettingNotFoundException& nfex) {
    cout<<"Setting StripInfoLocation = local ... "<<endl;
    SetStripInfoLocation("local");
    G4cerr << "No 'SetStripInfoLocation' setting in configuration file."
           << G4endl;
    G4cerr << "Setting the geometry to default value ...."
           << G4endl;
  }

  double ironThickness;
  bool cfgiron = cfg.lookupValue("iron_thickness",ironThickness);
  if(cfgiron) {
    cout<<"cfgiron "<<cfgiron<<endl;
    SetnFeThickness(ironThickness);
    cout<<"Setting IronThickness = "<<ironThickness<<" mm ... "<<endl;
  } else {
    SetnFeThickness(56.0);
    cout<<"cfgiron "<<cfgiron<<endl;
    cout << "No 'iron_thickness' setting in configuration file."
           << endl;
    cout << "Setting the iron thickness to default value 56 mm ...."
           << endl;
  }

  double airThickness;
  bool cfgair = cfg.lookupValue("air_thickness",airThickness);
  if(cfgair) {
    cout<<"cfgair "<<cfgair<<endl;
    SetnAirGap(airThickness);
    cout<<"Setting AirThickness = "<<airThickness<<" mm ... "<<endl;
  } else {
    SetnAirGap(40.0);
    cout<<"cfgair "<<cfgair<<endl;
    cout << "No 'air_thickness' setting in configuration file." << endl;
    cout << "Setting the air thickness to default value 40 mm ...." << endl;
  }

  int layern;
  bool cfglayern = cfg.lookupValue("n_layer",layern);
  if(cfglayern) {
    cout<<"cfglayern "<<cfglayern<<endl;
    SetnLayer(layern);
    cout<<"Setting number of layers = "<<layern<<" ... "<<endl;
  } else {
    SetnLayer(150);
    cout << "No 'n_layer' setting in configuration file."
           << endl;
    cout << "Setting the number of layers to default value 150 ...." << endl;

  }

  double strwdXY;
  bool cfgstrwdXY = cfg.lookupValue("strwd_xy",strwdXY);
  if(cfgstrwdXY) {
    cout<<"cfglayern "<<cfglayern<<endl;
    SetXYstrwd(strwdXY);
    cout<<"Setting X-Y stripwidth = "<<strwdXY<<" mm ... "<<endl;
  } else {
    SetXYstrwd(30.0);
    cout << "No 'strwd_xy' setting in configuration file." << endl;
    cout << "Setting the stripwidth to default value 30 mm ...." << endl;
  }

  double CorrTimeSmrPar;
  bool cfgCorrTimeSmr = cfg.lookupValue("corr_time_smr",CorrTimeSmrPar);
  if(cfgCorrTimeSmr) {
    cout<<"cfgCorrTimeSmr "<<cfgCorrTimeSmr<<endl;
    SetCorrTimeSmr(CorrTimeSmrPar);
    cout<<"Setting CorrTimeSmr = "<<CorrTimeSmrPar<<" ns ... "<<endl;
  } else {
    SetCorrTimeSmr(0.7);
    cout<<"cfgCorrTimeSmr "<<cfgCorrTimeSmr<<endl;
    cout << "No 'corr_time_smr' setting in configuration file." << endl;
    cout << "Setting the corr time smear to default value 0.7 ns ...." << endl;
  }

  double UnCorrTimeSmrPar;
  bool cfgUnCorrTimeSmr = cfg.lookupValue("uncorr_time_smr",UnCorrTimeSmrPar);
  if(cfgUnCorrTimeSmr) {
    cout<<"cfgUnCorrTimeSmr "<<cfgUnCorrTimeSmr<<endl;
    SetUnCorrTimeSmr(UnCorrTimeSmrPar);
    cout<<"Setting UnCorrTimeSmr = "<<UnCorrTimeSmrPar<<" ns ... "<<endl;
  } else {
    SetUnCorrTimeSmr(0.7);
    cout<<"cfgXUnCorrTimeSmr "<<cfgUnCorrTimeSmr<<endl;
    cout << "No 'uncorr_time_smr' setting in configuration file." << endl;
    cout << "Setting the uncorr time smear to default value 0.7 ns ...." << endl;
  }

  double SignalSpeedPar;
  bool cfgSignalSpeed = cfg.lookupValue("signal_speed",SignalSpeedPar);
  if(cfgSignalSpeed) {
    cout<<"cfgSignalSpeed "<<cfgSignalSpeed<<endl;
    SetSignalSpeed(SignalSpeedPar);
    cout<<"Setting SignalSpeedPar = "<<SignalSpeedPar<<" ns ... "<<endl;
  } else {
    SetSignalSpeed(0.15);
    cout<<"cfgSignalSpeed "<<cfgSignalSpeed<<endl;
    cout << "No 'signal_speed' setting in configuration file." << endl;
    cout << "Setting the signal speed to default value 0.7 ns ...." << endl;
  }

  double TimeToDigiConvPar;
  bool cfgTimeToDigiConv = cfg.lookupValue("time_to_digi_conv",TimeToDigiConvPar);
  if(cfgTimeToDigiConv) {
    cout<<"cfgTimeToDigiConv "<<cfgTimeToDigiConv<<endl;
    SetTimeToDigiConv(TimeToDigiConvPar);
    cout<<"Setting TimeToDigiConv = "<<TimeToDigiConvPar<<" ns ... "<<endl;
  } else {
    SetTimeToDigiConv(0.1);
    cout<<"cfgTimeToDigiConv "<<cfgTimeToDigiConv<<endl;
    cout << "No 'time_to_digi_conv' setting in configuration file." << endl;
    cout << "Setting the time to digi conv factor to default value 0.2 ns ...." << endl;
  }

//   // Get the file version
//   try {
//     string location = cfg.lookup("file_version");
//     //    cout <<"file_version = "<< location<<endl;
//     SetFileVersion(location);
//   } catch
//     (const SettingNotFoundException& nfex) {
//     SetFileVersion("v1.00");
//     G4cerr << "No 'file_version' setting in configuration file."
//            << G4endl;
//     G4cerr << "Setting the file version to v1.00 ...."
//            << G4endl;
//   }
}

ParameterMessenger::~ParameterMessenger()
{
  //dtor
}
