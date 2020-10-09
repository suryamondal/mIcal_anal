#ifndef PARAMETERMESSENGER_H
#define PARAMETERMESSENGER_H

#include "globals.hh"

class ParameterMessenger {
public:
  ParameterMessenger();
  ParameterMessenger(char* inf);
  virtual ~ParameterMessenger();
  
  static ParameterMessenger* AnPointer;
  
  //Getters
  G4String GetParameterLocation() {return parameterLocation;}
  G4String GetGeometryLocation() {return geometryLocation;}
  G4String GetStripInfoLocation() {return stripInfoLocation;}
  double GetnFeThickness() {return nFeThickness;}
  double GetnAirGap() {return nAirGap;}
  int GetnLayer() {return nLayer;}
  double GetXYstrwd() {return XYstrwd;}
  G4String GetFileVersion() {return fileVersion;}
  int GetDetectorType() {return DetectorType;}
  int GetInputOutput() {return InputOutput;}
  int GetVisualOutput() {return VisualOutput;}
  int GetXTermOutput() {return XTermOutput;}
  double GetCorrTimeSmr() {return CorrTimeSmr;}
  double GetUnCorrTimeSmr() {return UnCorrTimeSmr;}
  double GetSignalSpeed() {return SignalSpeed;}
  double GetTimeToDigiConv() {return TimeToDigiConv;}
  int GetCollatedIn() {return CollatedIn;}
  
  //Setters
  void SetParameterLocation(G4String value) {parameterLocation = value;}
  void SetGeometryLocation(G4String value) {geometryLocation = value;}
  void SetStripInfoLocation(G4String value) {stripInfoLocation = value;}
  void SetnFeThickness(double value) {nFeThickness = value;}
  void SetnAirGap(double value) {nAirGap = value;}
  void SetnLayer(int value) {nLayer = value;}
  void SetXYstrwd(double value) {XYstrwd = value;}
  void SetFileVersion(G4String value) {fileVersion = value;}
  void SetDetectorType(int value) {DetectorType = value;}
  void SetVisualOutput(int value) {VisualOutput = value;}
  void SetInputOutput(int value) {InputOutput = value;}
  void SetXTermOutput(int value) {XTermOutput = value;}
  void SetCorrTimeSmr(double value) {CorrTimeSmr = value;}
  void SetUnCorrTimeSmr(double value){UnCorrTimeSmr = value;} 
  void SetSignalSpeed(double value){SignalSpeed = value;} 
  void SetTimeToDigiConv(double value){TimeToDigiConv = value;}
  void SetCollatedIn(int value){CollatedIn = value;}
protected:

private:
  G4String parameterLocation;
  G4String geometryLocation;
  G4String stripInfoLocation;
  G4String fileVersion;
  int nLayer;
  double nFeThickness;
  double nAirGap;
  double XYstrwd;
  int DetectorType;
  int InputOutput;
  int VisualOutput;
  int XTermOutput;
  int CollatedIn;
  double CorrTimeSmr;
  double UnCorrTimeSmr;
  double SignalSpeed;
  double TimeToDigiConv;

};

#endif // ICAL0PARAMETERMESSENGER_H
