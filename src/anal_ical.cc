// Reconstruction Code of DIGI output of INOICAL simulation code

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <map>
#include <iomanip>
#include <utility>
#include "ParameterMessenger.hh"
#include "DetectorParameterDef.hh"
#include "FieldPropagator.hh"
#include "InoDigiAlg.hh"
#include "InoRecoAlg.hh"
#include "GeneralRecoInfo.hh"
#include "MultiSimAnalysisDigi.hh"
#include "vect_manager.h"


int main(int argc, char** argv) {
  
  cout<<"---------------------------------------"<<endl;
  bool runcode = true;

  ParameterMessenger* detectorConfig = new ParameterMessenger();
  int InputOutput = atoi(argv[2]);
  double TimeToDigi = atof(argv[3]);
  int printModulo = atoi(argv[5]);
  int collatedIn;
  if(InputOutput==0) {
    collatedIn = atoi(argv[4]);
  } else {
    collatedIn = 0;
  }
  detectorConfig->SetInputOutput(InputOutput);
  detectorConfig->SetTimeToDigiConv(TimeToDigi);
  detectorConfig->SetCollatedIn(collatedIn);
  DetectorParameterDef* paradef = new DetectorParameterDef;
  InoGeometry_Manager* geoManager;
  GeneralRecoInfo* greco;
  FieldPropagator *pfield;
  MultiSimAnalysisDigi *pAnalysis;

  char fileCorrName[300];
  if(InputOutput) {
    if(InputOutput==1) {
      greco = new GeneralRecoInfo();
    } else if(InputOutput==2) {
      sprintf(fileCorrName,"%s",argv[4]);
      greco = new GeneralRecoInfo(fileCorrName);
    }
    if(paradef->GetDetectorType()==0) {
      geoManager = new InoGeometry_Manager("icalgeo.gdml");
    } else if(paradef->GetDetectorType()==1) {
      geoManager = new InoGeometry_Manager("geo_mical_world.gdml");
    }  
    pfield = new  FieldPropagator();
  }
  
  char rootfiles[200];
  char ffoutname[200];
  
  sprintf(rootfiles, "%s",argv[1]);//q_mical_test_300k_m1.log"); //logFiles/w_test_sim1.log");//mical_sim.log");
  // cout<<"Enter filename: ";
  // cin >> rootfiles;

  int len1 = strlen(rootfiles);
  strncpy(ffoutname,rootfiles,len1-4);
  ffoutname[len1-4] = '\0';

  if(InputOutput) {
    greco->OpenRootFiles(ffoutname);
    pfield->PrintFieldMap();
  }
  
  char outfile[300];
  if(runcode) {
    pAnalysis = new MultiSimAnalysisDigi();
    pAnalysis->SetTimeToDigiConvVal(TimeToDigi);//detectorConfig->GetTimeToDigiConv());
    if(InputOutput) {
      if(InputOutput>1) {
	sprintf(outfile,"./recodata/%s_data",ffoutname);
      } else {
	sprintf(outfile,"./recodata/%s_reco",ffoutname);
      }
    } else {
      sprintf(outfile,"./digidata/%s_digi",ffoutname);
      // sprintf(outfile,"/home/apoorva/mical_20190829/%s_digi",ffoutname);
    }
    cout<<"outfile "<<outfile<<endl;
    pAnalysis->OpenOutputRootFiles(outfile);
    
    UInt_t numEntry_old = 0;
    UInt_t nfileRead = 0;
    ifstream file_db;
    file_db.open(rootfiles);  
    while(!(file_db.eof())) {
      char indatafile[300];
      char outfilx[300];
      char infile[300];
      int nentrymx = 0;
      int ini_ievt = 0; //205474;

      file_db >> indatafile>>nentrymx>>ini_ievt;
      if (strstr(indatafile,"#")) continue;
      if(file_db.eof()) break;
      if(InputOutput==1) {
	sprintf(infile, "./digidata/%s_digi", indatafile);
      } else if(InputOutput==2) {
	sprintf(infile, "../rredata/%s", indatafile);
	cout<<"ii "<<infile<<endl;
      } else {
	// sprintf(infile, "/media/surya/Surya_5/sim01Backup/Gobinda/IICHEP/mical_20190829/%s_sim", indatafile);
	sprintf(infile, "./simdata/%s_sim", indatafile);
      }
      cout<<"infile is "<<indatafile<<endl;
      cout<<"outfile is "<<outfile<<endl;
      pAnalysis->OpenInputRootFiles(infile);
      cout<<"main pAnalysis "<<pAnalysis<<endl;
      if(InputOutput==0) {
	pAnalysis->OpenCollatedRootFile();
      }
      pAnalysis->inputRootFile->cd();
      int numEntry = pAnalysis->inputEventTree->GetEntries();
      cout<<"Input file "<<infile<<" has "<<numEntry<<" entries....."<<endl;
      numEntry = min(numEntry,nentrymx);
      cout<<endl;
      cout<<"ini "<<ini_ievt<<" "<<numEntry<<" "<<numEntry+ini_ievt<<endl;
      // numEntry = 10;
      for(int ievt1=ini_ievt; ievt1<numEntry+ini_ievt; ievt1++) {
	// cout<<"numEntry "<<numEntry<<" "<<numEntry_old<<" "<<nfileRead<<endl;
	if(ievt1%printModulo==0) {
	  cout<<"---> Processing Event # "<<ievt1<<"..."<<endl;
	}
	if(InputOutput) {
	  InoRecoAlg RecoAlgINO(InputOutput);
	  RecoAlgINO.ReadEvent(ievt1);
	  pAnalysis->ievt2 = (numEntry_old*nfileRead) + ievt1;
	  // cout<<"panal "<<pAnalysis->ievt2<<endl;
	  // RecoAlgINO.ReadEvent(403);
	  pAnalysis->SaveGenVisFile();
	  RecoAlgINO.PerformTrackReconstruction();
	} else {
	  InoDigiAlg DigiAlgINO;
	  DigiAlgINO.ReadEvent(ievt1);
	  pAnalysis->ievt2 = (numEntry_old*nfileRead) + ievt1;
	  // cout<<"panal "<<pAnalysis->ievt2<<endl;
	  DigiAlgINO.DigitiseSimData();
	  // DigiAlgINO.NoiseGenLoop();
	  DigiAlgINO.CalculateTrigger();
	  DigiAlgINO.SaveDigiData();
	}
	// cout<<"---> Event # "<<ievt1<<" processed."<<endl;
	// cout<<endl;
      } // for(int ievt1=0;ievt1<numEntry;ievt1++) {
      numEntry_old = numEntry;
      nfileRead++;
      pAnalysis->CloseInputRootFiles();
    } // while(!(file_db.eof())) {
    pAnalysis->CloseOutputRootFiles();
    delete pAnalysis; pAnalysis = 0;
  } else {
    cout<<"No code run..."<<endl;
  }

  if(InputOutput==1) {
    if(pfield) {delete pfield; pfield=0;}
    if(geoManager) {delete geoManager; geoManager = 0;}
    if(greco) {
      greco->CloseRootFiles();
      delete greco; greco =0;
    }
  }
  
  delete paradef; paradef=0;
  delete detectorConfig; detectorConfig=0;
    
  cout<<"Bye.. Bye.."<<endl;
  cout<<endl;
  cout<<endl;
  
  
  return 1;
}

	
