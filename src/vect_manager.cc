#include "vect_manager.h"

InoTDCHitx_Manager ::InoTDCHitx_Manager() {    APointer = this;  }
InoTDCHitx_Manager:: ~InoTDCHitx_Manager(){ 
  for(int ij=0; ij<nlmx; ij++) {
    for(int jk=0; jk<ntmx; jk++) {
      xtdctiming[ij][jk].clear();
    }
  }
}
InoTDCHitx_Manager* InoTDCHitx_Manager::APointer;

InoTDCHity_Manager ::InoTDCHity_Manager() {    APointer = this;  }
InoTDCHity_Manager:: ~InoTDCHity_Manager(){ 
  for(int ij=0; ij<nlmx; ij++) {
    for(int jk=0; jk<ntmx; jk++) {
      ytdctiming[ij][jk].clear();
    }
  }
}
InoTDCHity_Manager* InoTDCHity_Manager::APointer;

InoCal0Hit_Manager ::InoCal0Hit_Manager() {    APointer = this;  }
InoCal0Hit_Manager:: ~InoCal0Hit_Manager(){}
InoCal0Hit_Manager* InoCal0Hit_Manager::APointer;

InoStrip_Manager ::InoStrip_Manager() {    APointer = this;  }
InoStrip_Manager:: ~InoStrip_Manager(){}
InoStrip_Manager* InoStrip_Manager::APointer;

InoStripX_Manager ::InoStripX_Manager() {    APointer = this;  }
InoStripX_Manager:: ~InoStripX_Manager(){}
InoStripX_Manager* InoStripX_Manager::APointer;

InoStripY_Manager ::InoStripY_Manager() {    APointer = this;  }
InoStripY_Manager:: ~InoStripY_Manager(){}
InoStripY_Manager* InoStripY_Manager::APointer;

InoHit_Manager ::InoHit_Manager() {    APointer = this;  }
InoHit_Manager:: ~InoHit_Manager(){}
InoHit_Manager *InoHit_Manager::APointer;

InoCluster_Manager ::InoCluster_Manager() {    APointer = this;  }
InoCluster_Manager:: ~InoCluster_Manager(){}
InoCluster_Manager *InoCluster_Manager::APointer;

InoTrack_Manager ::InoTrack_Manager() {    APointer = this;  }
InoTrack_Manager:: ~InoTrack_Manager(){}
InoTrack_Manager *InoTrack_Manager::APointer;

InoFittedTrack_Manager ::InoFittedTrack_Manager() {    APointer = this;  }
InoFittedTrack_Manager:: ~InoFittedTrack_Manager(){}
InoFittedTrack_Manager *InoFittedTrack_Manager::APointer;

InoTrackCand_Manager ::InoTrackCand_Manager() {    APointer = this;  }
InoTrackCand_Manager:: ~InoTrackCand_Manager(){}
InoTrackCand_Manager *InoTrackCand_Manager::APointer;

InoGeometry_Manager ::InoGeometry_Manager(G4String geoFiles) 
{    
  cout <<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx "<<endl;
  icalGeometry        = TGeoManager::Import(geoFiles);
  cout << "Imported geometry file " << geoFiles << " ..." << endl;
  APointer = this;
}
InoGeometry_Manager:: ~InoGeometry_Manager(){}
InoGeometry_Manager *InoGeometry_Manager::APointer;
/*
InoMuRange_Manager :: InoMuRange_Manager(){APointer = this;}
InoMuRange_Manager ::~InoMuRange_Manager(){}
InoMuRange_Manager *InoMuRange_Manager::APointer;
*/


InoRPCStrip_Manager :: InoRPCStrip_Manager(){APointer = this;}
InoRPCStrip_Manager ::~InoRPCStrip_Manager(){}
InoRPCStrip_Manager *InoRPCStrip_Manager::APointer;
