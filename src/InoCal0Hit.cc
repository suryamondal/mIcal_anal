#include "InoCal0Hit.hh"
using namespace std;

InoCal0Hit::InoCal0Hit() {
  pdgid=-25;
  edep = 0;
  toff = 1000000;
}

InoCal0Hit::~InoCal0Hit() {;}

InoCal0Hit::InoCal0Hit(const InoCal0Hit &right) : G4VHit() {
  pdgid  = right.pdgid;
  edep = right.edep;
  pos = right.pos;
  mom = right.mom;
  toff = right.toff;
  HitId = right.HitId;
}

const InoCal0Hit& InoCal0Hit::operator=(const InoCal0Hit &right) {
  pdgid = right.pdgid;
  edep = right.edep;
  pos = right.pos;
  mom = right.mom;
  toff = right.toff;
  HitId = right.HitId;
  return *this;
}

G4int InoCal0Hit::operator==(const InoCal0Hit &right) const {
  return (this==&right) ? 1 : 0;
}

void InoCal0Hit::Draw() {;}

void InoCal0Hit::Print() {
  cout<<"hit "<<HitId<<" "<<pos<<endl;
}


