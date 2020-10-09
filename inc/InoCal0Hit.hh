//

#ifndef InoCal0Hit_h
#define InoCal0Hit_h 1

#include "G4VHit.hh"
#include "G4ThreeVector.hh"

class InoCal0Hit : public G4VHit {
public:
  InoCal0Hit();
  ~InoCal0Hit();
  InoCal0Hit(const InoCal0Hit &right);
  const InoCal0Hit& operator=(const InoCal0Hit &right);
  G4int operator==(const InoCal0Hit &right) const;
  
  // inline void *operator new(size_t);
  // inline void operator delete(void *aHit);
  
  void Draw();
  void Print();
  
private:
  G4double edep;
  G4int pdgid;  //Particle ID
  G4ThreeVector pos;
  G4double localx;
  G4double localy;
  G4ThreeVector mom; //Momentum of track at earliest energy deposite, has meaning only for muon track, not usefull at all for hadronic shower
  G4double toff;
  G4double tofx;
  G4double tofy;
  unsigned long HitId;
public:
  inline void SetpdgId(G4int id) { pdgid = id; }
  inline void SetEdep(G4double de) { edep = de; }
  inline void AddEdep(G4double de) { edep +=de; }
  inline G4double GetEdep() { return edep; }
  inline void SetPos(G4ThreeVector xyz) { pos = xyz; }
  inline G4ThreeVector GetPos() { return pos; }
  inline void SetMom(G4ThreeVector xyz) { mom = xyz; }
  inline G4ThreeVector GetMom() { return mom; }
  inline void SetTime(G4double tf) { toff = tf; }
  inline G4double GetTime() { return toff; }
  inline void SetLocalXPos(G4double xyz) { localx = xyz; }
  inline G4double GetLocalXPos() { return localx; }
  inline void SetLocalYPos(G4double xyz) { localy = xyz; }
  inline G4double GetLocalYPos() { return localy; }
  inline void SetHitId (unsigned long id) { HitId = id; }
  inline unsigned long GetHitId() { return HitId; }
  inline G4int GetpdgId() { return pdgid; }
};

#endif
