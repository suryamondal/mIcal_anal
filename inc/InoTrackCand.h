#ifndef INOTRACKCAND_H
#define INOTRACKCAND_H
//CandTrackHandle
#include "TObject.h"
#include <map>
class InoTrack;
class InoShowerCand;
class InoVertex;
class InoHit;
class InoTrackCand //: public CandRecoHandle
{

public:
  InoTrackCand();
  InoTrackCand(const InoTrackCand &cdh);
  InoTrackCand(InoTrack *trk, bool forward);
  virtual ~InoTrackCand();
  //  virtual InoTrackCand *DupHandle() const;

  virtual void Trace(const char *c = "") const;

  //  static NavKey KeyFromSlice(const InoTrackCand *);

  void SetU(Int_t,Float_t);
  void SetV(Int_t,Float_t);
  void SetdS(Int_t,Float_t);
  void SetTrackPointXError(Int_t,Float_t);
  void SetTrackPointYError(Int_t,Float_t);
  void SetRange(Int_t plane,Float_t g_cm2);

  void Set2dS(Int_t,Float_t);
  void Set2Range(Int_t plane,Float_t g_cm2);

  int GetFCPC() const;
  void SetFCPC(int);

  void SetT(Int_t,Double_t);
  //  Bool_t BelongsWithTrack(InoTrackCand * trk, 
  //                          Double_t tolTPos2, Double_t tolZPos, Double_t tolTime){return true;};
  //  Bool_t BelongsWithShower(InoShowerCand* shw, 
  //                           Double_t tolTPos2, Double_t tolZPos, Double_t tolTime){return false;};
  void SetVtxTrace(Double_t);
  void SetVtxTraceZ(Double_t);
  void SetVtxnActiveUpstream(Int_t);
  void SetEndTrace(Double_t);
  void SetEndTraceZ(Double_t);
  void SetEndnActiveDownstream(Int_t);
  void SetVtxDistToEdge(Double_t);
  void SetEndDistToEdge(Double_t);

  void SetStraightLineSlopeX(Double_t v1) {StraightLineSlopeX=v1;};
  void SetStraightLineSlopeY(Double_t v1) {StraightLineSlopeY=v1;};
  void SetStraightLineInterceptX(Double_t v1) {StraightLineInterceptX=v1;};
  void SetStraightLineInterceptY(Double_t v1) {StraightLineInterceptY=v1;};
  void SetStraightLineChi2X(Double_t v1) {StraightLineChi2X=v1;};
  void SetStraightLineChi2Y(Double_t v1) {StraightLineChi2Y=v1;};
  void SetStraightLineNhitsX(Int_t v1) {StraightLineNhitsX=v1;};
  void SetStraightLineNhitsY(Int_t v1) {StraightLineNhitsY=v1;};
  void SetStraightLineXExpec(Int_t v1) {StraightLineXExpec=v1;};
  void SetStraightLineYExpec(Int_t v1) {StraightLineYExpec=v1;};

  void SetSimpleRadii(Double_t v1) {SimpleRadii=v1;};
  void SetSimpleCurv(Double_t v1) {SimpleCurv=v1;};
  void SetSimpleX0(Double_t v1) {SimpleX0=v1;};
  void SetSimpleZ0(Double_t v1) {SimpleZ0=v1;};
  void SetSimpleChi2Pos(Double_t v1) {SimpleChi2Pos=v1;};
  void SetSimpleChi2Neg(Double_t v1) {SimpleChi2Neg=v1;};
  void SetSimpleChi2Cndn(Double_t v1) {SimpleChi2Cndn=v1;};
  void SetSimpleAvgXPos(Double_t v1) {SimpleAvgXPos=v1;};
  void SetSimpleAvgXNeg(Double_t v1) {SimpleAvgXNeg=v1;};
  void SetSimpleAvgXCndn(Double_t v1) {SimpleAvgXCndn=v1;};
  void SetSimpleAvgXMeas(Double_t v1) {SimpleAvgXMeas=v1;};
  void SetSimpleNhits(Int_t v1) {SimpleNhits=v1;};
  void SetSimpleXExpec(Int_t v1) {SimpleXExpec=v1;};
 
  Double_t GetStraightLineSlopeX() {return StraightLineSlopeX;};
  Double_t GetStraightLineSlopeY() {return StraightLineSlopeY;};
  Double_t GetStraightLineInterceptX() {return StraightLineInterceptX;};
  Double_t GetStraightLineInterceptY() {return StraightLineInterceptY;};
  Double_t GetStraightLineChi2X() {return StraightLineChi2X;};
  Double_t GetStraightLineChi2Y() {return StraightLineChi2Y;};
  Int_t GetStraightLineNhitsX() {return StraightLineNhitsX;};
  Int_t GetStraightLineNhitsY() {return StraightLineNhitsY;};
  Double_t GetStraightLineXExpec() {return StraightLineXExpec;};
  Double_t GetStraightLineYExpec() {return StraightLineYExpec;};

  Double_t GetSimpleRadii() {return SimpleRadii;};
  Double_t GetSimpleCurv() {return SimpleCurv;};
  Double_t GetSimpleX0() {return SimpleX0;};
  Double_t GetSimpleZ0() {return SimpleZ0;};
  Double_t GetSimpleChi2Pos() {return SimpleChi2Pos;};
  Double_t GetSimpleChi2Neg() {return SimpleChi2Neg;};
  Double_t GetSimpleChi2Cndn() {return SimpleChi2Cndn;};
  Double_t GetSimpleAvgXPos() {return SimpleAvgXPos;};
  Double_t GetSimpleAvgXNeg() {return SimpleAvgXNeg;};
  Double_t GetSimpleAvgXCndn() {return SimpleAvgXCndn;};
  Double_t GetSimpleAvgXMeas() {return SimpleAvgXMeas;};
  Int_t GetSimpleNhits() {return SimpleNhits;};
  Double_t GetSimpleXExpec() {return SimpleXExpec;};
  
  virtual void ClearMaps(); // clears all STL maps, including fdS and fRange
  virtual void ClearUVT(); // same as ClearMap, left for backward compatibility

  //  Bool_t IsXPosValid(Int_t) const; // returns true if U and V positions have been set for this plane
  //  Bool_t IsXPosValid(Int_t) const; // returns true if U and V positions have been set for this plane
  Float_t GetU(Int_t) const; // U position at specified plane
  Float_t GetV(Int_t) const; // V position at specified plane
  Float_t GetZ(Int_t) const; // Z position at specified plane

  //  Int_t GetNTrackPlane(PlaneView::PlaneView_t = PlaneView::kUnknown) const;

  Float_t GetTrackPointXError(Int_t) const; // X-position error at specified plane
  Float_t GetTrackPointYError(Int_t) const; // Y-position error at specified plane
  Double_t GetT(Int_t) const; // time at specified plane
  //  Double_t GetT(Int_t,StripEnd::StripEnd_t) const; // time at specified plane
  //  Double_t GetT(StripEnd::StripEnd_t,Int_t) const; // time at specified plane
  Float_t GetdS(Int_t = -1) const; // travel distance from vertex
  Float_t GetRange(Int_t = -1) const; // g/cm**2 from vertex

  Float_t Get2dS(Int_t = -1) const; // travel distance from vertex
  Float_t Get2Range(Int_t = -1) const; // g/cm**2 from vertex



  Double_t GetVtxTrace() const;
  Double_t GetVtxTraceZ() const;
  Int_t GetVtxnActiveUpstream() const;
  Double_t GetEndTrace() const;
  Double_t GetEndTraceZ() const;
  Int_t GetEndnActiveDownstream() const;
  Double_t GetVtxDistToEdge() const;
  Double_t GetEndDistToEdge() const;

  virtual Double_t GetScore() const; // larger is better

  //  Int_t IsInShower(CandStripHandle *) const;
  //  void SetInShower(CandStripHandle *, Int_t);

  Double_t GetMomentum() const;
  void SetMomentum(Double_t);

  Int_t GetNFinderHits() const {return fNFinderHits;};
  void SetNFinderHits(Int_t vv) {fNFinderHits=vv;};

  Double_t GetCircleRadii() const {return fCircleRadii;};
  Double_t GetCircleMom() const {return fCircleMom;};
  Double_t GetCircleChisq() const {return fCircleChisq;};

  void SetCircleRadii(Double_t vv) {fCircleRadii = vv;};
  void SetCircleMom(Double_t vv) {fCircleMom = vv;};
  void SetCircleChisq(Double_t vv) {fCircleChisq = vv;};

  Double_t GetNewMomentum() const;
  void SetNewMomentum(Double_t);
 
  Int_t GetNewFitOut() const;
  void SetNewFitOut(Int_t);

  Double_t GetTheta() const;
  void SetTheta(Double_t);

  Double_t GetPhi() const;
  void SetPhi(Double_t);

  Double_t GetdSExtra() const;
  void SetdSExtra(Double_t);

  Double_t GetRangeExtra() const;
  void SetRangeExtra(Double_t);

  Int_t GetFitType() const {return fType;}
  void SetFitType(Int_t typ) {fType = typ;}

  

  //  Bool_t IsUnphysical(Float_t trkFrac=0.660,Float_t asymCut=0.500,
  //                      Float_t xtalkFrac=0.660,Float_t xtalkCut=2.0){return false;};
  Bool_t IsContained();

  void SetNTrackStrip(Int_t);
  void SetNTrackDigit(Int_t);
  void SetNTimeFitDigit(Int_t);
  void SetTimeFitChi2(Double_t);
  void SetTimeForwardFitRMS(Double_t);
  void SetTimeForwardFitNDOF(Int_t);
  void SetTimeBackwardFitRMS(Double_t);
  void SetTimeBackwardFitNDOF(Int_t);

  Int_t GetNTrackStrip() const;
  Int_t GetNTrackDigit() const;
  Int_t GetNTimeFitDigit() const;
  Double_t GetTimeFitChi2() const;
  Double_t GetTimeForwardFitRMS() const;
  Int_t GetTimeForwardFitNDOF() const;
  Double_t GetTimeBackwardFitRMS() const;
  Int_t GetTimeBackwardFitNDOF() const;

  mutable map<const InoHit*, Int_t> fInShower;
  mutable map<Int_t,Float_t> fUPos;
  mutable map<Int_t,Float_t> fVPos;
  mutable map<Int_t,Float_t> fdS;                           // in meters
  mutable map<Int_t,Float_t> fRange;                       // in g/cm**2
  mutable map<Int_t,Float_t> fXPosError;
  mutable map<Int_t,Float_t> fYPosError;  
  mutable map<Int_t,Double_t> fTime[2];

  // in meters for being points, which is not used to 
  //update track parameter, but included in tracks
  mutable map<Int_t,Float_t> f2dS;
  mutable map<Int_t,Float_t> f2Range;    // in g/cm**2

  Double_t fVtxTrace;
  Double_t fVtxTraceZ;
  Double_t fEndTrace;
  Double_t fEndTraceZ;
  Double_t fVtxDistToEdge;
  Double_t fEndDistToEdge;
  Int_t fVtxnActiveUpstream;
  Int_t fEndnActiveDownstream;
  
  Int_t fNTrackStrip;               // # of strips that have InShower<=1
  Int_t fNTrackDigit;               // # of digits that have InShower<=1
  Int_t fNTimeFitDigit;          // # of digits used to determine timing
  Double_t fTimeFitChi2;
  Double_t fTimeForwardFitRMS;
  Int_t fTimeForwardFitNDOF;
  Double_t fTimeBackwardFitRMS;
  Int_t fTimeBackwardFitNDOF;

  Double_t fTimeSlope;
  Double_t fTimeOffset;
  Double_t fMomentum;
  Double_t fNewMomentum;
  Int_t fNewFitOut;
  Double_t fTheta;
  Double_t fPhi;

  Int_t fNFinderHits;

  Double_t fCircleRadii;	// in mm ?
  Double_t fCircleMom;		// in GeV with charge
  Double_t fCircleChisq;	// 

  Double_t NewTimeEndPlaneExp;
  Double_t NewTimeSlope;
  Double_t NewTimeIntercept;

  Double_t fdSExtra;
  Double_t fRangeExtra; 

  Double_t StraightLineSlopeX;
  Double_t StraightLineSlopeY;
  Double_t StraightLineInterceptX;
  Double_t StraightLineInterceptY;
  Double_t StraightLineChi2X;
  Double_t StraightLineChi2Y;
  Int_t StraightLineNhitsX;  
  Int_t StraightLineNhitsY;
  Double_t StraightLineXExpec;
  Double_t StraightLineYExpec;
    
  Double_t SimpleRadii;
  Double_t SimpleCurv;
  Double_t SimpleX0;
  Double_t SimpleZ0;
  Double_t SimpleChi2Pos;
  Double_t SimpleChi2Neg;
  Double_t SimpleChi2Cndn;
  Double_t SimpleAvgXPos;
  Double_t SimpleAvgXNeg;
  Double_t SimpleAvgXCndn;
  Double_t SimpleAvgXMeas;
  Int_t SimpleNhits;
  Double_t SimpleXExpec;
  
  InoTrack*     fTrack; 
  InoVertex*    fVertex;
  InoVertex*    fTerm;

  
  public : 
  Int_t GetNDaughters() const; 
  //////////////////////////////////////////////////
  //From  CandRecoHandle
  Int_t GetNStrip(Int_t i) const; 

  // intersection between two candreco objects
  //  Int_t GetNStrip(const CandRecoHandle *, 
  //                  PlaneView::PlaneView_t = PlaneView::kUnknown) const;

  Int_t GetNDigit(Int_t) const;
  
  Int_t GetNPlane(Int_t) const;
  
  Int_t GetBegPlane(Int_t) const;
  Int_t GetEndPlane(Int_t) const;   

  void SetVtxU(Double_t);
  Double_t GetVtxU() const;
  
  void SetVtxV(Double_t);
  Double_t GetVtxV() const;

 void SetVtxZ(Double_t);
  Double_t GetVtxZ() const;

  void SetVtxT(Double_t);
  Double_t GetVtxT() const;

  void SetVtxPlane(Int_t);
  Int_t GetVtxPlane() const;

  void SetEndU(Double_t);
  Double_t GetEndU() const;
  
  void SetEndV(Double_t);
  Double_t GetEndV() const;

  void SetEndZ(Double_t);
  Double_t GetEndZ() const;

  void SetEndT(Double_t);
  Double_t GetEndT() const;

  void SetEndPlane(Int_t);
  Int_t GetEndPlane() const;

  void SetVtxDirCosU(Double_t);
  Double_t GetVtxDirCosU() const;

  void SetVtxDirCosV(Double_t);
  Double_t GetVtxDirCosV() const;

  void SetVtxDirCosZ(Double_t);
  Double_t GetVtxDirCosZ() const;

  void SetEndDirCosU(Double_t);
  Double_t GetEndDirCosU() const;

  void SetEndDirCosV(Double_t);
  Double_t GetEndDirCosV() const;

  void SetEndDirCosZ(Double_t);
  Double_t GetEndDirCosZ() const;

  void SetDirCosU(Double_t);
  Double_t GetDirCosU() const;

  void SetDirCosV(Double_t);
  Double_t GetDirCosV() const;

  void SetDirCosZ(Double_t);
  Double_t GetDirCosZ() const;
 
  void SetTimeSlope(Double_t);
  Double_t GetTimeSlope() const;

  void SetTimeOffset(Double_t);
  Double_t GetTimeOffset() const;

  void SetNewTimeEndPlaneExp(Double_t ival1) {NewTimeEndPlaneExp=ival1;}
  Double_t GetNewTimeEndPlaneExp() {return NewTimeEndPlaneExp;};
  void SetNewTimeSlope(Double_t ival1) {NewTimeSlope=ival1;}
  Double_t GetNewTimeSlope() {return NewTimeSlope;};
  void SetNewTimeIntercept(Double_t ival1) {NewTimeSlope=ival1;}
  Double_t GetNewTimeIntercept() {return NewTimeIntercept;};


  //  CalTimeType::CalTimeType_t GetCalTimeType() const;

  //  void CalibrateSigMapped(UInt_t fEncoded,Float_t);
  //  void CalibrateMIP(UInt_t fEncoded,Float_t);

  Double_t GetPulse() const; // CalStripType::CalStripType_t = CalStripType::kMIP) const;
  Double_t GetPlanePulse(Int_t iplane) const; //, CalStripType::CalStripType_t = CalStripType::kMIP) const;
  //  Double_t GetStripCharge(const InoHit*) const;
  //  Double_t GetStripCharge(const InoHit*) const; //, CalStripType::CalStripType_t, StripEnd::StripEnd_t = StripEnd::kWhole) const;
  //  Double_t GetStripCharge(const CandStripHandle *, StripEnd::StripEnd_t, CalStripType::CalStripType_t = CalStripType::kMIP) const;

  //InoFitterTrackHandle

  void SetFinderMomentum(double);
  double GetFinderMomentum() const;

  void SetMomentumdS(double);
  double GetMomentumdS() const;

  void SetMomentumRange(double);
  double GetMomentumRange() const;
  
  void SetMomentumCurve(double);
  double GetMomentumCurve() const;

  void SetEndMomentumCurve(double);
  double GetEndMomentumCurve() const;

  double GetEMCharge() const;
  void SetEMCharge(double);

  double GetVtxQPError() const;
  void SetVtxQPError(double);

  double GetVtxUError() const;
  void SetVtxUError(double);
  
  double GetVtxVError() const;
  void SetVtxVError(double);
  
  double GetVtxdUError() const;
  void SetVtxdUError(double);
  
  double GetVtxdVError() const;
  void SetVtxdVError(double);  

  void SetBave(Double_t);
  Double_t GetBave() const;

  void SetEndQP(Double_t);
  void SetPlaneChi2(Int_t,Double_t);
  void SetPlaneQP(Int_t,Double_t);
  void SetPlaneStateVector(int,int,double);
  double GetPlaneStateVector(int,int) const;
  void SetPlaneCovMatrix(int,int,int,double);
  double GetPlaneCovMatrix(int,int,int) const;
  void SetEndUError(Double_t);
  void SetEndVError(Double_t);
  void SetEnddUError(Double_t);
  void SetEnddVError(Double_t);
  void SetEndQPError(Double_t);
  void SetNSwimFail(Int_t);

  Double_t GetEndQP() const;
  Float_t GetPlaneChi2(Int_t) const;
  Float_t GetPlaneQP(Int_t) const;
  Double_t GetEndUError() const;
  Double_t GetEndVError() const;
  Double_t GetEnddUError() const;
  Double_t GetEnddVError() const;
  Double_t GetEndQPError() const;
  Int_t GetNSwimFail() const;

  Int_t GetNhitsEndPlane() const {return NhitsEndPlane;};
  Int_t GetNhitsEndPlaneM1() const {return NhitsEndPlaneM1;};

  void SetNhitsEndPlane(Int_t v1) {NhitsEndPlane=v1;};
  void SetNhitsEndPlaneM1(Int_t v1) {NhitsEndPlaneM1=v1;};

  Double_t GetChi2() const;
  Double_t Getcval() const {return mclight;}
  void SetChi2(Double_t);
  void Setcval(Double_t c) {mclight= c;}

  Int_t GetNDOF() const;
  void SetNDOF(Int_t);

  Double_t GetRangeBiasedQP() const;
  void SetRangeBiasedQP(Double_t qp);

  Int_t GetNIterate() const;
  void SetNIterate(Int_t nit);

  vector<InoCluster*>ClustsInTrack;

  vector<float> xfitpos1;
  vector<float> yfitpos1;
  vector<float> zfitpos1;
  vector<int> zfitlay1;

  vector<float> filteredx0;
  vector<float> filteredx1;
  vector<float> filteredx2;
  vector<float> filteredx3;
  vector<float> filteredx4;
  vector<float> filteredmom;
  vector<float> filteredthe;
  vector<float> filteredphi;
  vector<int> filteredlay;

  vector<float> extrapolx0;
  vector<float> extrapolx1;
  vector<float> extrapolmom;
  vector<float> momvecdiff1;
  vector<float> radialdiff1;
  
  vector<InoHit*> HitsNotInTrack;

  unsigned int GetEntries() const {return ClustsInTrack.size();}
  private :
    
  double mFinderMomentum;  
  double mMomentumdS;
  double mMomentumRange;
  double mMomentumCurve;
  double mEndMomentumCurve;
  double mEMCharge;
  double mVtxQPError;
  
  int FCPC;
  Int_t mNDOF;
  double mVtxUError;
  double mVtxVError;
  double mVtxdUError;
  double mVtxdVError;

  int NhitsEndPlane;
  int NhitsEndPlaneM1;
  
  Double_t mEndQP;
  Double_t mBave;
  Float_t mPlaneChi2[500];
  Float_t mPlaneQP[500];
  double mPlaneStateVector[500][5];
  double mPlaneCovMatrix[500][5][5];
  Double_t mEndUError;
  Double_t mEndVError;
  Double_t mEnddUError;
  Double_t mEnddVError;
  Double_t mEndQPError;
  Int_t mNSwimFail;
  Double_t mChi2;
  Double_t mclight;
  Double_t fQP_rangebiased;
  Int_t fNIterate;
  Int_t fType;
  
};

#endif                                              // INOTRACKCAND_H
