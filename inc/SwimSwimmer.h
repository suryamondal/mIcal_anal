// Swimming through a particle, forward or backward
//
// Units are: position  : meters
//            momentum  : GeV/c
//            mass      : GeV/c^2
//            stepmax   : meters
//            stepmin   : meters
//
#ifndef SWIMSWIMMER_H
#define SWIMSWIMMER_H
#include "TVector3.h"
#include <cassert>
#include <cmath>
#include "vect_manager.h"
#include "TGeoManager.h"
#include "FieldPropagator.hh"
#include "DetectorParameterDef.hh" //AAR:these variables added to include variable airgap

//class BField;
//class VldContext;
class SwimParticle;

//......................................................................

class SwimSwimmer
{
public:
    // this constructor should be used for all reco purposes!!!!!!
    SwimSwimmer( double dist, double halfgap);
    SwimSwimmer( int plane, double dist, double halfgap);// const VldContext& vldc);

    // these two constructors are ONLY for special testing - they allow to force
    ~SwimSwimmer();
    void         SetNmaxStep(int n) { assert(n>0); fNmaxStep = n; }
    bool         SetStepper(const char* name = 0);
    void SetBPlane(int n) {BPlane=n;};
    bool         SwimForward(SwimParticle& particle, int& nextplane, double& b_ave);//, SwimCondition& c);
    bool         SwimBackward(SwimParticle& particle, int& nextplane, double& b_ave); //, SwimCondition& c);
    double         Swim(SwimParticle& particle, int& nextplane); //, SwimCondition& c);
    bool         SwimForward(SwimParticle& particle, double& b_ave);//, SwimCondition& c);
    bool         SwimBackward(SwimParticle& particle,double& b_ave); //, SwimCondition& c);
    double         Swim(SwimParticle& particle); //, SwimCondition& c);
    double         SwimExtrapolate(SwimParticle& particle);

    bool         SwimForwardExtrapolate(SwimParticle& particle, double& b_ave);//, SwimCondition& c);
    bool         SwimBackwardExtrapolate(SwimParticle& particle,double& b_ave); //, SwimCondition& c);

    //  SwimStepper* GetStepper() { return fStepper; }

    //inline void SwimStepData::SetStepper(SwimStepper* stepper)
    //{ fStepper = stepper; }

    inline void SetIsForward(bool isForward)
    { fIsForward = isForward; }

    //inline void SwimStepData::SetSwimMaterial(SwimGeo::SwimMaterial_t material)
    //{ fSwimMaterial = material; }

    inline void SetStepSize(double stepSize)
    { fStepSize = stepSize; }

    inline void SetSPI(int n)
    { fSPI = n; }

    TVector3 getCrossingShift() {return lastCrossingShift;}

    //protected:
    //  SwimStepper* CreateDefaultStepper();
    //  void         DeleteBfield();
    //  void         DeleteSwimGeo();
    //  void         DeleteStepper();

private:
    FieldPropagator *pFieldMap;
    void anal_getnrot(double*, double*);
    void anal_getarot(double, double, double*);
    void trace_track_planef(double*, double*, double, double*, double&);//AAR: 1st argument added to pass the Field array.
    void anal_rotme(double*, double*, double*);
    void anal_rotmet(double*, double*, double*);
    void track_move_pt_align(double*,double*, double, double*, double*, double&); //AAR: 1st argument added to pass the Field array.

    //  SwimStepData*   fStepData;   // contains stepper stepsize, geo, dir (owner)
    //  BField*         fMagField;   // The magnetic field
    //  SwimGeo*        fSwimGeo;    // The Swimmer's geometry
    //  SwimStepper*    fStepper;    // The integration driver
    double          fStepMax;    // maximum step size
    double          fStepMin;    // minimum step size
    double          fAcc;        // Requested fractional accuracy
    int             fNmaxStep;   // Maximum number of steps to allow
    bool            fOwnBfield;  // Own magnetic field
    bool            fOwnSwimGeo; // Own Swim geometry
    int BPlane;
    bool                    fIsForward;    // Swim direction in time
    //  SwimGeo::SwimMaterial_t fSwimMaterial; // material in the step
    double                  fStepSize;     // Stepsize
    int                     fSPI;          // Next SwimPlaneInterface to be crossed
    double          fDistance;  // Distance to travel in Z-direction
    double          fHalfLayerThickness;
    double fHalfAirGap;                //AAR:these variables added to include variable airgap
    DetectorParameterDef* paradef;//AAR:these variables added to include variable airgap
    TVector3    lastCrossingShift;
    TVector3    startposXYZ;
    int        Plane;
    int nbfield;
    double b_ave;
    TGeoManager* icalGeometry;
};

#endif // SWIMSWIMMER_H
