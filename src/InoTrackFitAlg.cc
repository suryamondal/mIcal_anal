//ASM 11Nov2009 : How track size increases with time  : event 55 of /home/gobinda/ilc/simul/inoical0/test_output_digi_051005.root
//    : Why cluster of so large distance is inside track : event 3

// GMA  Stripwidht = 2.00e-2 Need to put from database
// GMA Need great care for hits in only one plane
// GMA Here postion is taken from cluster position, but need to
//     Modify it to only hits position and choose hits in a layer
//     Where chisquare is minimum.
//     Also identify the criteria for hits with single strip and both strips

//GMAA  Zincreaseswith time ? Use timing infor doe this criteria in trackfinder
//GMAA 1. PossibleJoins[View].push_back(Seg0); in InoTrackFinder::FormTracks()
//GMAA 2. Missing layer in the time of fit !!!!!!!!!!!!!!!!1
//GMA PEthresh=0.00001; should be same as it was in trackfinder

#include <cassert>
#include <cmath>
#include "TMath.h"
#include "TVector3.h"
#include "InoTrackFitAlg.h"
#include "TRandom.h"
#include "vect_manager.h"
#include "InoNewFitAlg.h"
#include "SwimSwimmer.h"
#include "SwimParticle.h"
#include <sys/time.h>
#include "G4SystemOfUnits.hh"
#include "Fcnsg.h"
#define MINLAYER 3

#include "mystuff.h"
#include "data.h"
#include "circle.h"
#include "../inc/Utilities.cpp"
#include "../inc/CircleFitByChernovLesort.cpp"
#include "TVector3.h"

InoTrackFitAlg::InoTrackFitAlg()
{  // cout<<"  InoTrackFitAlg::InoTrackFitAlg() " <<endl;
  pAnalysis = MultiSimAnalysisDigi::AnPointer;   
  icalGeometry= (InoGeometry_Manager::APointer)->icalGeometry;
  pFieldMap = FieldPropagator::FdPointer;

  debug_fit = false;
  debug_fcpc = false;
  debug_new = false;//true;
  debug_mical = false;
  timecheck = true;
  TrkFitterDebug=0;//11;
  inoTrackCand_pointer = new InoTrackCand_Manager();

  DetectorParameterDef* paradef = DetectorParameterDef::AnPointer;

  StripXWidth = (1/m)*paradef->GetXStrwd();
  StripYWidth = (1/m)*paradef->GetYStrwd();
  nLayer      = paradef->GetnLayer();
  LayerThickness = (1/m)*2*(paradef->GetParlay(2)+paradef->GetParirlay(2));
  ShiftInX = (1/m)*(paradef->GetRPCShift(0) + paradef->GetStackShift(0));
  ShiftInY = (1/m)*(paradef->GetRPCShift(1) + paradef->GetStackShift(1));
  ShiftInZ = (1/m)*paradef->GetRPCShift(2);
  // cout<<"paradef->GetStackShift(2) "<<paradef->GetStackShift(2)<<endl;
  for (int ijk=0; ijk<nLayer; ijk++) {
    // cout<<"ShiftInZ("<<ijk<<") "<<(paradef->GetShiftInZ(ijk) + paradef->GetStackShift(2))/m<<endl;
    ZPosLayer[ijk] = (paradef->GetRPCLayerPosZ(ijk) + paradef->GetShiftInZ(ijk) + paradef->GetStackShift(2))/m; //AAR need to correct this   8/12/2011 for layer thickness 5.6cm this is fine
    // cout<<"ijk "<< ijk<<" "<<ZPosLayer[ijk]<<endl;
  }
  
  UseGeoSwimmer = 0;
  if(pAnalysis->isXtermOut==1){
    cout <<"\n InoTrackFitAlg::InoTrackFitAlg() "<<StripXWidth<<" "<<StripYWidth<<" "<<LayerThickness<<endl;
  }
}


InoTrackFitAlg::~InoTrackFitAlg() {
  //GMA need to clear after evergy events;
  for (unsigned int i=0; i< inoTrackCand_pointer->InoTrackCand_list.size(); i++) {
    if (inoTrackCand_pointer->InoTrackCand_list[i]) { 
      delete inoTrackCand_pointer->InoTrackCand_list[i]; inoTrackCand_pointer->InoTrackCand_list[i]=0;}}
  inoTrackCand_pointer->InoTrackCand_list.clear();
  if (inoTrackCand_pointer) {delete inoTrackCand_pointer; inoTrackCand_pointer=0;}

  //  for (unsigned int i=0; i< inoFittedTrack_pointer->InoFittedTrack_list.size(); i++) {
  //    if (inoFittedTrack_pointer->InoFittedTrack_list[i]) 
  //{ delete inoFittedTrack_pointer->InoFittedTrack_list[i]; inoFittedTrack_pointer->InoFittedTrack_list[i]=0;}}
  
  //  inoFittedTrack_pointer->InoFittedTrack_list.clear();
  //  if (inoFittedTrack_pointer) { delete inoFittedTrack_pointer; inoFittedTrack_pointer=0;}
  
  //  delete fFinderTrack;
  
}


void InoTrackFitAlg::RunAlg( ) {

  // cout<<"BsmrGaussPar = "<<BsmrGaussPar<<endl;
  // for (unsigned ijk=0; ijk<nLayer; ijk++) {
  //   // double smr11 = gRandom->Gaus(0.0,BsmrGaussPar);
  //   // if(ijk==0 || ijk==5 || ijk==10) {
  //   //   pFieldMap->BShift[ijk] = 0.0;
  //   // } else {
  //   //   pFieldMap->BShift[ijk] = 0.1;
  //   // }
  //   // cout<<"ijk "<<ijk<<" "<<pFieldMap->BShift[ijk]<<endl;
  //   // pAnalysis->pctbshift->Fill(pFieldMap->BShift[ijk]);
  // }

  // for (unsigned ijk=0; ijk<11; ijk++) {
  //   double tmBx,tmBy;
  //   double dxyz[3] = {1000.0,1000.0,20.0};
  //   pFieldMap->ElectroMagneticField(dxyz, tmBx,tmBy,ijk);
  //   tmBx *= 1000; tmBy *= 1000;
  //   pAnalysis->pctbshiftx->Fill(tmBx);
  //   pAnalysis->pctbshifty->Fill(tmBy);
  //   // cout<<ijk<<" "<<dxyz[0]<<" "<<dxyz[1]<<" "<<dxyz[2]<<" "<<tmBx<<" "<<tmBy<<endl;
  // }

  InoTrack_Manager *ptrackCollection = InoTrack_Manager::APointer; 
  if(debug_fit)
    cout <<"----------------------InoTrackFitAlg::RunAlg( ) "<<ptrackCollection->InoTrack_list.size()<<"----------------------"<<endl;
  if (ptrackCollection) {
    for (unsigned int itrk=0; itrk<ptrackCollection->InoTrack_list.size(); itrk++) {
      
      /*
	nbfield=0;
	bave=0;
	fFinderTrack = ptrackCollection->InoTrack_list[itrk];
	fTrackCand = new InoTrackCand(ptrackCollection->InoTrack_list[itrk]);
      */
      
      /*
	InoTrackCand* junkcand =new InoTrackCand();
	junkcand->SetVtxTrace(0.01);
	
	InoTrackCand junkcandx;
	junkcandx.SetVtxTrace(10.01);
	cout <<"junk "<<junkcandx.GetVtxTrace()<<" "<<junkcand->GetVtxTrace()<<endl;
      */
      
      // Initialisations
      
      //      for (int ij=0; ij<2; ij++) { 
      //GMA use time information to go for only one direction

      for (int ij=1; ij<2; ij++) {
	switch (ij)
	  {
	  case 0 : ZIncreasesWithTime=false; break; 
	  case 1 : ZIncreasesWithTime=true; break;
          default: ZIncreasesWithTime=false;
	  }
	//    if( fTrackCand->GetEndZ() > fTrackCand->GetVtxZ() )  {ZIncreasesWithTime=true;}
	//    else {ZIncreasesWithTime=false;}
	
	
	nbfield=0;
	bave=0;
 
	//MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;            //asm	
	fFinderTrack = ptrackCollection->InoTrack_list[itrk];
        
	// pAnalysis->nhits_finder[itrk]=fFinderTrack->ClustsInTrack.size();  //asm 
	ZIncreasesWithTime= DirectionFromFinderHits(ptrackCollection->InoTrack_list[itrk]);

	ZIncreasesWithTime=false;

	fTrackCand = new InoTrackCand(ptrackCollection->InoTrack_list[itrk], ZIncreasesWithTime); //VALGRIND

        if(debug_fit) {
	  cout <<"----------------------trksize "<<ij << " "<<ptrackCollection->InoTrack_list[itrk]->GetEntries()<<" "<<fTrackCand->GetEntries()<<" "<<ZIncreasesWithTime<<endl;
	  
	  for(unsigned int kpc=0; kpc<fTrackCand->ClustsInTrack.size(); kpc++) {
	    cout<<"clust "<<kpc<<" = ("<<fTrackCand->ClustsInTrack[kpc]->GetXPos()<<","<<fTrackCand->ClustsInTrack[kpc]->GetYPos()<<","<<fTrackCand->ClustsInTrack[kpc]->GetZPlane()<<");"<<endl;
	  }
	}
	fTrackCand->SetFitType((ZIncreasesWithTime) ? 1 : 0);
	
	SaveData=false;
	SwimThroughShower=false;
	PassTrack=true;
	
	MaxPlane=-20;
	MinPlane=nLayer;
	
	DeltaZ=-99;
	DeltaPlane=-99;
	ShowerEntryPlane=-99;
	
       	NIter=0;TotalNSwimFail=0; NumFinderStrips=0;
	x_k[5]=0; x_k_minus[5]=0;
	for (unsigned int i=0; i<5; ++i) {
	  x_k[i]=0; x_k_minus[i]=0; H_k[0][i]=0; H_k[1][i]=0; K_k[i][0]=0;K_k[i][1]=0;
	  VtxCov[i]=-999; EndCov[i]=-999;
	  for (unsigned int j=0; j<5; ++j) {
	    C_k[i][j]=0; C_k_minus[i][j]=0; C_k_intermediate[i][j]=0;
	    F_k[i][j]=0; F_k_minus[i][j]=0;
	    Q_k[i][j]=0; Q_k_minus[i][j]=0;
	    Identity[i][j]=0;
	  }
	}
	
	Identity[0][0]=1; Identity[1][1]=1; Identity[2][2]=1; Identity[3][3]=1; Identity[4][4]=1; 
	
	// Set initial parameters
	x_k_minus[0]=fTrackCand->GetVtxU();
	x_k_minus[1]=fTrackCand->GetVtxV();
	if(fTrackCand->GetVtxDirCosZ()!=0.) {
	  x_k_minus[2]=fTrackCand->GetVtxDirCosU()/fTrackCand->GetVtxDirCosZ();
	  x_k_minus[3]=fTrackCand->GetVtxDirCosV()/fTrackCand->GetVtxDirCosZ();
	}
	x_k_minus[4]=0.;
	x_k_minus[5]=0.;

	x_k4_biased=0;
	//	cout <<"Initial vector ";
	//	for (int i=0; i<6; i++) {
	//	  cout <<" "<<x_k_minus[i]<<" "<< x_k[i];
	//	}
	//	cout<<endl;
	
	// Run the high level methods
	InitialFramework_new();// slice,cx);

	//a	cout <<"Initial fameworkr "<<endl;
	//	GetFitData(MinPlane,MaxPlane);

	int status = RunCircleFit();

	RunTheFitter_new();	//VALFRIND
	  

        if (pAnalysis->ihist < pAnalysis->nhistmx-1 && ij==1 && pAnalysis->isVisOut>=2) {
          for (unsigned int i=0; i<nLayer; ++i) {
            for (unsigned int j=0; j<FilteredData[i].size(); j++) {
              //	      cout <<"FilteredData[i].size() "<<i<<" "<<j<<" "<<FilteredData[i].size()<<endl;
              if (InitTrkClustData[i].size()>0) {
                pAnalysis->gens_list[5][pAnalysis->ihist]->Fill(FilteredData[i][j].x_k0,
                                                                FilteredData[i][j].x_k1,
                                                                ZPosLayer[i]+0.05); //InitTrkClustData[i][0].csh->GetZPos()+0.05);

                vectGr  tmpgr;
                tmpgr.x = FilteredData[i][j].x_k0;
                tmpgr.y = FilteredData[i][j].x_k1;
                tmpgr.z = ZPosLayer[i]+0.05; //SlcClustData[i][0].csh->GetZPos()+0.05;
                tmpgr.dx = 0.0;
                tmpgr.dy = 0.0;
                tmpgr.dz = 0.0;
                //	      pAnalysis->fitr_vect.push_back(tmpgr);
                if (pAnalysis->isVisOut==3) pAnalysis->gens_vect[5].push_back(tmpgr);
              }
            }
          }
        }
	//a	cout <<"Initial fameworkxxr "<<pAnalysis->ihist<<endl;


        // Clear up
	//a	cout <<"after fitter "<<itrk<<" "<<ptrackCollection->InoTrack_list.size()<<" "<<fTrackCand->GetMomentumRange()<<" "<<fTrackCand->GetMomentum()<<" "<<fFinderTrack->GetBegXDir()<<endl;
        if (fTrackCand->GetNDOF()>0 && fTrackCand->GetNDOF()<1000 && fTrackCand->GetEntries()>2 && PassTrack) {
	  inoTrackCand_pointer->InoTrackCand_list.push_back(fTrackCand);
	} else {
	  fTrackCand=0;
	} 
	//a	cout <<" fittrer "<< inoTrackCand_pointer->InoTrackCand_list.size()<<endl;


	for (unsigned int i=0; i<doubleLa; ++i) {
	  InitTrkClustData[i].clear();
	  SlcClustData[i].clear();
	  TrkClustsData[i].clear();
	  FilteredData[i].clear();
	  ExtraPolData[i].clear();
	}
	
      }
    } 
    //-------------------------------------------------------------------------------------------------------------------------
    //asm: this part of the code was inserted to change the order of the tracks in trackCand vector.
    //Now the vectors are placed so that the segment close to the vertex part will taken placs as zeroth element of the vector.
    //I think this will not take care of hadron track.
    bool RepeatIteration=false;
    vector <InoTrackCand*>::iterator it;
    vector <InoTrack*>::iterator it_find;



    if(inoTrackCand_pointer->InoTrackCand_list.size()>1&&inoTrackCand_pointer->InoTrackCand_list.size()<3){
      InoTrackCand* tempTrk;
      InoTrack* tempFind;

      //cout<<"\n inoTrackCand_pointer->InoTrackCand_list.size(): "<< inoTrackCand_pointer->InoTrackCand_list.size()<<endl;

      for(unsigned int jj=0; jj<inoTrackCand_pointer->InoTrackCand_list.size()-1;jj++){
	//for(unsigned int jj=0; jj<1;jj++){

	for(unsigned int ii=jj+1; ii<inoTrackCand_pointer->InoTrackCand_list.size();ii++){


	  if(abs(inoTrackCand_pointer->InoTrackCand_list[jj]->GetVtxZ()-inoTrackCand_pointer->InoTrackCand_list[ii]->GetVtxZ())>0.5){
	    //cout<< abs(inoTrackCand_pointer->InoTrackCand_list[jj]->GetVtxZ()-inoTrackCand_pointer->InoTrackCand_list[ii]->GetVtxZ())<<endl;

	    if(inoTrackCand_pointer->InoTrackCand_list[jj]->GetFitType()==1 &&inoTrackCand_pointer->InoTrackCand_list[ii]->GetFitType()==1){
	      //cout<<"up"<<endl;
	      if(inoTrackCand_pointer->InoTrackCand_list[jj]->GetVtxZ()>inoTrackCand_pointer->InoTrackCand_list[ii]->GetVtxZ())
		{
		  //cout<< "1b4" <<inoTrackCand_pointer->InoTrackCand_list[jj]->GetVtxZ() <<" " <<inoTrackCand_pointer->InoTrackCand_list[ii]->GetVtxZ();
		  tempTrk = inoTrackCand_pointer->InoTrackCand_list[ii];
		  it=inoTrackCand_pointer->InoTrackCand_list.begin();
		  inoTrackCand_pointer->InoTrackCand_list.erase(it+ii);
		  inoTrackCand_pointer->InoTrackCand_list.insert(it+jj,tempTrk);

		  tempFind = ptrackCollection->InoTrack_list[ii];
		  it_find = ptrackCollection->InoTrack_list.begin();
		  ptrackCollection->InoTrack_list.erase(it_find+ii);
		  ptrackCollection->InoTrack_list.insert(it_find+jj,tempFind);

		  //cout<< inoTrackCand_pointer->InoTrackCand_list[jj]->GetVtxZ() <<" " <<inoTrackCand_pointer->InoTrackCand_list[ii]->GetVtxZ();

		}
	    }
	    else if ( inoTrackCand_pointer->InoTrackCand_list[jj]->GetFitType()==0 && inoTrackCand_pointer->InoTrackCand_list[ii]->GetFitType()==0){
	      //cout<<"down"<<endl;
	      if(inoTrackCand_pointer->InoTrackCand_list[jj]->GetVtxZ()<inoTrackCand_pointer->InoTrackCand_list[ii]->GetVtxZ())
		{
		  //cout<<"0b4"<< inoTrackCand_pointer->InoTrackCand_list[jj]->GetVtxZ() <<" " <<inoTrackCand_pointer->InoTrackCand_list[ii]->GetVtxZ();

		  tempTrk = inoTrackCand_pointer->InoTrackCand_list[ii];
		  it=inoTrackCand_pointer->InoTrackCand_list.begin();
		  inoTrackCand_pointer->InoTrackCand_list.erase(it+ii);
		  inoTrackCand_pointer->InoTrackCand_list.insert(it+jj,tempTrk);

		  tempFind =ptrackCollection->InoTrack_list[ii];
		  it_find =ptrackCollection->InoTrack_list.begin();
		  ptrackCollection->InoTrack_list.erase(it_find+ii);
		  ptrackCollection->InoTrack_list.insert(it_find+jj,tempFind);

		}

	    }

	  }
	}

      }

      vector <InoTrack*>::iterator it2;
      double c=-0.06;double m=0.9;double maxDelta=5; // co efficients obtained from pathlength analysis
      for(unsigned int jj=0; jj<inoTrackCand_pointer->InoTrackCand_list.size()-1;jj++){

	//cout<< inoTrackCand_pointer->InoTrackCand_list[jj]->GetVtxZ() << " :"<< inoTrackCand_pointer->InoTrackCand_list[jj]->GetEndZ() <<endl;
	//cout<< inoTrackCand_pointer->InoTrackCand_list[jj+1]->GetVtxZ() << " : "<< inoTrackCand_pointer->InoTrackCand_list[jj+1]->GetEndZ() <<endl;
	//cout<< endl;
	double dzz=inoTrackCand_pointer->InoTrackCand_list[jj+1]->GetVtxZ()-inoTrackCand_pointer->InoTrackCand_list[jj]->GetEndZ();
	double dxx=pow(pow(inoTrackCand_pointer->InoTrackCand_list[jj+1]->GetVtxU()-inoTrackCand_pointer->InoTrackCand_list[jj]->GetEndU(),2)+
		       pow(inoTrackCand_pointer->InoTrackCand_list[jj+1]->GetVtxV()-inoTrackCand_pointer->InoTrackCand_list[jj]->GetEndV(),2)+
		       pow(inoTrackCand_pointer->InoTrackCand_list[jj+1]->GetVtxZ()-inoTrackCand_pointer->InoTrackCand_list[jj]->GetEndZ(),2),0.5);
	double dpp=(inoTrackCand_pointer->InoTrackCand_list[jj+1]->GetMomentumCurve()-inoTrackCand_pointer->InoTrackCand_list[jj]->GetEndMomentumCurve());
	double del =   fabs(dpp)- fabs(m*dxx+c);
	pAnalysis->trk_gap->Fill(del);
	if(inoTrackCand_pointer->InoTrackCand_list[jj]->GetFitType()==1 && inoTrackCand_pointer->InoTrackCand_list[jj+1]->GetFitType()==1){

	  if(dzz>0&& fabs(del)<maxDelta){
	    ptrackCollection->InoTrack_list[jj]->AddTrack(ptrackCollection->InoTrack_list[jj+1]);
	    it2=ptrackCollection->InoTrack_list.begin();
	    ptrackCollection->InoTrack_list.erase(it2+jj+1);
	    RepeatIteration=true;
	  } 
	  else if(dzz<0){
	    cout<<"FINAL NO OF TRACKS"<<  inoTrackCand_pointer->InoTrackCand_list.size() << endl;
	  }
	}
	else if (inoTrackCand_pointer->InoTrackCand_list[jj]->GetFitType()==0 &&inoTrackCand_pointer->InoTrackCand_list[jj+1]->GetFitType()==0){
	  if(dzz<0&& fabs(del)<maxDelta){
	    ptrackCollection->InoTrack_list[jj+1]->AddTrack(ptrackCollection->InoTrack_list[jj]);
	    it2=ptrackCollection->InoTrack_list.begin();
	    ptrackCollection->InoTrack_list.erase(it2+jj);
	    RepeatIteration=true;
	  };
	}
      }
    }

    RepeatIteration = false;

    if(RepeatIteration){
      cout<<"Repeating iteration...."<<endl;
      inoTrackCand_pointer->InoTrackCand_list.clear();
      //cout<<"repeatIteration "<< ptrackCollection->InoTrack_list.size() <<endl;
      for (unsigned int itrk=0; itrk<ptrackCollection->InoTrack_list.size(); itrk++) {

	for (int ij=1; ij<2; ij++) {
	  switch (ij)
	    {
	    case 0 : ZIncreasesWithTime=false; break;
	    case 1 : ZIncreasesWithTime=true; break;
	    default: ZIncreasesWithTime=false;
	    }


	  nbfield=0;
	  bave=0;



	  fFinderTrack = ptrackCollection->InoTrack_list[itrk];
	  pAnalysis->nhits_finder[itrk]=fFinderTrack->ClustsInTrack.size();  //asm
	  //AAR: don't need to repeat this
	  ZIncreasesWithTime= DirectionFromFinderHits(ptrackCollection->InoTrack_list[itrk]);
	  cout<<"Zincetime "<<ZIncreasesWithTime<<endl;

	
	  fTrackCand = new InoTrackCand(ptrackCollection->InoTrack_list[itrk], ZIncreasesWithTime); //VALGRIND
	  //	  if(debug_fit)
	    cout <<"----------------------trksize "<<ij << " "<<ptrackCollection->InoTrack_list[itrk]->GetEntries()<<" "<<fTrackCand->GetEntries()<<"----------------------"<<endl;
	  
	  fTrackCand->SetFitType((ZIncreasesWithTime) ? 1 : 0);

	  SaveData=false;
	  SwimThroughShower=false;
	  PassTrack=true;

	  MaxPlane=-20;
	  MinPlane=nLayer;

	  DeltaZ=-99;
	  DeltaPlane=-99;
	  ShowerEntryPlane=-99;

	  NIter=0;TotalNSwimFail=0; NumFinderStrips=0;
	  x_k[5]=0; x_k_minus[5]=0;
	  for (unsigned int i=0; i<5; ++i) {
	    x_k[i]=0; x_k_minus[i]=0; H_k[0][i]=0; H_k[1][i]=0; K_k[i][0]=0;K_k[i][1]=0;
	    VtxCov[i]=-999; EndCov[i]=-999;
	    for (unsigned int j=0; j<5; ++j) {
	      C_k[i][j]=0; C_k_minus[i][j]=0; C_k_intermediate[i][j]=0;
	      F_k[i][j]=0; F_k_minus[i][j]=0;
	      Q_k[i][j]=0; Q_k_minus[i][j]=0;
	      Identity[i][j]=0;
	    }
	  }

	  Identity[0][0]=1; Identity[1][1]=1; Identity[2][2]=1; Identity[3][3]=1; Identity[4][4]=1;

	  // Set initial parameters
	  x_k_minus[0]=fTrackCand->GetVtxU();
	  x_k_minus[1]=fTrackCand->GetVtxV();
	  if(fTrackCand->GetVtxDirCosZ()!=0.) {
	    x_k_minus[2]=fTrackCand->GetVtxDirCosU()/fTrackCand->GetVtxDirCosZ();
	    x_k_minus[3]=fTrackCand->GetVtxDirCosV()/fTrackCand->GetVtxDirCosZ();
	  }
	  x_k_minus[4]=0.;
	  x_k_minus[5]=0.;

	  x_k4_biased=0;

	  // Run the high level methods
	  InitialFramework_new();// slice,cx);
	  int status = RunCircleFit();
	  RunTheFitter_new();	//VALFRIND

	  if (pAnalysis->ihist < pAnalysis->nhistmx-1 && ij==1 && pAnalysis->isVisOut>=2) {
	    for (unsigned int i=0; i<nLayer; ++i) {
	      for (unsigned int j=0; j<FilteredData[i].size(); j++) {
		//	      cout <<"FilteredData[i].size() "<<i<<" "<<j<<" "<<FilteredData[i].size()<<endl;
		if (InitTrkClustData[i].size()>0) {
		  pAnalysis->gens_list[5][pAnalysis->ihist]->Fill(FilteredData[i][j].x_k0,
								  FilteredData[i][j].x_k1,
								  ZPosLayer[i]+0.05); //InitTrkClustData[i][0].csh->GetZPos()+0.05);

		  vectGr  tmpgr;
		  tmpgr.x = FilteredData[i][j].x_k0;
		  tmpgr.y = FilteredData[i][j].x_k1;
		  tmpgr.z = ZPosLayer[i]+0.05; //SlcClustData[i][0].csh->GetZPos()+0.05;
		  tmpgr.dx = 0.0;
		  tmpgr.dy = 0.0;
		  tmpgr.dz = 0.0;
		  //	      pAnalysis->fitr_vect.push_back(tmpgr);
		  if (pAnalysis->isVisOut==3) pAnalysis->gens_vect[5].push_back(tmpgr);
		}
	      }
	    }
	  }


	  // Clear up
	  //a	cout <<"after fitter "<<itrk<<" "<<ptrackCollection->InoTrack_list.size()<<" "<<fTrackCand->GetMomentumRange()<<" "<<fTrackCand->GetMomentum()<<" "<<fFinderTrack->GetBegXDir()<<endl;
	  if (fTrackCand->GetNDOF()>0 && fTrackCand->GetNDOF()<1000) {
	    inoTrackCand_pointer->InoTrackCand_list.push_back(fTrackCand);
	  } else {
	    fTrackCand=0;
	  }


	  for (unsigned int i=0; i<doubleLa; ++i) {
	    InitTrkClustData[i].clear();
	    SlcClustData[i].clear();
	    TrkClustsData[i].clear();
	    FilteredData[i].clear();
	    ExtraPolData[i].clear();
	  }

	}
      }
      //cout<<"OUT"<<endl;
    }

  }
}


void InoTrackFitAlg::InitialFramework_new( ) {

  if(debug_fit) 
    cout<<"----------------------InoTrackFitAlg::InitialFramework_new()----------------------" <<endl;
  // Store InoHit and make the strips accessible by plane number
  double MisalignmentError=0.0; //2.5e-5;

  // double MisalignmentError=1e-8; //1e-6; //1e-4; //4e-4; //1e-8; //4e-6;  // GMA need number from INO : Squared error for misalignment of strips
  //  double strXwd = StripXWidth; // 0.0196;  //GMA use common variable, this is in metre (NOT IN CM)
  //  double XposErrorSq = pow(strXwd/pow(12.,0.5),2.);
  //  double strYwd = StripYWidth; //0.0196;  //GMA use common variable, this is in metre (NOT IN CM)
  //  double YposErrorSq = pow(strYwd/pow(12.,0.5),2.);

  int SlicePlane; 

  //a  cout <<"inside InitialFramework "<<endl;
  // Store all clusters

  InoCluster_Manager *pinoclust = InoCluster_Manager::APointer;
  //a  cout <<"1size "<< pinoclust->InoCluster_list.size()<<endl;
  for (unsigned i=0; i<pinoclust->InoCluster_list.size(); i++) {
    SlicePlane=pinoclust->InoCluster_list[i]->GetZPlane();
    ClustStruct temp;
    temp.csh=pinoclust->InoCluster_list[i];
    SlcClustData[SlicePlane].push_back(temp);
  }

  int TrackPlane; 
  //  Store all track clusters found,
  //a  cout<<"InoTrackFitAlg::InitialFramework total hitlistx "<<pinoclust->InoCluster_list.size()<<" Hists in track "<<fFinderTrack->ClustsInTrack.size()<<endl;
  for (unsigned i=0; i<fFinderTrack->ClustsInTrack.size(); i++) {
    
    SlicePlane=fFinderTrack->ClustsInTrack[i]->GetZPlane();
    ClustStruct temp;
    temp.csh=fFinderTrack->ClustsInTrack[i];
    //    InitTrkClustData[SlicePlane].push_back(temp);

    TrackPlane=temp.csh->GetZPlane();

    TrkDataStruct tempdata;
    tempdata.numInList = i;
    tempdata.cltime =temp.csh->GetTime(); 

    tempdata.ZPos=ZPosLayer[TrackPlane];
    tempdata.PlaneView = temp.csh->GetView();

    tempdata.XPos=temp.csh->GetXPos();
    tempdata.XPosErrSq = pow(temp.csh->GetXPosErr(),2.0) + MisalignmentError ; // + XposErrorSq; // pow(temp.csh->GetXPosErr(),2.0);

    tempdata.YPos=temp.csh->GetYPos();
    tempdata.YPosErrSq = pow(temp.csh->GetYPosErr(),2.0) + MisalignmentError; //  + YposErrorSq; // pow(temp.csh->GetYPosErr(),2.0);
    
    tempdata.Straight = fFinderTrack->ClustsInTrack[i]->GetStraight();
    int ishift = (fFinderTrack->ClustsInTrack[i]->GetStraight()) ? 0 : shiftLa;
    TrkClustsData[TrackPlane+ishift].push_back(tempdata);
    InitTrkClustData[SlicePlane+ishift].push_back(temp);
    
    if (ishift >0) { 
      temp.csh->SetStraight(false);
    } else {
      temp.csh->SetStraight(true);
    }


    // Identify ends of initial track
    if (TrackPlane>MaxPlane) {MaxPlane=TrackPlane;}
    if (TrackPlane<MinPlane) {MinPlane=TrackPlane;}
  }

  //a cout <<"Exiting InitialFramework "<<MinPlane<<" "<<MaxPlane<<endl;
} 


void InoTrackFitAlg::ShowerStrips() {
  //a  cout <<"InoTrackFitAlg : ShowerStrips, Look for large vertex shower" << endl;
     
  // Initialisations
  int Increment; int NumberOfHits;
  int Plane; int NewPlane;

  int VtxShwWindow=8; 
  int HitsForShw=4;
  double PEThreshold=0.00001; //GMA  .1;

  if(ZIncreasesWithTime==true) {Plane=MinPlane; Increment=1;}
  else {Plane=MaxPlane; Increment=-1;}
  NewPlane=Plane;

  // Identify any vertex showers
  while(abs(Plane-NewPlane)<=VtxShwWindow && NewPlane>=MinPlane && NewPlane<=MaxPlane) {

    if(SlcClustData[NewPlane].size()>0) {
      NumberOfHits=0;

      // Set the number of hits on a plane required for the plane to be identified as 'in the 
      // shower'. We account for the gradient of the track, with the factor of 0.25 representing
      // the approximate ratio of strip thickness to strip width.
      if(FilteredData[NewPlane].size()>0) {  //GMA Need optimisation
	// GMA 07/02/2009 Excluding layer, which does not have any hit points, 
	// but in track extrapolation, it is stored.
	if (FilteredData[NewPlane][0].x_k4 !=0.0) {  
	  if(SlcClustData[NewPlane][0].csh->GetView()==2) { //GMA what is this view
	    HitsForShw=max(min(7,int( 4+(0.25*fabs(FilteredData[NewPlane][0].x_k2)) )),
			   min(7,int( 4+(0.25*fabs(FilteredData[NewPlane][0].x_k3)) )));
	  } else if (SlcClustData[NewPlane][0].csh->GetView()==0) {
	    HitsForShw=min(7,int( 4+(0.25*fabs(FilteredData[NewPlane][0].x_k2)) ));
	  } else {
	    HitsForShw=min(7,int( 4+(0.25*fabs(FilteredData[NewPlane][0].x_k3)) ));
	  }
	}
      } else {
	HitsForShw=4;
      }


      // Count number of strips on plane with greater than 2PEs
      for(unsigned int j=0; j<SlcClustData[NewPlane].size(); ++j) {
        if(SlcClustData[NewPlane][j].csh->GetPulse()>PEThreshold) {NumberOfHits++;}
      }
      
      // If a vertex shower is found, note that we should use the Swimmer
      // to find the most likely track strips inside the shower
      if(NumberOfHits>=HitsForShw) {ShowerEntryPlane=NewPlane; SwimThroughShower=true; break;}
      
      NewPlane+=Increment;
    }
    else {NewPlane+=Increment;}
  }
  
  // Find the plane at which the 'clean' section of track enters the shower
  if(SwimThroughShower==true) {
    NewPlane=ShowerEntryPlane+Increment;
    int PlanesSinceLastHit=0;
    int PlaneWindow=4;

    while(PlanesSinceLastHit<PlaneWindow && NewPlane>=MinPlane && NewPlane<=MaxPlane) {
      if(SlcClustData[NewPlane].size()>0) {
        NumberOfHits=0;
	
        // Account for gradient of track, as before
        if(FilteredData[NewPlane].size()>0) { 
	  // GMA 07/02/2009 Excluding layer, which does not have any hit points, 
	  // but in track extrapolation, it is stored.
	  if (FilteredData[NewPlane][0].x_k4 !=0.0) {  
	    if(SlcClustData[NewPlane][0].csh->GetView()==2) {
	      HitsForShw=max(min(7,int(4+(0.25*fabs(FilteredData[NewPlane][0].x_k2)) )),
			     min(7,int(4+(0.25*fabs(FilteredData[NewPlane][0].x_k3)) )));
	    } else if (SlcClustData[NewPlane][0].csh->GetView()==0) {
	      HitsForShw= min(7,int(4+(0.25*fabs(FilteredData[NewPlane][0].x_k2)) ));
	    } else {
	      HitsForShw=min(7,int(4+(0.25*fabs(FilteredData[NewPlane][0].x_k3)) ));
	    }
	  }
        } else {
	  HitsForShw=4;
	}


        // Count number of strips on plane with greater than 2PEs
        for(unsigned int j=0; j<SlcClustData[NewPlane].size(); ++j) {
          if(SlcClustData[NewPlane][j].csh->GetPulse()>PEThreshold) {NumberOfHits++;}
        }
        if(NumberOfHits>=HitsForShw) {
          ShowerEntryPlane=NewPlane; NewPlane+=Increment; PlanesSinceLastHit=0;
        }
        else {PlanesSinceLastHit++; NewPlane+=Increment;}
        
      }
      else {PlanesSinceLastHit++; NewPlane+=Increment;}
    }
  }
}


int InoTrackFitAlg::RunCircleFit() {

  int entries = fTrackCand->GetEntries();
  vector<TVector3> xyzpos, xyzerr;
  xyzpos.reserve(entries);
  xyzerr.reserve(entries);
  cout<<" entries "<<entries<<endl;
  if(entries<=MINLAYER) {return 1;}
  
  for (int ijk=entries-1; ijk>=0; ijk--) {
    // cout<<" ijk "<<ijk<<endl;
    int i = fTrackCand->ClustsInTrack[ijk]->GetZPlane();
    if (i <=int(nLayer)) {
      xyzpos.push_back({fTrackCand->ClustsInTrack[ijk]->GetXPos(),
			fTrackCand->ClustsInTrack[ijk]->GetYPos(),
			ZPosLayer[i]});
      xyzerr.push_back({pow(fTrackCand->ClustsInTrack[ijk]->GetXPosErr(),2),
			pow(fTrackCand->ClustsInTrack[ijk]->GetYPosErr(),2),
			0.});
      cout << " " << xyzpos.back().X()
	   << " " << xyzpos.back().Y()
	   << " " << xyzpos.back().Z()
	   << " " << xyzerr.back().X()
	   << " " << xyzerr.back().Y()
	   << " " << xyzerr.back().Z()
	   << endl;
    }
  }

  double halfLayerThickness = (ZPosLayer[1] - ZPosLayer[0])*0.5;

  double pos1[3];
  pos1[0] = xyzpos.front().X()*1000;
  pos1[1] = xyzpos.front().Y()*1000;
  pos1[2] = (xyzpos.front().Z() + halfLayerThickness)*1000;
  double Bx,By;
  pFieldMap->ElectroMagneticField( pos1, Bx,By, 0);
  cout<<" Bx "<<Bx<<" By "<<By<<endl;
  
  // for(int jk=0;jk<ndfi;jk++) {
  //   xyzpos.push_back(inPoints.xyzpos[ndfi-1-jk]);
  //   TVector3 xxx(inPoints.xyerr[ndfi-1-jk].X(),
  // 		 inPoints.xyerr[ndfi-1-jk].Y(),0.);
  //   xyzerr.push_back(xxx);
  //   // cout << " " << jk
  //   //      << " " << xyzpos[jk].X()/strpwidth
  //   //      << " " << xyzpos[jk].Y()/strpwidth
  //   //      << " " << xyzpos[jk].Z()/strpwidth
  //   //      << " " << xyzerr[jk].X()
  //   //      << " " << xyzerr[jk].Y()
  //   //      << " " << xyzerr[jk].Z()
  //   //      << endl;
  //   if(jk>1 &&
  //      // calPointDist(xyzpos[0],xyzpos[jk])>circleLen) {
  //      calPointDist(xyzpos[0],xyzpos[jk])>=calPointDist(xyzpos[0],xyzpos[1])*circlePt*1.1) {
  //     break;}
  // }
  // const int ndfi3 = int(xyzpos.size());
  // // cout<<" ndfi3 "<<ndfi3<<endl;
	       
  // TVector3 MomIniDir(xyzpos[1].X()-xyzpos[0].X(),
  // 		     xyzpos[1].Y()-xyzpos[0].Y(),
  // 		     xyzpos[1].Z()-xyzpos[0].Z());
  // MomIniDir *= 1./MomIniDir.Mag();
  // // cout << " MomIniDir " << MomIniDir.X()
  // //      << " " << MomIniDir.Y()
  // //      << " " << MomIniDir.Z() << endl;
	       
  // TVector3 MagField = GetMagneticField(xyzpos[0].X(),
  // 				       xyzpos[0].Y(),
  // 				       xyzpos[0].Z()+(airGap+ironThickness)*0.5);
  // // cout << " MagField " << MagField.X()
  // //      << " " << MagField.Y()
  // //      << " " << MagField.Z() << endl;
	
  // TVector3 MomAlongMagField = MomIniDir;
  // MomAlongMagField.RotateUz(MagField.Unit());
  // double pZ = MomAlongMagField.Z();
  // double pT = MomAlongMagField.Perp();
  // // cout << " pZ " << pZ << " pT " << pT << endl;
	       
  // TVector3 forceDirection = MomIniDir.Cross(MagField); // F = Q*PxB
  // double rotAngle = -forceDirection.Phi();
  // forceDirection.RotateZ(rotAngle); // Rotating around Z, making Y component zero
  // // cout << " rotAngle " << rotAngle*180./TMath::Pi() << endl;
  // // cout << " forceDirection " << forceDirection.X()
  // //      << " " << forceDirection.Y()
  // //      << " " << forceDirection.Z() << endl;
	       
  // for(int jk=0;jk<ndfi3;jk++) {
  //   xyzpos[jk].RotateZ(rotAngle);
  //   xyzerr[jk].RotateZ(rotAngle);
  //   // cout << " " << jk
  //   //      << " " << xyzpos[jk].X()/strpwidth
  //   //      << " " << xyzpos[jk].Y()/strpwidth
  //   //      << " " << xyzpos[jk].Z()/strpwidth
  //   //      << " " << xyzerr[jk].X()
  //   //      << " " << xyzerr[jk].Y()
  //   //      << " " << xyzerr[jk].Z()
  //   //      << endl;
  // }


  // vector<double> DataX3, DataY3, ErrY3;
  // DataX3.reserve(entries);
  // DataY3.reserve(entries);
  // ErrY3.reserve(entries);
    
  // for(int jk=0;jk<ndfi3;jk++) {
  //   DataX3[jk] = xyzpos[jk].Z();
  //   DataY3[jk] = xyzpos[jk].X();
  //   ErrY3[jk]  = fabs(xyzerr[jk].X());
  //   // cout << " " << jk << " " << DataX[jk] << " " << DataY[jk] << endl;
  // }
  // reals LambdaIni3=0.00001;
  // Data data13(ndfi3,DataX3,DataY3,ErrY3);
  // // Circle circle3,circleIni3(0.5,DataY3[0],2.);
  // // int codet = CircleFitByChernovLesort(data13,circleIni3,LambdaIni3,circle3);
  // // if(codet!=0) {continue;}
  // // cout<<" "<<iev<<" cradii "<<circle3.r<<" "<<circle3.a/strpwidth<<" "<<circle3.b/strpwidth<<" s "<<circle3.s<<endl;
  // // double cradii3 = circle3.r;
  // // double ca3 = circle3.a;
  // // double cb3 = circle3.b;
	       
  // double ttang = atan(-(DataX3[0]-DataX3[ndfi3-1])/(DataY3[0]-DataY3[ndfi3-1]));
  // Circle circle31,circleIni31(0.5*(DataX3[0]+DataX3[ndfi3-1])+4.*cos(ttang),0.5*(DataY3[0]+DataY3[ndfi3-1])+4.*sin(ttang),4.);
  // int code31 = CircleFitByChernovLesort(data13,circleIni31,LambdaIni3,circle31);
  // Circle circle32,circleIni32(0.5*(DataX3[0]+DataX3[ndfi3-1])-4.*cos(ttang),0.5*(DataY3[0]+DataY3[ndfi3-1])-4.*sin(ttang),4.);
  // int code32 = CircleFitByChernovLesort(data13,circleIni32,LambdaIni3,circle32);
	       
  // // if(fabs(circle31.s-circle32.s)>0.0001) {
  // // 	 cout<<" "<<iev<<endl;
  // // 	 cout<<" circle31 "<<circle31.r<<" "<<circle31.a/strpwidth<<" "<<circle31.b/strpwidth<<" s "<<circle31.s<<endl;
  // // 	 cout<<" circle32 "<<circle32.r<<" "<<circle32.a/strpwidth<<" "<<circle32.b/strpwidth<<" s "<<circle32.s<<endl;
  // // }
	       
  // int codet = (circle31.s<circle32.s)?code31:code32;
  // if(codet!=0) {continue;}
	       
  // double cradii3 = (circle31.s<circle32.s)?circle31.r:circle32.r;
  // double ca3 = (circle31.s<circle32.s)?circle31.a:circle32.a;
  // double cb3 = (circle31.s<circle32.s)?circle31.b:circle32.b;
  // double cs3 = (circle31.s<circle32.s)?circle31.s:circle32.s;
  // // cout<<" "<<iev<<" cradii "<<cradii3<<" "<<ca3/strpwidth<<" "<<cb3/strpwidth<<" s "<<cs3<<endl;
	       
  // double pXZ = 0.3*uniformField*cradii3;
  // double parMom = pXZ*sqrt(1. + (pZ/pT)*(pZ/pT));
  // // cout<<" "<<iev<<" pXZ "<<pXZ<<" parMom "<<parMom<<endl;
  // MomIniDir.SetMag(parMom);
	       
  // // TVector3 startPt(xyzpos.back().X()-cb3,
  // // 		 0.,xyzpos.back().Z()-ca3);
  // // TVector3   endPt(xyzpos.front().X()-cb3,
  // // 		 0.,xyzpos.front().Z()-ca3);
  // // double skangle = endPt.Angle(startPt);
  // // cout << " " << iev
  // //      << " skangle " << skangle*180./TMath::Pi()
  // //      << endl;
	       
  // TVector3   endPt(xyzpos.back().X()-cb3,
  // 		   xyzpos.back().Z()-ca3,0.);
  // TVector3 startPt(xyzpos.front().X()-cb3,
  // 		   xyzpos.front().Z()-ca3,0.);
  // double skangle = endPt.Phi()-startPt.Phi();
  // // cout << " 0 " << iev
  // //      << " skangle " << skangle*180./TMath::Pi()
  // //      << endl;
  // if(skangle<-TMath::Pi()) {
  //   skangle += 2.*TMath::Pi();
  // } else if(skangle>TMath::Pi()) {
  //   skangle -= 2.*TMath::Pi();}
  // // cout << " 1 " << iev
  // //      << " skangle " << skangle*180./TMath::Pi()
  // //      << endl;
  // int chargeVal = skangle>0.?1.:-1.;
  // // cout<<" chargeVal "<<chargeVal<<endl;
	       
  // csFill = cs3;
  // crFill = cradii3;
  // cmomFill = parMom*chargeVal;


  return 0;
  
}


void InoTrackFitAlg::RunTheFitter_new( ) {
  
  if(debug_fit) 
    cout<< "---------------------- InoTrackFitAlg::RunTheFitter_new( ) ----------------------" <<endl;
  //MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;
  //cout <<"InoTrackFitAlg : RunTheFitter_new, Call methods in the appropriate order" << endl;
  GetInitialCovarianceMatrix(true);
  //cout <<"InoTrackFitAlg : RunTheFitter_new GetInitialCovarianceMatrix"<<endl;
  const bool GoForward=true;
  const bool GoBackward=false;
  double StateVector[6]; double Prediction[6]={0,0,0,0,0,0};
  
  // Control the iterations backwards and forwards
  //  Detector::Detector_t detector = vldc->GetDetector();

  // Control iterations over a track for which ZIncreasesWithTime

  int niteration=5;
  double chisq_old=-100;
  int ndof_old = -100;
  int ndifchi = -1;

  if(ZIncreasesWithTime==true)
    {

      //First iteration
      NIter++;
      
      //Vtx to End

      SaveData=true;
      // cout<<"Storing Filtered Data ab1"<<endl;
      StoreFilteredData(MinPlane); 

      LastIteration=false;
      GoForwards_new(false); // true); 
      //GoForwards();


      ResetCovarianceMatrix();

      //a      
      cout <<"true============ResetCovarianceMatrix============================="<<endl;
      
      if(SwimThroughShower==true) {RemoveTrkHitsInShw();}
      
      for (unsigned int i=0; i<doubleLa; ++i) {
	FilteredData[i].clear();
	ExtraPolData[i].clear();
      }
      // cout<<"Storing Filtered Data ab2"<<endl;
      StoreFilteredData(MaxPlane); 
      GoBackwards_new(false);
      //GoBackwards();

      if(SwimThroughShower==true) {ShowerSwim();}

      ResetCovarianceMatrix();

      bool ClusterFound = true; //FindTheStrips(false); //cth,false);   //asm: here we can set some check to decide when no cluster found


      if(ClusterFound==true) { // Guard against finding no strips
        for(int nint=0;nint<=niteration;nint++) {  


	  // cout<<"-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-"<<endl;
	  // cout<<"nint = "<<nint<<" "<<niteration<<" ZIncreasesWithTime "<<ZIncreasesWithTime<<endl;
	  // cout<<"-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-"<<endl;
	  
	  double chisq_new = fTrackCand->GetChi2();
	  int ndof_new = fTrackCand->GetNDOF(); 	
	  
	  if (ndof_old == ndof_new && abs(chisq_old - chisq_new) < 0.01) { ndifchi=0;}  //Valgrind comments : Conditional jump or move depends on uninitialised value //VALGRIND 
	  
	  ndof_old = ndof_new;
	  chisq_old = chisq_new;
	  
	  //GMA stop loop if there is no change in hit points
          //Keep in mind SaveData=true (Data from only last iteration is stored);
	  NIter++;
          if((nint==niteration || ndifchi==0) && nint>1) LastIteration = true;        
	  if(debug_new) cout<<"nint "<<nint<<" "<<niteration<<" "<<NIter<<" "<<LastIteration<<endl;

	  if (nint>0) {
	    // cout<<"Get Fit Data new 1"<<endl;
	    GetFitData_new(MinPlane,MaxPlane);
	  }
	  
          if (MinPlane > MaxPlane) { cout<<" PassTrack 3" << MinPlane<< "<"<< MaxPlane<<" ievt "<<pAnalysis->ievt2<<endl;PassTrack=false; break;}
	  //a	  cout <<"2true============ResetCovarianceMatrix============================="<<endl;
          SaveData=false; 
          //if (nint==0) { GoForwards_new(true);} else { GoForwards_new(false);}
	  
	  fTrackCand->f2dS.clear();
	  fTrackCand->f2Range.clear();
	  
	  GoForwards_new(true);
          //GoForwards();
	  
	  if (LastIteration) {
	    if(debug_fcpc)  cout<<"XXXXXXXXXXXXXXXXXX"<<endl;
	    double fcpc_statevector[6];
	    if(FilteredData[MaxPlane].size()>0) {
	      fcpc_statevector[0] = FilteredData[MaxPlane][0].x_k0;
	      fcpc_statevector[1] = FilteredData[MaxPlane][0].x_k1;
	      fcpc_statevector[2] = FilteredData[MaxPlane][0].x_k2;
	      fcpc_statevector[3] = FilteredData[MaxPlane][0].x_k3;
	      fcpc_statevector[4] = FilteredData[MaxPlane][0].x_k4;
	      fcpc_statevector[5] = FilteredData[MaxPlane][0].x_k5;
	    }
	    FCorPCForward = CheckFCPCUpOrDn(fcpc_statevector, true, MaxPlane, true);
	    if(debug_fcpc)  cout<<"1FCorPCForward "<<FCorPCForward<<endl;
	    if(debug_fcpc)  cout<<"-------------------------------"<<endl;
	  }
	  
	  ResetCovarianceMatrix();
	  //Look on this
          
          //End back to vtx again
	  for (unsigned int i=0; i<doubleLa; ++i) { 
	    for (unsigned jk=0; jk<FilteredData[i].size(); jk++) {
	      if (FilteredData[i][jk].x_k5==1) {
		FilteredData[i].erase(FilteredData[i].begin()+jk);
		jk--;
	      }
	    }
	    
	    for (unsigned jk=0; jk<ExtraPolData[i].size(); jk++) {
	      if (ExtraPolData[i][jk].x_k5==1) {
		ExtraPolData[i].erase(ExtraPolData[i].begin()+jk);
		jk--;
	      }
	    }
	    
	    //	    for (unsigned jk=1; jk<FilteredData[i].size(); jk++) {
	    //	      if (FilteredData[i][jk].x_k5==0) {
	    //		FilteredData[i].erase(FilteredData[i].begin()+jk);
	    //		jk--;
	    //	      }
	    //	    }
	  } 
	  //	    if (TrkClustsData[i].size()>0) FilteredData[i].clear();}

	  //          if(nint==niteration || ndifchi==0 ) 
	  SaveData=true;
	  // cout<<"Storing Filtered Data ab3"<<endl;
	  StoreFilteredData(MaxPlane);    

	  fTrackCand->fdS.clear();
	  fTrackCand->fRange.clear();


	  GoBackwards_new(false);    
	  //	  GoBackwards();

	  if (LastIteration) {
	    if(debug_fcpc)  cout<<"XXXXXXXXXXXXXXXXXX"<<endl;
	    double fcpc_statevector[6];
	    if(FilteredData[MaxPlane].size()>0) {
	      fcpc_statevector[0] = FilteredData[MinPlane][0].x_k0;
	      fcpc_statevector[1] = FilteredData[MinPlane][0].x_k1;
	      fcpc_statevector[2] = FilteredData[MinPlane][0].x_k2;
	      fcpc_statevector[3] = FilteredData[MinPlane][0].x_k3;
	      fcpc_statevector[4] = FilteredData[MinPlane][0].x_k4;
	      fcpc_statevector[5] = FilteredData[MinPlane][0].x_k5;
	    }
	    FCorPCBackward = CheckFCPCUpOrDn(fcpc_statevector, false, MinPlane, false);
	    if(debug_fcpc)  cout<<"2FCorPCBackward "<<FCorPCBackward<<endl;
	    if(debug_fcpc)  cout<<"-------------------------------"<<endl;
	  }

	  // if (nint>0) {
	  //   // cout<<"Get Fit Data new 1"<<endl;
	  //   GetFitData_new(MinPlane,MaxPlane);
	  // }

	  ResetCovarianceMatrix();

          if(nint==0) x_k4_biased= x_k[4];

	  if ((nint == niteration || ndifchi==0) && nint>1) {

	    for(int i=0; i<6; ++i) {StateVector[i]=x_k_minus[i];}
	    
	    double loc_zend;
	    if(MinPlane>0) {
	      loc_zend = ZPosLayer[MinPlane] - ZPosLayer[MinPlane-1];
	    }
	    bool GetPrediction=Swim(StateVector, Prediction, MinPlane, loc_zend, GoBackward);
	    if(pAnalysis->isXtermOut){	  
	      if (GetPrediction) {
		for (int i=0; i<6; i++) { cout <<i<<" "<<StateVector[i]<<" "<< Prediction[i]<<endl;}
		cout <<" end "<< 1/StateVector[4]<<" "<<1/Prediction[4]<<endl;
	      }
	    }
	    break;
	  }
	  
	  //	  if (ndifchi==0) break;
        }
      }
      else {cout<<" PassTrack 3.1 " <<" ievt "<<pAnalysis->ievt2<<endl; PassTrack=false;}
    }
  
  
  // Control iterations over a track for which ZDecreasesWithTime
  else
    {
      //cout<< "iMinErr: ZIncreasesWithTime==false"<< endl;
      // First iteration
      NIter++;

      // Vtx to End

      SaveData=true;
      // cout<<"Storing Filtered Data ab4"<<endl;
      StoreFilteredData(MaxPlane); 

      LastIteration=false;

      GoBackwards_new(false); //VALGRIND
      //      GoBackwards(); 

      ResetCovarianceMatrix();
      // cout<<"ResetCovarianceMatrix "<<endl;
      // End back to vtx, swimming through any large vtx shower and
      // identifying any strips in ND spectrometer
      //GMAA for the time being do not use showerswim
      //      ShowerStrips();  // Try to identify vtx showers, now we have an idea of gradient
      if(SwimThroughShower==true) {RemoveTrkHitsInShw();}

      for (unsigned int i=0; i<doubleLa; ++i) {
	FilteredData[i].clear();
	//	for (unsigned jk=0; jk<FilteredData[i].size(); jk++) {
	//	  if (FilteredData[i][jk].x_k5==1) FilteredData[i].erase(FilteredData[i].begin()+jk);
	//	}	
      }
      //	if (TrkClustsData[i].size()>0) FilteredData[i].clear();}
      // cout<<"Storing Filtered Data ab5"<<endl;
      StoreFilteredData(MinPlane); 


      GoForwards_new(false);
      //      GoForwards();
      // cout<<"I am here..."<<endl;
      if(SwimThroughShower==true) {ShowerSwim();}

      ResetCovarianceMatrix();

      bool ClusterFound = true; //GMAA FindTheStrips(false); // cth,false);
      //      bool ClusterFound = FindTheStrips(false); // cth,false);
      // Second iteration

      if(ClusterFound==true) { // Guard against finding no strips
        for(int nint=0;nint<=niteration;nint++) {  

	  // cout<<"-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-"<<endl;
	  // cout<<"nint = "<<nint<<" "<<niteration<<" ZIncreasesWithTime "<<ZIncreasesWithTime<<endl;
	  // cout<<"-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-"<<endl;

	  double chisq_new = fTrackCand->GetChi2();
	  int ndof_new = fTrackCand->GetNDOF(); 	
	  
	  if (ndof_old == ndof_new && abs(chisq_old - chisq_new) < 0.01) { ndifchi=0;} //VALGRIND
	  
	  ndof_old = ndof_new;
	  chisq_old = chisq_new;

          if((nint==niteration || ndifchi==0) && nint>1) LastIteration = true;
	  NIter++;

	  //          for (unsigned int i=0; i<doubleLa; ++i) {if (TrkClustsData[i].size()>0) TrkClustsData[i].clear();}
	  //a	  cout <<"111true=====ResetCovarianceMatrix=== "<<NIter<<endl;

	  if (nint>0) {
	    // cout<<"Get Fit Data new 2"<<endl;
	    GetFitData_new(MinPlane,MaxPlane);
	  }

          if (MinPlane > MaxPlane) { cout<<" PassTrack 3.2 " << MinPlane<< "<"<< MaxPlane<<" ievt "<<pAnalysis->ievt2<<endl; PassTrack=false; break;}
          SaveData=false; 

	  fTrackCand->f2dS.clear();
	  fTrackCand->f2Range.clear();
	  //	  if (nint==0) { GoBackwards_new(true); } else { GoBackwards_new(false); }
	  GoBackwards_new(true);
	  //	  GoBackwards();
	  bool fcpc_col = false;
	  if (LastIteration && fcpc_col) {
	    if(debug_fcpc)  cout<<"XXXXXXXXXXXXXXXXXX"<<endl;
	    double fcpc_statevector[6];
	    cout<<"FilteredData[MaxPlane].size() "<<FilteredData[MaxPlane].size()<<" "<<MaxPlane<<endl;
	    if(FilteredData[MaxPlane].size()>0) {
	      cout<<"hhh..."<<endl;
	      fcpc_statevector[0] = FilteredData[MinPlane][0].x_k0;
	      cout<<"hhh...222"<<endl;
	      fcpc_statevector[1] = FilteredData[MinPlane][0].x_k1;
	      fcpc_statevector[2] = FilteredData[MinPlane][0].x_k2;
	      fcpc_statevector[3] = FilteredData[MinPlane][0].x_k3;
	      fcpc_statevector[4] = FilteredData[MinPlane][0].x_k4;
	      fcpc_statevector[5] = FilteredData[MinPlane][0].x_k5;
	    }
	    cout<<"Hhhh >>>>"<<endl;
	    FCorPCForward = CheckFCPCUpOrDn(fcpc_statevector, false, MinPlane, false);
	    if(debug_fcpc)  cout<<"3FCorPCForward "<<FCorPCForward<<endl;
	    if(debug_fcpc)  cout<<"-------------------------------"<<endl;
	  }

	  ResetCovarianceMatrix();

          // End to Vtx again
	  for (unsigned int i=0; i<doubleLa; ++i) {
	    for (unsigned jk=0; jk<FilteredData[i].size(); jk++) {
	      if (FilteredData[i][jk].x_k5==1) {
		FilteredData[i].erase(FilteredData[i].begin()+jk);
		jk--;
	      }
	    }
	    for (unsigned jk=0; jk<ExtraPolData[i].size(); jk++) {
	      if (ExtraPolData[i][jk].x_k5==1) {
		ExtraPolData[i].erase(ExtraPolData[i].begin()+jk);
		jk--;
	      }
	    }
	  }
	  //if (TrkClustsData[i].size()>0) FilteredData[i].clear();}
	  //          if(nint==niteration || ndifchi==0 )
	  SaveData=true;  
	  // cout<<"Storing Filtered Data ab6"<<endl;
          StoreFilteredData(MinPlane);

	  fTrackCand->fdS.clear();
	  fTrackCand->fRange.clear();

	  GoForwards_new(false); 
	  //	  GoForwards();
	  // cout<<"I am here...GoFOrwards..."<<NIter<<" "<<nint<<endl;
	  
	  if (LastIteration && fcpc_col) {
	    if(debug_fcpc)  cout<<"XXXXXXXXXXXXXXXXXX"<<endl;
	    double fcpc_statevector[6];
	    if(FilteredData[MaxPlane].size()>0) {
	      fcpc_statevector[0] = FilteredData[MaxPlane][0].x_k0;
	      fcpc_statevector[1] = FilteredData[MaxPlane][0].x_k1;
	      fcpc_statevector[2] = FilteredData[MaxPlane][0].x_k2;
	      fcpc_statevector[3] = FilteredData[MaxPlane][0].x_k3;
	      fcpc_statevector[4] = FilteredData[MaxPlane][0].x_k4;
	      fcpc_statevector[5] = FilteredData[MaxPlane][0].x_k5;
	    }
	    FCorPCBackward = CheckFCPCUpOrDn(fcpc_statevector, true, MaxPlane, true);
	    if(debug_fcpc)  cout<<"4FCorPCBackward "<<FCorPCBackward<<endl;
	    if(debug_fcpc)  cout<<"-------------------------------"<<endl;
	  }
	  
	  // if (nint>0) {
	  //   // cout<<"Get Fit Data new 2"<<endl;
	  //   GetFitData_new(MinPlane,MaxPlane);
	  // }

	  ResetCovarianceMatrix();
          if(nint==0) x_k4_biased= x_k[4];
	  
	  if ((nint == niteration || ndifchi==0 )&&nint>1) {
	    for(int i=0; i<6; ++i) {StateVector[i]=x_k_minus[i];}
	    
	    double loc_zend = 0.0;
	    if(MaxPlane<nLayer-1) {
	      loc_zend = ZPosLayer[MaxPlane+1] - ZPosLayer[MaxPlane];
	    }
	    bool GetPrediction=Swim(StateVector, Prediction, MaxPlane, loc_zend, GoForward);
	    
	    if(pAnalysis->isXtermOut){
              if (GetPrediction) {
		for (int i=0; i<6; i++) { cout <<i<<" "<<StateVector[i]<<" "<< Prediction[i]<<endl;}
		cout <<" endf "<< 1/StateVector[4]<<" "<<1/Prediction[4]<<endl;
	      }
            }
            
	    break;
	  }

	  //	  if (ndifchi==0) break;
        }
      }
      else {cout<<" PassTrack 3.3 " <<" ievt "<<pAnalysis->ievt2<< endl;PassTrack=false;}
    }
  bool fcpc_col = false;
  FCorPC = 0;
  if (LastIteration && fcpc_col) {
    FCorPC += FCorPCBackward;
    FCorPC <<=4;
    FCorPC += FCorPCForward;
    // FCorPC = (((FCorPCBackward<<4)&0x0000f0)|(FCorPCForward&0x0f));
    // cout<<"FCorPC = "<<FCorPC<<", checkfcorpc = "<<checkfcorpc<<endl;
    if(debug_fcpc) {
      cout<<" FCPC = "<<", FCorPC = "<<FCorPC<<", FCorPCForward = "<<FCorPCForward<<", FCorPCBackward = "<<FCorPCBackward<<endl;
      if(FCorPC>255 || FCorPC<0) {
	cout<<"Error Error Error Error Error.... 111 "<<FCorPC<<endl;
      }
    }
  }
  
  // Organise the output
  
  if(pAnalysis->isXtermOut){
    cout<<endl; 
    for (int ij=MinPlane; ij<=MaxPlane; ij++) { 
      cout<<"pln "<<ij<<" size "<<FilteredData[ij].size()<<endl;
      for (unsigned i=0; i<FilteredData[ij].size(); i++) {
	cout<<"iMax pl "<<ij<<" si "<<i<<" ("
	    <<FilteredData[ij][i].x_k5<<","      
	    <<FilteredData[ij][i].x_k0<<","
	    <<FilteredData[ij][i].x_k1<<","
	    <<FilteredData[ij][i].x_k2<<","
	    <<FilteredData[ij][i].x_k3<<","
	    <<1./(FilteredData[ij][i].x_k4)<<","
	    <<FilteredData[ij][i].x_k5<<endl;
      }
    }
  }  

  if(pAnalysis->isXtermOut==1){
    if (MaxPlane >=0 && MaxPlane <int(doubleLa)) {

      cout<<endl;
      cout<<"-------------------------------------------------------------------"<<endl; 
      cout<< " indx "<<" Pln no " <<"FilteredData.x_k5 "<<
	"x_k0 " << "x_k1 " << "x_k2 " << "x_k3 " << "1./x_k4 "<< "x_k5 "<<endl;
      cout<<endl;

      for (unsigned i=0; i<FilteredData[MaxPlane].size(); i++) {
	cout<<"iMax "<<i<<" "<<MaxPlane<<" "
	    <<FilteredData[MaxPlane][i].x_k5<<" "      
	    <<FilteredData[MaxPlane][i].x_k0<<" "
	    <<FilteredData[MaxPlane][i].x_k1<<" "
	    <<FilteredData[MaxPlane][i].x_k2<<" "
	    <<FilteredData[MaxPlane][i].x_k3<<" "
	    <<1./(FilteredData[MaxPlane][i].x_k4)<<" "
	    <<FilteredData[MaxPlane][i].x_k5<<endl;
      }
    }

    if (MinPlane >=0 && MinPlane <int(doubleLa)) {  
      for (unsigned i=0; i<FilteredData[MinPlane].size(); i++) {
	cout<<"iMin "<<i<<" "<<MinPlane<<" "
	    <<FilteredData[MinPlane][i].x_k5<<" "    
	    <<FilteredData[MinPlane][i].x_k0<<" "
	    <<FilteredData[MinPlane][i].x_k1<<" "
	    <<FilteredData[MinPlane][i].x_k2<<" "
	    <<FilteredData[MinPlane][i].x_k3<<" "
	    <<1./(FilteredData[MinPlane][i].x_k4)<<" "
	    <<FilteredData[MinPlane][i].x_k5<<endl;
      }
    }
  }
  // If the fit was successful
  if(x_k[4]!=0. && PassTrack==true) {
    
    //JAM modify tweak following range bias removal
    // Tweak q/p to remove offset
    //    x_k[4]*=1.01+(0.1*fabs(x_k[4]));
    
    x_k4_biased *=1.01+(0.1*fabs(x_k[4]));
    x_k[4] *=1.013;  //GMA how do we get 13% increase in Q/P
    
    // Find final strips and add them to the fitted track
    FillGapsInTrack();
    bool FinalClusterFound = true; // GMAA 
    
    //    bool FinalClusterFound = FindTheStrips(true); // GMA chech this cth,true);
    
    // If final strips found, set the fitted track properties
    if(FinalClusterFound==true) {

      //      for (unsigned int i=0; i<doubleLa; i++) {
      //	for (unsigned int j=0; j<InitTrkClustData[i].size(); j++) {
      //	  fTrackCand->ClustsInTrack.push_back(InitTrkClustData[i][j].csh);
      //	}
      //      }
      
      int NumInUView = fTrackCand->GetNPlane(0);
      int NumInVView = fTrackCand->GetNPlane(1);
	

      // cout <<"numview "<< NumInUView <<" "<<NumInVView<<" "<<fTrackCand->GetEntries()<<endl;

      if(NumInUView>1 && NumInVView>1) {
	//        cth.SetPass(1);
        //cout <<" "<<Prediction[0]<<" "<<Prediction[1]<<" "<<Prediction[2]<<" "<<Prediction[3]<<" "<<Prediction[4]<<endl;
	
	// cout<<"Setting track prop1"<<endl;
        SetTrackProperties(Prediction); //cth);

      } else { // else of if(NumInUView>1 && NumInVView>1) {
	cout<<" PassTrack 4" << " ievt "<<pAnalysis->ievt2<<endl; //asm: PassTrack 4 was noted frequently.
	PassTrack=false;
      } // if(NumInUView>1 && NumInVView>1) {
    } else { // else of if(FinalClusterFound==true) {
      // Otherwise fail the track at this final stage
      cout<<" PassTrack 6" <<" ievt "<<pAnalysis->ievt2<< endl;PassTrack=false;
    } // if(FinalClusterFound==true) {
  }
  
  
  // If the fit has failed (e.g. q/p is zero and/or u, v are nonsense)
  if(x_k[4]==0. || PassTrack==false) {

    //GMA need to reactivate again  
    // SetPropertiesFromFinderTrack(); //cth);
    
    /*
    // Remove any existing strips in the failed fitted track
    vector<CandStripHandle*> Daughters;

    TIter FitTrackStripItr = cth.GetDaughterIterator();
    while(CandStripHandle* FitTrackStrip = dynamic_cast<CandStripHandle*>(FitTrackStripItr()))
    {Daughters.push_back(FitTrackStrip);}

    for(unsigned int i=0; i<Daughters.size(); ++i) {cth.RemoveDaughter(Daughters[i]);}
    Daughters.clear();


    // Put strips from track finder in failed fitted track
    TIter TrkStripItr = track->GetDaughterIterator();
    while(CandStripHandle* TrkStrip = dynamic_cast<CandStripHandle*>(TrkStripItr()))
    {cth.AddDaughterLink(*TrkStrip);}
    
    // Set position/direction properties using the finder track
    cth.SetPass(0); 
    cth.SetMomentumCurve(0.); cth.SetEMPulse(0); cth.SetVtxQPError(-999.);
    SetPropertiesFromFinderTrack(cth);
    */
  }
  // cout<<"Complete... RunTheFitter_new..."<<endl;
}  


void InoTrackFitAlg::RemoveTrkHitsInShw() {
  // If the 'clean' section of track is large enough, remove the track finding
  // data for planes before the ShowerEntryPlane
  if(debug_fit)
    cout <<"----------------------InoTrackFitAlg : RemoveTrkHitsInShw, Discard track finding data in shower ----------------------" << endl;

  int NumTrackHitsLeft=0;
  
  if(ZIncreasesWithTime==true) {
    for(int i=ShowerEntryPlane; i<=MaxPlane; ++i) {
      if(TrkClustsData[i].size()>0) {NumTrackHitsLeft++;}
    }
  }
  else if(ZIncreasesWithTime==false) {
    for(int i=MinPlane; i<=ShowerEntryPlane; ++i) {
      if(TrkClustsData[i].size()>0) {NumTrackHitsLeft++;}
    }
  }
  
  // Carry out removal if there will be 6 or more strips left afterwards
  if(NumTrackHitsLeft>5) { 
    if(ZIncreasesWithTime==true) {
      for(int i=MinPlane; i<=ShowerEntryPlane; ++i) {TrkClustsData[i].clear();}
    }   
    else if(ZIncreasesWithTime==false) {
      for(int i=ShowerEntryPlane; i<=MaxPlane; ++i) {TrkClustsData[i].clear();    }
    }
  } 
  // Otherwise note that we should not run the ShowerSwim method
  else {
    cout <<"InoTrackFitAlg : RemoveTrkHitsInShw, not enough hits after removal. Must use all finder data." << endl;
    SwimThroughShower=false;
  }
  
  // Find the new max and min planes
  MaxPlane=-20; MinPlane=5000;
  for (int i=0; i<(int)nLayer; ++i) {   
    if(TrkClustsData[i].size()>0) {
      if(i>MaxPlane) {MaxPlane=i;}
      if(i<MinPlane) {MinPlane=i;}
    }
  }
}

void InoTrackFitAlg::ShowerSwim() {
  // Method is called if we have a large shower near the track vertex
  //
  // The Swimmer is used to find the most likely track strip in the shower
  // and this strip is added to the fit

  if(debug_fit){ 
    cout<<" =========================================================="<<endl;
    cout <<"InoTrackFitAlg : ShowerSwim, improved track finding in shower" << endl;
    cout<<" =========================================================="<<endl;
  }
  // Initialisations
  int Plane; int NewPlane;
  double StateVector[6]; double NState[6];
  bool GoForward; bool SwimBack;
  int PlanesSinceLastHit=0;
  int PlaneView;
  int Increment;
  
  double StripXDistance=0; double StripYDistance=0; 
  double MinXDistanceToStrip=999; double MinYDistanceToStrip=999.;
  //  double StripXWidth=2.00e-2; double StripYWidth=2.00e-2;// 4.108e-2;

  if(ZIncreasesWithTime==true) {
    GoForward=false; Plane=MinPlane; Increment=-1;
  } else {
    GoForward=true; Plane=MaxPlane; Increment=1;
  }

  NewPlane=Plane+Increment;

  // Continue until we reach a 4 plane window with no likely hit or we reach 
  // the end of the detector
  while(PlanesSinceLastHit<4 && NewPlane>0 && NewPlane<=(int)nLayer-5){ //145) { //GMA Put those number from database
    if(SlcClustData[NewPlane].size()>0) {
      
      //      if(SlcClustData[NewPlane][0].csh->GetPlaneView()==2) {PlaneView=0;}
      //      else {PlaneView=1;}
      PlaneView = SlcClustData[NewPlane][0].csh->GetView();
      
      
      for(int i=0; i<6; ++i) {StateVector[i]=x_k_minus[i];}

      SwimBack=Swim(StateVector, NState, Plane, NewPlane, GoForward);
      if(!SwimBack){break;}
      for(int i=0; i<6; ++i) {x_k[i]=NState[i];}
 
      // Find the closest strip (within a distance 'MinDistanceToStrip') and
      // temporarily store CandStripHandle
      // Results are very sensitive to value of MinDistanceToStrip
      InoCluster* CurrentClust=0;
      //GMA Original 0.0055, but do not have much clue about it
      // Is it (0.01*gap/stripwidth) ? Then for INO it is 0.0426
      //Was put 0.0852 also why donot remember now (30/01/08)

      MinXDistanceToStrip=(1.5*StripXWidth)+ fabs(0.0055*x_k[2]); //Original
      MinYDistanceToStrip=(1.5*StripYWidth)+ fabs(0.0055*x_k[3]);
      
      for(unsigned int j=0; j<SlcClustData[NewPlane].size(); ++j) {
        if (PlaneView%2==0) StripXDistance=fabs(SlcClustData[NewPlane][j].csh->GetXPos()-x_k[0]);
	if (PlaneView>0) StripYDistance=fabs(SlcClustData[NewPlane][j].csh->GetYPos()-x_k[1]);
        
        if(StripXDistance<MinXDistanceToStrip && StripYDistance<MinYDistanceToStrip) {
          if (PlaneView%2==0) MinXDistanceToStrip=StripXDistance;
	  if (PlaneView>0)    MinYDistanceToStrip=StripYDistance;
          CurrentClust=SlcClustData[NewPlane][j].csh;
        }
      }            
      
      // If we find a likely track strip, add it to the fit data and call the Kalman
      // update equations before repeating process to find next track strips in the shower

      //      cout <<"CurrentClust "<<int(CurrentClust)<<" "<<NewPlane<<" "<<SlcClustData[NewPlane].size()<<endl;
      //      cout <<"CurrentClust  "<<NewPlane<<" "<<SlcClustData[NewPlane].size()<<endl;
      if(CurrentClust) {
        ClustStruct temp;
        temp.csh = CurrentClust;
        InitTrkClustData[NewPlane].push_back(temp);
        
        // Convert the strip to data required for Kalman fit
	// cout<<"Get Fit Data new 3"<<endl;
        GetFitData_new(NewPlane,NewPlane);

        // Carry out the Kalman fit
	for (int i=0; i<2; i++) {
	  for (int j=0; j<5; j++) {
	    H_k[i][j]=0;
	  }
	}
        if (PlaneView%2==0) {H_k[0][0]=1;}
	if (PlaneView   >0) {H_k[1][1]=1;}
	//	  cout <<"InoTrackFitAlg.Showerswim : WARNING : PlaneView for hits is not matching with 0/1/2"<<endl;

	double tmp_ds, tmp_drange;
        bool CombiPropagatorOk=GetCombiPropagator(Plane,NewPlane,GoForward,&tmp_ds,&tmp_drange);
        
        if(CombiPropagatorOk ) {
          GetNoiseMatrix(Plane,NewPlane);
          ExtrapCovMatrix();
          CalcKalmanGain(NewPlane);
          UpdateStateVector(Plane,NewPlane,GoForward);
          UpdateCovMatrix();
          MoveArrays();
	  // cout<<"Storing Filtered Data ab7"<<endl;
          StoreFilteredData(NewPlane);
          
          if(ZIncreasesWithTime) {MinPlane=NewPlane; Plane=MinPlane;}
          else {MaxPlane=NewPlane; Plane=MaxPlane;}
          NewPlane=Plane+Increment;
          
          PlanesSinceLastHit=0;
        }
      } else {
	NewPlane+=Increment; PlanesSinceLastHit++;
      }
    } else {
      NewPlane+=Increment; PlanesSinceLastHit++;
    }
  }
  // Note that shower swim is complete
  SwimThroughShower=false; //GMA why false ?
}


void InoTrackFitAlg::GetFitData_new(int& Plane1, int& Plane2) {
  // Loop over the initial track strip data and create the final data for fitting
  // if(debug_fit)
  if(TrkFitterDebug>10) {
    cout <<"----------------------InoTrackFitAlg : GetFitData_new---------------------- 1,2 "<<Plane1<<" "<<Plane2<<endl ;
  }
  // Initialsmislaiisations
  double MisalignmentError=0.0; //2.5e-5; //  double MisalignmentError=1e-8; //1e-6; //1e-4; //4e-4; //1e-8; //4e-6;  // GMA need number from INO : Squared error for misalignment of strips


  //  double strXwd = StripXWidth; // 0.0196;  //GMA use common variable, this is in metre (NOT IN CM)
  //  double XposErrorSq = pow(strXwd/pow(12.,0.5),2.);
  //  double strYwd = StripYWidth; //0.0196;  //GMA use common variable, this is in metre (NOT IN CM)
  //  double YposErrorSq = pow(strYwd/pow(12.,0.5),2.);

  //  cout <<"errorsquare "<<XposErrorSq<<" "<<YposErrorSq<<endl;

  // Get the data for region between the planes specified

  //  cout<<"data "<< fTrackCand->ClustsInTrack.size()<<endl;
 
  Plane1=10000;
  Plane2 = -20;
  fTrackCand->ClustsInTrack.clear();
  for (unsigned int i=0; i<doubleLa; ++i) {
    InitTrkClustData[i].clear();
    TrkClustsData[i].clear();
  }
  
  for (unsigned ijk=0; ijk<doubleLa; ijk++) {
    int misLayer=1;
    if (FilteredData[ijk].size()==0) { misLayer=misLayer+1;continue; };
    // cout <<"ijlayer "<<ijk<<" "<<FilteredData[ijk].size()<<endl;
    
    for (unsigned kl=0; kl<FilteredData[ijk].size(); kl++) {
      // cout<<"ijk "<<ijk<<" "<<kl<<endl;
      double x1 = FilteredData[ijk][kl].x_k0;
      double y1 = FilteredData[ijk][kl].x_k1;
      int ijk_1=0;
  
      double x2,y2;
      double ctheta = -100;
      double dmn = 0.07; //Maximum 10 cm               asm: note this 270711
      // cout<<"ijk "<<ijk<<" "<<kl<<endl;
      if(ijk>0){
	if(FilteredData[ijk-misLayer].size()){
	  ijk_1=ijk-misLayer;
     
	  x2 = FilteredData[ijk_1][0].x_k0;
	  y2 = FilteredData[ijk_1][0].x_k1;
	  // cout<<x2<<" "<<y2<<" "<<ijk<<endl;
	  // cout<< pow(pow(x1-x2,2)+pow(y1-y2,2)+ pow(0.096,2),0.5)<<endl;
	  ctheta = 0.096/pow(pow(x1-x2,2)+pow(y1-y2,2)+ pow(0.096,2),0.5);
	  if(abs(ctheta)<0.3)dmn = misLayer*0.12;
	  else {dmn=misLayer*0.07;}
	  dmn = 0.07*misLayer/2;
	}
      }
      //int ij = ijk ;// (ijk >=shiftLa) ? ijk-shiftLa : ijk;
      int ij = (ijk >=shiftLa) ? ijk-shiftLa : ijk;
      int jkx = -1;
      for(unsigned int jk=0; jk<SlcClustData[ij].size(); ++jk) {

	// if(ij==3 || ij==4 || ij==6) {
	//   for(unsigned int sq1=0; sq1<SlcClustData[ij][jk].csh->HitsInCluster.size(); ++sq1) {
	//     if (SlcClustData[ij][jk].csh->HitsInCluster[sq1]->GetXPosErr() < 100) {
	//       cout<<"xxx... HitesInCluster["<<sq1<<"] = "<<SlcClustData[ij][jk].csh->HitsInCluster[sq1]->GetXStrip()->GetId()<<":"<<SlcClustData[ij][jk].csh->HitsInCluster[sq1]->GetYStrip()->GetId()<<":"<<SlcClustData[ij][jk].csh->HitsInCluster[sq1]->GetXStripNum()<<":"<<SlcClustData[ij][jk].csh->HitsInCluster[sq1]->GetYStripNum()<<endl;
	//     }
	//     if (SlcClustData[ij][jk].csh->HitsInCluster[sq1]->GetYPosErr() < 100) {
	//       cout<<"yyy... HitsInCluster["<<sq1<<"] = "<<SlcClustData[ij][jk].csh->HitsInCluster[sq1]->GetXStrip()->GetId()<<":"<<SlcClustData[ij][jk].csh->HitsInCluster[sq1]->GetYStrip()->GetId()<<":"<<SlcClustData[ij][jk].csh->HitsInCluster[sq1]->GetXStripNum()<<":"<<SlcClustData[ij][jk].csh->HitsInCluster[sq1]->GetYStripNum()<<endl;
	//     }
	//   }
	// }
	
	if(SlcClustData[ij][jk].csh->HitsInCluster.size()>4) continue;

	// if(SlcClustData[ij][jk].csh->GetXEntries()>3) continue;
	// if(SlcClustData[ij][jk].csh->GetYEntries()>3) continue;
	if(SlcClustData[ij][jk].csh->GetNXStripsInClust()>3) continue;
	if(SlcClustData[ij][jk].csh->GetNYStripsInClust()>3) continue;
	
	double dx = fabs(SlcClustData[ij][jk].csh->GetXPos()-x1); // /SlcClustData[ij][jk].csh->GetXPosErr();
	double dy = fabs(SlcClustData[ij][jk].csh->GetYPos()-y1); // /SlcClustData[ij][jk].csh->GetYPosErr();
	if(LastIteration==true) {
	  pAnalysis->hxpos_ext_kalman[SlcClustData[ij][jk].csh->GetZPlane()]->Fill((SlcClustData[ij][jk].csh->GetXPos()-x1)/StripXWidth);
	  pAnalysis->hypos_ext_kalman[SlcClustData[ij][jk].csh->GetZPlane()]->Fill((SlcClustData[ij][jk].csh->GetYPos()-y1)/StripYWidth);
	}
	// cout<<"ij,jk = "<<ij<<","<<jk<<endl;
	// cout <<"xmeaspos, xext, dx "<< SlcClustData[ij][jk].csh->GetXPos()<<" "<<x1<<" "<<dx<<endl;
	// cout<<"ymeaspos, yext, dy "<<SlcClustData[ij][jk].csh->GetYPos()<<" "<<y1<<" "<<dy<<endl;
        pAnalysis->trk_edge->Fill(misLayer, pow(dx*dx+dy*dy,0.5));
	// cout<<"dmn "<<jk<<" "<<dmn<<" "<<pow(dx*dx+dy*dy,0.5)<<endl;
        if(dmn > pow(dx*dx+dy*dy,0.5)){dmn = pow(dx*dx+dy*dy,0.5); jkx=jk;
	  // cout<<"jkx "<<jkx<<" "<<jk<<endl;
	}
	//else{ mn_gt= pow(dx*dx+dy*dy,0.5); }
      } 


      if (jkx>=0){
	// cout<<"jkx "<<jkx<<endl;
	fTrackCand->ClustsInTrack.push_back(SlcClustData[ij][jkx].csh);
	//	InitTrkClustData[ij].push_back(SlcClustData[ij][jkx]);
	//      const InoCluster* tempcls = SlcClustData[ij][jkx].csh;
	
	int TrackPlane= SlcClustData[ij][jkx].csh->GetZPlane();
	
	TrkDataStruct tempdata;
	tempdata.numInList = fTrackCand->GetEntries()-1;

	tempdata.cltime =SlcClustData[ij][jkx].csh->GetTime(); 
	tempdata.ZPos=ZPosLayer[TrackPlane];
	tempdata.PlaneView =  SlcClustData[ij][jkx].csh->GetView();
	
	tempdata.XPos=  SlcClustData[ij][jkx].csh->GetXPos();
	tempdata.XPosErrSq = pow(SlcClustData[ij][jkx].csh->GetXPosErr(),2.0) + MisalignmentError; //  + XposErrorSq; // pow(temp.csh->GetXPosErr(),2.0);
	
	tempdata.YPos= SlcClustData[ij][jkx].csh->GetYPos();
	tempdata.YPosErrSq = pow(SlcClustData[ij][jkx].csh->GetYPosErr(),2.0) + MisalignmentError ; // + YposErrorSq; // pow(temp.csh->GetYPosErr(),2.0);
	
	// cout <<"tmpdata "<< ij<<" "<< jkx<<" "<<tempdata.cltime <<" "<<SlcClustData[ij][jkx].csh->GetTime()<<" "<<tempdata.XPos<<" "<<tempdata.YPos<<" "<<tempdata.ZPos<< endl;

	
	tempdata.Straight = FilteredData[ijk][kl].x_k6;
	int ishift = (FilteredData[ijk][kl].x_k6) ? 0 : shiftLa;
	TrkClustsData[TrackPlane+ishift].push_back(tempdata);
	InitTrkClustData[ij+ishift].push_back(SlcClustData[ij][jkx]);
	if (ishift >0) { 
	  SlcClustData[ij][jkx].csh->SetStraight(false);
	} else {
	  SlcClustData[ij][jkx].csh->SetStraight(true);
	}

	//	TrkClustsData[TrackPlane].push_back(tempdata);
	
	//      cout <<"TrackPlane1 "<< TrackPlane<<" "<<Plane1<<" "<<Plane2<<endl; 

        //if (TrackPlane>Plane2 & ishift==0) Plane2 = TrackPlane;
        //if (TrackPlane<Plane1 & ishift==0) Plane1 = TrackPlane;
        if (TrackPlane>Plane2) Plane2 = TrackPlane;
        if (TrackPlane<Plane1) Plane1 = TrackPlane;


	// cout <<"TrackPlane2 "<< TrackPlane<<" "<<Plane1<<" "<<Plane2<<endl;

        // cout<< " dmn   :" << ijk <<" "<< dmn<<" "<<kl<< endl;
      } else {
        //pAnalysis->trk_egde->Fill(pow(pow(x1-x2,2)+pow(y1-y2,2)+ pow(0.096,2),0.5), );
        // cout<< "else dmn   :" << ijk <<" "<< dmn<<" "<<kl<<" "<<FilteredData[ijk].size()<<" "<<ExtraPolData[ijk].size()<< endl;
	
	//	FilteredData[ijk].clear(); 
        FilteredData[ijk].erase(FilteredData[ijk].begin()+kl);
	ExtraPolData[ijk].erase(ExtraPolData[ijk].begin()+kl);
	kl-- ;
      }
    }
      
  }

}



void InoTrackFitAlg::FillGapsInTrack() {
  // If there is no filtered data for a plane (between MinPlane and MaxPlane),
  // but this plane has hits in the slice, we interpolate from the nearest 
  // state vectors
  //
  // As with all filtered data, the interpolated data will be compared to
  // strip positions in the FindTheStrips method
  // if (debug_fit) 
  if(TrkFitterDebug>100) {

    cout <<"----------------------InoTrackFitAlg :: FillGapsInTrack----------------------" << endl;
  }
  
  int CurrentPlane; int ForwardsPlane; int BackwardsPlane;
  int Plane; int NewPlane;  bool GoForward;
  double StateVector[6]; double Prediction[6]; bool GetPrediction;
  
  for (int i=MinPlane; i<=MaxPlane; ++i) {   
    if(SlcClustData[i].size()>0) {
      if(FilteredData[i].size()==0) {
	
        // Find nearest filtered state vectors (within two planes) and ZPos differences
        // Forwards
        CurrentPlane=i+1; ForwardsPlane=-99;
	
        while(CurrentPlane<=MaxPlane && CurrentPlane<=(i+2)) {
          if(FilteredData[CurrentPlane].size()>0) {
	    ForwardsPlane=CurrentPlane; break;
	  } else {
	    CurrentPlane++;
	  }
        }
	
        // Backwards
        CurrentPlane=i-1; BackwardsPlane=-99;
	
        while(CurrentPlane>=MinPlane && CurrentPlane>=(i-2) ) {
          if(FilteredData[CurrentPlane].size()>0) {
	    BackwardsPlane=CurrentPlane; break;
	  } else {
	    CurrentPlane--;
	  }
        }
	
        // Find and store possible new filtered data, range and dS
        if(ForwardsPlane!=-99 && BackwardsPlane!=-99) {
	  
          // Swimmer method
          GetPrediction=false;
          NewPlane=i;
          if(ZIncreasesWithTime==true) {Plane=ForwardsPlane; GoForward=false;}
          else{Plane=BackwardsPlane; GoForward=true;}
          if(FilteredData[Plane].size()>0) {
            StateVector[0] = FilteredData[Plane][0].x_k0;
            StateVector[1] = FilteredData[Plane][0].x_k1;
            StateVector[2] = FilteredData[Plane][0].x_k2;
            StateVector[3] = FilteredData[Plane][0].x_k3;
            StateVector[4] = FilteredData[Plane][0].x_k4;           
	    StateVector[5] = (double)FilteredData[Plane][0].x_k5;    	    
            GetPrediction=Swim(StateVector, Prediction, Plane, NewPlane, GoForward);
            
            if(GetPrediction==true) {
              // Store possible new state vector
              FiltDataStruct temp;
              temp.x_k0 = Prediction[0];
              temp.x_k1 = Prediction[1];
              temp.x_k2 = Prediction[2]; 
              temp.x_k3 = Prediction[3];
              temp.x_k4 = Prediction[4];              
	      temp.x_k5 = int(Prediction[5]);   
	      temp.x_k6 = true; //11Nov2009
	      //	      FilteredData[i].clear(); // 110809
              FilteredData[i].push_back(temp);
              ExtraPolData[i].push_back(temp);
            }
          }
        }
      }
    }
  }
  if(TrkFitterDebug>100) {
    cout<<"----} FillGapsInTrack()"<<endl;
  }
}

        


bool InoTrackFitAlg::GetCombiPropagator(const int Plane, const int NewPlane, const bool GoForward, double* ddS, double* ddRange) {
  // Combination propagator, essentially same as SR propagator, but
  // generation of matrix reduces calls to swimmer by 80%

  if (debug_new) 
    cout <<" ----------------------InoTrackFitAlg : GetCombiPropagator ----------------------" << endl;

  for (int i=0; i<5; ++i) {
    for (int j=0; j<5; ++j) {
      F_k_minus[i][j]=0;  
    }
  }

  F_k_minus[0][0]=1; F_k_minus[1][1]=1; 
  F_k_minus[2][2]=1; F_k_minus[3][3]=1;
  
  DeltaZ=fabs(TrkClustsData[NewPlane][0].ZPos-TrkClustsData[Plane][0].ZPos);
  DeltaPlane=abs(NewPlane-Plane);

  // Swimmer section for last column
  double PState[6];  double NState[6];  double StateVector[6];
  double Increment=0.01;
  bool SwimInc=false; bool SwimDec=false;
  int nswimfail=0;

  
  if(GoForward==true) {F_k_minus[0][2]=DeltaZ; F_k_minus[1][3]=DeltaZ;}
  else if(GoForward==false) {F_k_minus[0][2]=-DeltaZ; F_k_minus[1][3]=-DeltaZ;}
  
  
  // Give swimmer fixed number of opportunities for successful swim
  while((SwimInc==false || SwimDec==false) && nswimfail<=10) {

    Increment=0.05*fabs(x_k_minus[4]); 
    if(Increment<.01) {Increment=.01;}
    
    //GMAA
    //    if (x_k_minus[4]==0) Increment=0.2;


    for(int j=0; j<6; ++j) {StateVector[j]=x_k_minus[j];}  
    
    // Increment then swim
    StateVector[4]+=Increment;
    SwimInc=Swim(StateVector, NState, Plane, NewPlane, GoForward, ddS, ddRange);
    
    StateVector[4]=x_k_minus[4];
    
    // Decrement then swim
    StateVector[4]-=Increment;
    SwimDec=Swim(StateVector, PState, Plane, NewPlane, GoForward, ddS, ddRange);
    
    // If swim failed, double momentum and swim again
    if(SwimInc==false || SwimDec==false) {

      // GMA reopen this      cout <<" InoTrackFitAlg : GetCombiPropagator, Swim failed - Double momentum and swim again" << endl;
      x_k_minus[4]*=0.5;
      nswimfail++; TotalNSwimFail++;
      break;  //AAR
    }
    
    // Form last row of propagator matrix.  Need to transpose to get proper Kalman F_k_minus
    else {

      if(Increment!=0.) {     //this should  be set to Increment >0.1
        for(int j=0; j<5; ++j) {
          F_k_minus[j][4] = (NState[j]-PState[j]) / (2*Increment);
	}

      }
      else {F_k_minus[4][4]=1;}
      //cout<<"CombipropagatorF error "<< Plane<<" "<<" "<<NewPlane<<endl;
      //cout<<" "<<PState[0]<<" "  <<PState[1]<<" "<<PState[2]<<" "<<PState[3]<<" "<<PState[4]<<endl;
      //cout<<" "<<NState[0]<<" "  <<NState[1]<<" "<<NState[2]<<" "<<NState[3]<<" "<<NState[4]<<endl;

    }

  } // End while statement

  if(nswimfail>=10) {cout <<" InoTrackFitAlg : GetCombiPropagator, nswimfail>10, fail track" << endl; return false;}




  // Display
  if(debug_fit) {
    cout<<"------------------------------------------------------------"<<endl;
    cout << "Combi F_k_minus" << endl;
    for(int i=0; i<5; ++i) {  
      for(int j=0; j<5; ++j) {
	cout << F_k_minus[i][j] << " ";
      }
      cout << endl;
    }
    cout<<"------------------------------------------------------------"<<endl;
  }
  if(debug_new) cout<<"combi propogator------------------------------------------------------------"<<endl;
  return true;
}

bool InoTrackFitAlg::Swim(double* StateVector, double* Output, const int Plane, 
			  const int NewPlane,const bool GoForward_loc, double* dS, double* Range, double* dE) {

  if(debug_new)   cout <<" ----------------------InoTrackFitAlg : Swim ----------------------"<< Plane<<" "<<NewPlane<<endl;
  double loc_layerthickness = 0.096; // in meter
  if(Plane>NewPlane) {
    loc_layerthickness = ZPosLayer[Plane] - ZPosLayer[NewPlane];
  } else if(Plane<NewPlane) {
    loc_layerthickness = ZPosLayer[NewPlane] - ZPosLayer[Plane];
  }
  double halfLayerThickness;
  if(ZIncreasesWithTime==true && GoForward_loc==false) {
    if(NewPlane>0) {
      halfLayerThickness= (ZPosLayer[NewPlane] - ZPosLayer[NewPlane-1])*0.5;
    }
  } else if(ZIncreasesWithTime==false && GoForward_loc==true) {
    if(NewPlane<nLayer-1) {
      halfLayerThickness = (ZPosLayer[NewPlane] - ZPosLayer[NewPlane+1])*0.5;
    }
  }


  // Initialisations
  // customize for bfield scaling.
  //  BField * bf = new BField(*vldc,-1,0);
  SwimSwimmer* myswimmer = new SwimSwimmer(loc_layerthickness, halfLayerThickness); //*vldc,bf);
  //if(debug_fit)
  int cplane = (ZIncreasesWithTime)? Plane+1 : Plane;
  myswimmer->SetBPlane(cplane);
  if(debug_new) {
    // cout <<" InoTrackFitAlg : Swim "<< Plane<<" "<<NewPlane<<endl;
    // cout <<" InoTrackFitAlg : Swim StateVector "<<StateVector[0]<<" "<<StateVector[1]<<" "<<StateVector[2]<<" "<<StateVector[3]<<" "<<StateVector[4]<<" "<<endl;
  }
  //  if(UseGeoSwimmer) GeoSwimmer::Instance()->Initialize(*vldc);
 
  //  double invSqrt2 = pow(1./2.,0.5);
  double charge = 0.;
  bool done = false;
    
  if(fabs(StateVector[4])>1.e-10) {

    // if(debug_new) cout<<"Swim( 1st option"<<endl;
    
    double modp = fabs(1./StateVector[4]);
    
    // Fix, to account for fact the cosmic muons could move in direction of negative z
    if(ZIncreasesWithTime==false) {modp=-modp;}
    
    double dsdz = pow((1.+pow(StateVector[2],2)+pow(StateVector[3],2)),0.5);

    double angle = 0;
    double ct=cos(angle);
    double st=sin(angle);

    //    double dxdz = invSqrt2*(StateVector[2]-StateVector[3]);
    //    double dydz = invSqrt2*(StateVector[2]+StateVector[3]);
    double dxdz = ct*StateVector[2]-st*StateVector[3];
    double dydz = st*StateVector[2]+ct*StateVector[3];

    // Set up current muon details
    if(StateVector[4]>0.) charge = 1.;
    else if(StateVector[4]<0.) charge = -1.;
    
    TVector3 position(ct*StateVector[0]-st*StateVector[1],
                      st*StateVector[0]+ct*StateVector[1],
                      ZPosLayer[Plane]); //SlcClustData[Plane][0].csh->GetZPos());

    TVector3 momentum(modp*(dxdz/dsdz),
                      modp*(dydz/dsdz),
                      modp/dsdz);

    //    TVector3 bfield = bf->GetBField(position);
    //TVector3 bfield(1.,1.,0.); //GMA-magnetic field  //AAR: commented out
    //    TVector3 bfield(1.5,0.,0.); //GMA-magnetic field
    //bave += TMath::Sqrt(bfield[0]*bfield[0]+bfield[1]*bfield[1]+bfield[2]*bfield[2]); //AAR: commented out
    //    bave += pow(bfield[0]*bfield[0]+bfield[1]*bfield[1]+bfield[2]*bfield[2],0.5);

    //nbfield++;   //AAR: commented out
    
    SwimParticle muon(position,momentum);
    muon.SetCharge(charge);
    //    cout <<"charge === "<<charge<<" "<<momentum.X()<<" "<<momentum.Y()<<" "<<momentum.Z()<<" "<<position.X()<<" "<<position.Y()<<" "<<position.Z()<<" "<<dxdz<<" "<<dydz<<" "<<dsdz<<" st "<<StateVector[0]<<" "<<StateVector[1]<<" "<<StateVector[2]<<" "<<StateVector[3]<<" "<<StateVector[4]<<" "<<muon.GetMomentum().Z()<<endl;
    // cout<<"SwimFunc: Plane "<<Plane<<" "<<StateVector[0]<<" "<<StateVector[1]<<" "<<muon.GetMomentumModulus()<<endl;
    //GMA    SwimZCondition zc(ZPosLayer[NewPlane]); //SlcClustData[NewPlane][0].csh->GetZPos());
    // Do the swim, accounting for direction of motion w.r.t time too
    if( (GoForward_loc==true && ZIncreasesWithTime==true)  || (GoForward_loc==false && ZIncreasesWithTime==false) ) {
      if(UseGeoSwimmer) {
	//        done = GeoSwimmer::Instance()->SwimForward(muon,ZPosLayer[NewPlane]); //SlcClustData[NewPlane][0].csh->GetZPos());
      } else {
        done = myswimmer->SwimForward(muon,t_bave);
      }
    } else if( (GoForward_loc==true && ZIncreasesWithTime==false)  || (GoForward_loc==false && ZIncreasesWithTime==true) ) {
      if(UseGeoSwimmer) {
	//        done = GeoSwimmer::Instance()->SwimBackward(muon,ZPosLayer[NewPlane]); //SlcClustData[NewPlane][0].csh->GetZPos());
      } else {
        done = myswimmer->SwimBackward(muon,t_bave);
      }
    }
 
     
    bave += t_bave;  //AAR:  added
    nbfield++;       //AAR:  added
    if(done==true) {

      double angle = 0;
      double ct=cos(angle);
      double st=sin(angle);
      if(muon.GetDirection().Z()!=0. && muon.GetMomentumModulus()!=0.) {
        Output[0]=(st*muon.GetPosition().Y()+ct*muon.GetPosition().X());
        Output[1]=(ct*muon.GetPosition().Y()-st*muon.GetPosition().X());
        Output[2]=(st*(muon.GetDirection().Y()/muon.GetDirection().Z())+ct*(muon.GetDirection().X()/muon.GetDirection().Z()));
        Output[3]=(ct*(muon.GetDirection().Y()/muon.GetDirection().Z())-st*(muon.GetDirection().X()/muon.GetDirection().Z()));
        Output[4]=muon.GetCharge()/muon.GetMomentumModulus();
	  if(debug_new) cout<<"SwimFunc: NewPlane "<<NewPlane<<" "<<Output[0]<<" "<<Output[1]<<" "<<muon.GetMomentumModulus()<<endl;
	// cout<<"Outmuon Info: ct : st : X : Y : Z "<<ct<<" : "<<st<<" : "<<muon.GetDirection().X()<<" : "<<muon.GetDirection().Y()<<" : "<<muon.GetDirection().Z()<<endl;
	// cout<<"Y/Z : X/Z : Output[2]=st*Y/Z + ct*X/Z "<<(muon.GetDirection().Y()/muon.GetDirection().Z())<<" : "<<(muon.GetDirection().X()/muon.GetDirection().Z())<<" : "<<Output[2]<<endl;
	// cout<<"Y/Z : X/Z : Output[3]=ct*Y/Z - st*X/Z "<<(muon.GetDirection().Y()/muon.GetDirection().Z())<<" : "<<(muon.GetDirection().X()/muon.GetDirection().Z())<<" : "<<Output[3]<<endl;
	// cout<<"atan2(Y/Z,X/Z) : theta : phi "<<180*TMath::ATan2(Output[3],Output[2])/3.14<<" : "<< muon.GetDirection().Theta()*180/3.14<<" : "<<muon.GetDirection().Phi()*180/3.14<<endl;

	Output[5]= StateVector[5];
        // Get range and dS from the Swimmer
        if(dS) {*dS=muon.GetS();} 
	if(Range) {*Range=muon.GetRange();} 
	if(dE){*dE=muon.GetMomentumModulus()-momentum.Mag();} 
	// cout<<""<<endl;
	//GMA put this more elegantly
	// fTrackCand->fdS[NewPlane] =muon.GetS();
	// fTrackCand->fRange[NewPlane] =muon.GetRange();
	
      } else {done=false;}
    }
  
  } else {

    // If infinite momentum, use straight line extrapolation
    
    // if(debug_new) cout<<"Swim( 2nd option"<<endl;

    double delz = LayerThickness;
    //    cout <<"delz "<< delz<<endl;
    if (SlcClustData[NewPlane].size()>0 && SlcClustData[Plane].size()>0) {
      //      delz = (SlcClustData[NewPlane][0].csh->GetZPos()-SlcClustData[Plane][0].csh->GetZPos());
      delz = ZPosLayer[NewPlane] - ZPosLayer[Plane];
    }
    
    //    cout <<"delz "<< delz<<endl;

    Output[0]=StateVector[0] + StateVector[2]*delz;
    Output[1]=StateVector[1] + StateVector[3]*delz;
    Output[2]=StateVector[2];
    Output[3]=StateVector[3];
    Output[4]=StateVector[4];
    Output[5]=StateVector[5];

    done=true;
  }

  if(debug_new) {
    cout <<" StateVector "<< StateVector[0]<<" "<<StateVector[1]<<" "<<StateVector[2]<<" "<<StateVector[3]<<" "<<StateVector[4]<<" "<<1/StateVector[4]<<endl;
    cout <<" Output "<< Output[0]<<" "<<Output[1]<<" "<<Output[2]<<" "<<Output[3]<<" "<<Output[4]<<" "<<1/Output[4]<<" "<<done<<endl;
  }
  //cout <<" Input S1 "<< StateVector[0]<<" "<<StateVector[1]<<" "<<StateVector[2]<<" "<<StateVector[3]<<" "<<StateVector[4]<<endl;
  //cout <<" Output S1 "<< Output[0]<<" "<<Output[1]<<" "<<Output[2]<<" "<<Output[3]<<" "<<Output[4]<<endl;
  delete myswimmer;
  //  delete bf;
  return done;
}


bool InoTrackFitAlg::Swim(double* StateVector, double* Output, const int Plane, 
			  const double zend, const bool GoForward_loc, double* dS, double* Range, double* dE) {
  double loc_layerthickness = 0.096; // in meter
  if(ZIncreasesWithTime==true && GoForward_loc==false) {
    if(Plane>0) {
      loc_layerthickness = ZPosLayer[Plane] - ZPosLayer[Plane-1];
    }
  } else if(ZIncreasesWithTime==false && GoForward_loc==true) {
    if(Plane<nLayer-1) {
      loc_layerthickness = ZPosLayer[Plane+1] - ZPosLayer[Plane];
    }
  }

  SwimSwimmer* myswimmer = new SwimSwimmer(zend, 0.5*loc_layerthickness); //*vldc,bf);
 
  double charge = 0.;
  bool done = false;
    
  if(fabs(StateVector[4])>1.e-10) {
    double modp = fabs(1./StateVector[4]);
    
    if(ZIncreasesWithTime==false) {modp=-modp;}
    
    double dsdz = pow((1.+pow(StateVector[2],2)+pow(StateVector[3],2)),0.5);

    double angle = 0;
    double ct=cos(angle);
    double st=sin(angle);

    double dxdz = ct*StateVector[2]-st*StateVector[3];
    double dydz = st*StateVector[2]+ct*StateVector[3];

    // Set up current muon details
    if(StateVector[4]>0.) charge = 1.;
    else if(StateVector[4]<0.) charge = -1.;
    
    TVector3 position(ct*StateVector[0]-st*StateVector[1],
                      st*StateVector[0]+ct*StateVector[1],
                      ZPosLayer[Plane]); //SlcClustData[Plane][0].csh->GetZPos());

    TVector3 momentum(modp*(dxdz/dsdz),
                      modp*(dydz/dsdz),
                      modp/dsdz);

    //TVector3 bfield(1.,1.,0.); //GMA-magnetic field  //AAR  commented out
    //    TVector3 bfield(1.5,0.,0.); //GMA-magnetic field
    //bave += //TMath::Sqrt(bfield[0]*bfield[0]+bfield[1]*bfield[1]+bfield[2]*bfield[2]); //AAR commented out
    //    bave += pow(bfield[0]*bfield[0]+bfield[1]*bfield[1]+bfield[2]*bfield[2],0.5);
    //nbfield++; //AAR commented out

    SwimParticle muon(position,momentum);
    muon.SetCharge(charge);

    //GMA    SwimZCondition zc(ZPosLayer[NewPlane]); //SlcClustData[NewPlane][0].csh->GetZPos());
    // Do the swim, accounting for direction of motion w.r.t time too
    if( (GoForward_loc==true && ZIncreasesWithTime==true)  || (GoForward_loc==false && ZIncreasesWithTime==false) ) {
      done = myswimmer->SwimForward(muon,t_bave);
    }
    else if( (GoForward_loc==true && ZIncreasesWithTime==false)  || (GoForward_loc==false && ZIncreasesWithTime==true) ) {
      done = myswimmer->SwimBackward(muon,t_bave);
    }
    bave += t_bave;  //AAR:  added
    nbfield++;       //AAR:  added
    if(done==true) {
      double angle = 0;
      double ct=cos(angle);
      double st=sin(angle);
      if(muon.GetDirection().Z()!=0. && muon.GetMomentumModulus()!=0.) {
        Output[0]=(st*muon.GetPosition().Y()+ct*muon.GetPosition().X());
        Output[1]=(ct*muon.GetPosition().Y()-st*muon.GetPosition().X());
        Output[2]=(st*(muon.GetDirection().Y()/muon.GetDirection().Z())+ct*(muon.GetDirection().X()/muon.GetDirection().Z()));
        Output[3]=(ct*(muon.GetDirection().Y()/muon.GetDirection().Z())-st*(muon.GetDirection().X()/muon.GetDirection().Z()));
        Output[4]=muon.GetCharge()/muon.GetMomentumModulus();
	Output[5]= StateVector[5];

        // Get range and dS from the Swimmer
        if(dS) {*dS=muon.GetS();} 
	if(Range) {*Range=muon.GetRange();} 
	if(dE){*dE=muon.GetMomentumModulus()-momentum.Mag();} 
	
	//GMA put this more elegantly
	fTrackCand->SetdSExtra(muon.GetS());
	fTrackCand->SetRangeExtra(muon.GetRange());
	
      }
      else {done=false;}
    }
  }

  else{
    // If infinite momentum, use straight line extrapolation
    Output[0]=StateVector[0] + StateVector[2]*zend;
    Output[1]=StateVector[1] + StateVector[3]*zend;
    Output[2]=StateVector[2];
    Output[3]=StateVector[3];
    Output[4]=StateVector[4];
    Output[5]=StateVector[5];

    done=true;
  }    
  //cout <<" Output "<< Output[0]<<" "<<Output[1]<<" "<<Output[2]<<" "<<Output[3]<<" "<<Output[4]<<endl;
  delete myswimmer;
  //  delete bf;
  return done;
}

bool InoTrackFitAlg::Swim(double* StateVector, double* Output, const double zzz, 
			  const int NewPlane,const bool GoForward_loc, double* dS, double* Range, double* dE) {
  if(debug_fit) cout <<" ----------------------InoTrackFitAlg : Swim, specified starting Z ----------------------" << endl;
  
  // Initialisations
  // customize for bfield scaling.
  //  BField * bf = new BField(*vldc,-1,0);
  
  // GMA Need to extrace proper Z values for a plane
  SwimSwimmer* myswimmer = new SwimSwimmer(fabs(LayerThickness*NewPlane-zzz), 0.5*LayerThickness);
  
  //  if(UseGeoSwimmer) GeoSwimmer::Instance()->Initialize(*vldc);
  
  //  double invSqrt2 = pow(1./2.,0.5);
  double charge = 0.;
  bool done = false;
  
  if(fabs(StateVector[4])>1.e-10) {
    double modp = fabs(1./StateVector[4]);
    
    // Fix, to account for fact the cosmic muons could move in direction of negative z
    if(ZIncreasesWithTime==false) {modp=-modp;}
    
    double dsdz = pow((1.+pow(StateVector[2],2)+pow(StateVector[3],2)),0.5);
    double angle=0;
    double ct=cos(angle);
    double st=sin(angle);
    double dxdz = ct*StateVector[2]-st*StateVector[3];
    double dydz = st*StateVector[2]+ct*StateVector[3];
    
    // Set up current muon details
    if(StateVector[4]>0.) charge = 1.;
    else if(StateVector[4]<0.) charge = -1.;
    
    TVector3 position(ct*StateVector[0]-st*StateVector[1],
                      st*StateVector[0]+ct*StateVector[1],
                      zzz);

    TVector3 momentum(modp*(dxdz/dsdz),
                      modp*(dydz/dsdz),
                      modp/dsdz);
    SwimParticle muon(position,momentum);
    muon.SetCharge(charge);
    //GMA    SwimZCondition zc(ZPosLayer[NewPlane]); //SlcClustData[NewPlane][0].csh->GetZPos());

   
    // Do the swim, accounting for direction of motion w.r.t time too
    if( (GoForward_loc==true && ZIncreasesWithTime==true)  || (GoForward_loc==false && ZIncreasesWithTime==false) ) {
      if(UseGeoSwimmer) {
	//        done = GeoSwimmer::Instance()->SwimForward(muon,ZPosLayer[NewPlane]); //SlcClustData[NewPlane][0].csh->GetZPos());
      } else {
        done = myswimmer->SwimForward(muon,t_bave);
      } 
    }
    else if( (GoForward_loc==true && ZIncreasesWithTime==false)  || (GoForward_loc==false && ZIncreasesWithTime==true) ) {
      if(UseGeoSwimmer) {
	//        done = GeoSwimmer::Instance()->SwimBackward(muon,ZPosLayer[NewPlane]); //SlcClustData[NewPlane][0].csh->GetZPos());

      } else {
        done = myswimmer->SwimBackward(muon,t_bave);
      }     
    }
    if(done==true) {
      double angle=0;
      double ct=cos(angle);
      double st=sin(angle);
      
      if(muon.GetDirection().Z()!=0. && muon.GetMomentumModulus()!=0.) {
        Output[0]=(st*muon.GetPosition().Y()+ct*muon.GetPosition().X());
        Output[1]=(ct*muon.GetPosition().Y()-st*muon.GetPosition().X());
        Output[2]=(st*(muon.GetDirection().Y()/muon.GetDirection().Z())+ct*(muon.GetDirection().X()/muon.GetDirection().Z()));
        Output[3]=(ct*(muon.GetDirection().Y()/muon.GetDirection().Z())-st*(muon.GetDirection().X()/muon.GetDirection().Z()));
        Output[4]=muon.GetCharge()/muon.GetMomentumModulus();
	Output[5]=StateVector[5];
        // Get range and dS from the Swimmer
        if(dS) {*dS=muon.GetS();} if(Range) {*Range=muon.GetRange();} if(dE){*dE=muon.GetMomentumModulus()-momentum.Mag();} 
      }
      else {done=false;}
    }
    
  }

  else{
    // If infinite momentum, use straight line extrapolation
    double delz = (ZPosLayer[NewPlane] -zzz); //SlcClustData[NewPlane][0].csh->GetZPos()-z);
    Output[0]=StateVector[0] + StateVector[2]*delz;
    Output[1]=StateVector[1] + StateVector[3]*delz;
    Output[2]=StateVector[2];
    Output[3]=StateVector[3];
    Output[4]=StateVector[4];
    Output[5]=StateVector[5];
    done=true;
  }    
  //cout <<" Output "<< Output[0]<<" "<<Output[1]<<" "<<Output[2]<<" "<<Output[3]<<" "<<Output[4]<<endl;
  delete myswimmer;
  //  delete bf;

  return done;
}

/*
//GMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
bool InoTrackFitAlg::Swim(double* StateVector, double* Output, const int Plane, 
                          const double Zend,const bool GoForward_loc, double* dS, double* Range, double* dE)
{
  cout <<" InoTrackFitAlg : Swim, specified end Z" << endl;

  // Initialisations
  // customize for bfield scaling.

  //  BField * bf = new BField(*vldc,-1,0);
  //GMA use proper calcualtion for Z-position using Plane number
  SwimSwimmer* myswimmer = new SwimSwimmer(fabs(LayerThickness*Plane-Zend), 0.5*LayerThickness);

  //  if(UseGeoSwimmer) GeoSwimmer::Instance()->Initialize(*vldc);

  //  double invSqrt2 = pow(1./2.,0.5);
  double charge = 0.;
  bool done = false;
    
  if(fabs(StateVector[4])>1.e-10) {
    double modp = fabs(1./StateVector[4]);
    
    // Fix, to account for fact the cosmic muons could move in direction of negative z
    if(ZIncreasesWithTime==false) {modp=-modp;}
  
    double angle=0;
    double ct=cos(angle);
    double st=sin(angle);
  
    double dsdz = pow((1.+pow(StateVector[2],2)+pow(StateVector[3],2)),0.5);
    double dxdz = ct*StateVector[2]-st*StateVector[3];
    double dydz = st*StateVector[2]+ct*StateVector[3];
  
    // Set up current muon details
    if(StateVector[4]>0.) charge = 1.;
    else if(StateVector[4]<0.) charge = -1.;
    
    TVector3 position(ct*StateVector[0]-st*StateVector[1],
                      st*StateVector[0]+ct*StateVector[1],
                      ZPosLayer[Plane]); //SlcClustData[Plane][0].csh->GetZPos());


    TVector3 momentum(modp*(dxdz/dsdz),
                      modp*(dydz/dsdz),
                      modp/dsdz);
    SwimParticle muon(position,momentum);
    muon.SetCharge(charge);
    //GMA    SwimZCondition zc(Zend);
 
    // Do the swim, accounting for direction of motion w.r.t time too
    if( (GoForward_loc==true && ZIncreasesWithTime==true)  || (GoForward_loc==false && ZIncreasesWithTime==false) ) {
      if(UseGeoSwimmer){
	//        done = GeoSwimmer::Instance()->SwimForward(muon,Zend);
      } else {
	done = myswimmer->SwimForward(muon);
      }   
    }
    else if( (GoForward_loc==true && ZIncreasesWithTime==false)  || (GoForward_loc==false && ZIncreasesWithTime==true) ) {
      if(UseGeoSwimmer){
	//        done = GeoSwimmer::Instance()->SwimBackward(muon,Zend);
      } else {
	done = myswimmer->SwimBackward(muon);
      }    
    }
    if(done==true) {
      double angle = 0;
      double ct=cos(angle);
      double st=sin(angle);

      if(muon.GetDirection().Z()!=0. && muon.GetMomentumModulus()!=0.) {
        Output[0]=(st*muon.GetPosition().Y()+ct*muon.GetPosition().X());
        Output[1]=(ct*muon.GetPosition().Y()-st*muon.GetPosition().X());
        Output[2]=(st*(muon.GetDirection().Y()/muon.GetDirection().Z())+ct*(muon.GetDirection().X()/muon.GetDirection().Z()));
        Output[3]=(ct*(muon.GetDirection().Y()/muon.GetDirection().Z())-st*(muon.GetDirection().X()/muon.GetDirection().Z()));
        Output[4]=muon.GetCharge()/muon.GetMomentumModulus();
	Output[5]=StateVector[5];
        // Get range and dS from the Swimmer
        if(dS) {*dS=muon.GetS();} if(Range) {*Range=muon.GetRange();} if(dE){*dE=muon.GetMomentumModulus()-momentum.Mag();}
      }
      else {done=false;}
    }
  } else {
    // If infinite momentum, use straight line extrapolation
    double delz = (Zend-ZPosLayer[Plane]); //SlcClustData[Plane][0].csh->GetZPos());
    Output[0]=StateVector[0] + StateVector[2]*delz;
    Output[1]=StateVector[1] + StateVector[3]*delz;
    Output[2]=StateVector[2];
    Output[3]=StateVector[3];
    Output[4]=StateVector[4];
    Output[5]=StateVector[5]; 
    done=true;
  }    

  delete myswimmer;
  //  delete bf;
  return done;
}
*/
//GMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

void InoTrackFitAlg::ResetCovarianceMatrix()
{
  // Simple method reset variables/arrays to allow propagation again
  
  DeltaPlane=0; DeltaZ=0; 
  GetInitialCovarianceMatrix(false);
}



void InoTrackFitAlg::GetInitialCovarianceMatrix(const bool FirstIteration)
{
 
  if (debug_fit)
    cout <<"---------------------- InoTrackFitAlg : GetInitialCovarianceMatrix---------------------- " << FirstIteration << endl;
  
  if(FirstIteration==true) {
    
    for(int i=0; i<5; ++i) {
      for(int j=0; j<5; ++j) {
        C_k_minus[i][j]=0.;
      }
    }
    
    //GMA use proper values for all these covariant matrix parameters
    //GMAA play with these parameters
    // Diagonal terms
    C_k_minus[0][0]=0.25; C_k_minus[1][1]=0.25; 
    C_k_minus[2][2]=100.; C_k_minus[3][3]=100.; 
    C_k_minus[4][4]=1.;
    
    // Off diagonal terms. Taken from SR - Origin of this?
    if(ZIncreasesWithTime==true) {
      C_k_minus[0][4]=7.5e-5;  C_k_minus[1][4]=7.5e-5;
      C_k_minus[4][0]=7.5e-5;  C_k_minus[4][1]=7.5e-5;
    }
    else if(ZIncreasesWithTime==false) {
      C_k_minus[0][4]=-7.5e-5;  C_k_minus[1][4]=-7.5e-5;
      C_k_minus[4][0]=-7.5e-5;  C_k_minus[4][1]=-7.5e-5;
    }
  } else if (FirstIteration==false) {
    // Results are very sensitive to this multiplication. A large number means
    // that further iterations start with the same uncertainties as the first,
    // albeit with improved "track finder" strips

    for(int i=0; i<5; ++i) {C_k_minus[i][i]*=100;}
    
    // Make sure not larger than very first covariance elements
    C_k_minus[0][0]=min(C_k_minus[0][0],0.25); C_k_minus[1][1]=min(C_k_minus[1][1],0.25);
    C_k_minus[2][2]=min(C_k_minus[2][2],100.); C_k_minus[3][3]=min(C_k_minus[3][3],100.);
    C_k_minus[4][4]=min(C_k_minus[4][4],1.);
    
    double cov_xqp = 7.5e-5;             // Taken from SR - Origin of this?
    
    for(int i=0; i<2; ++i){
      if(fabs(C_k_minus[i][4])>cov_xqp) C_k_minus[i][4] = (C_k_minus[i][4] > 0 ? cov_xqp : -cov_xqp);
      if(fabs(C_k_minus[4][i])>cov_xqp) C_k_minus[4][i] = (C_k_minus[4][i] > 0 ? cov_xqp : -cov_xqp);
    }

    cov_xqp /= 0.06;                     // Taken from SR - Origin of this?

    for(int i=2; i<4; ++i){
      if(fabs(C_k_minus[i][4])>cov_xqp) C_k_minus[i][4] = (C_k_minus[i][4] > 0 ? cov_xqp : -cov_xqp);
      if(fabs(C_k_minus[4][i])>cov_xqp) C_k_minus[4][i] = (C_k_minus[4][i] > 0 ? cov_xqp : -cov_xqp);
    }
  }
  
  // Display
  if(debug_fit) {
  cout<<"---------------------------------------------------------- "<<endl; 
   cout << "Initial covariance matrix" << endl;
    for(int p=0; p<5; ++p){
      for(int q=0; q<5; ++q){
        cout << C_k_minus[p][q] << " ";
      }  
      cout << endl;
    }
  cout<<"---------------------------------------------------------- "<<endl; 
  }

}

void InoTrackFitAlg::GetNoiseMatrix(const int Plane, const int NewPlane) {
  if(debug_new) cout<<"InoTrackFitAlg::GetNoiseMatrix(const int "<<Plane<<", const int "<<NewPlane<<")"<<endl;
  
  double LayerThickness_local;
  if(NewPlane>Plane) {
  LayerThickness_local = ZPosLayer[NewPlane] - ZPosLayer[Plane];
} else if(NewPlane<Plane) {
  LayerThickness_local = ZPosLayer[Plane] - ZPosLayer[NewPlane];
}

  double B[6];
  double x[3];
  double Bx,By;

  x[0]= (TrkClustsData[NewPlane][0].XPos+TrkClustsData[Plane][0].XPos)*0.5*1000;
  x[1]= (TrkClustsData[NewPlane][0].YPos+TrkClustsData[Plane][0].YPos)*0.5*1000;
  x[2]= (TrkClustsData[NewPlane][0].ZPos+TrkClustsData[Plane][0].ZPos)*0.5*1000;
  int cplane = (NewPlane>Plane)? Plane+1 : Plane-1;
  pFieldMap->ElectroMagneticField(x, Bx,By,cplane);

  // double tmppos_ll[3] = {x[0]*100, x[1]*100, x[2]*100};
  // double tmpdir_ll[3] = {0.0,0.0,1.0};
  // icalGeometry->InitTrack(tmppos_ll,tmpdir_ll );
  // // cout<<"Noise Matrix vol "<<icalGeometry->GetCurrentVolume()->GetName()<<endl;
  // cout<<"poscheck "<<ZPosLayer[NewPlane]<<" "<<ZPosLayer[Plane]<<" "<<0.5*(ZPosLayer[NewPlane]+ZPosLayer[Plane])<<endl;
  // cout<<"Noise Matrix vol "<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<icalGeometry->GetCurrentVolume()->GetName()<<" "<<Bx<<" "<<By<<endl;
  B[0] =Bx*1000;B[1]= By*1000; B[2]=0;
 // This method is essentially the same as in SR fitter
 //   cout <<" InoTrackFitAlg : GetNoiseMatrix" << endl;

  for (int p=0; p<5; ++p) {
    for (int q=0; q<5; ++q) {
      Q_k_minus[p][q]=0;    }
  }
  
  // Get gradients, etc from x_k_minus
  double dsdzSquared = 1.+pow(x_k_minus[2],2)+pow(x_k_minus[3],2);
  double dsdz = pow(dsdzSquared,0.5);

  // Implement noise matrix as in SR
  if (DeltaPlane!=-99 && DeltaZ!=-99) {  
    double qp = x_k_minus[4];
    if(fabs(qp)<0.01) qp = (qp>0 ? 0.01 : -0.01);
    //    int izdir = ((NewPlane-Plane)>0 ? 0 : 1);
  
    double LocalRadiationLength=(dsdz * double(DeltaPlane) * 3.41); 
    // GMA 3.41 radiation lengths per iron plane (Get it from database)
    
    double SigmaMS=(0.0136 * fabs(qp) * pow(LocalRadiationLength,0.5)
                    * (1. + (0.038 * log(LocalRadiationLength)) ));
    double SigmaMSSquared=pow(SigmaMS,2);

    double Sigma33Squared=0.5*SigmaMSSquared*dsdzSquared*(1.+pow(x_k_minus[2],2));

    double Sigma34Squared=0.5*SigmaMSSquared*dsdzSquared*(x_k_minus[2]*x_k_minus[3]);

    double Sigma44Squared=0.5*SigmaMSSquared*dsdzSquared*(1.+pow(x_k_minus[3],2));;

    //GMAAAAA
    // Treat steel as discrete scatter
    //    /*
    //    SwimGeo *spil = SwimGeo::Instance(*(const_cast<VldContext*>(vldc))); // Get edges of steel
    //    double z1 = spil->GetSteel(TrkClustsData[Plane][0].ZPos,izdir)->GetZ();
    //    double z2 = spil->GetSteel(z1,izdir)->GetZ();
    //    double dzscatter = fabs(TrkClustsData[NewPlane][0].ZPos-0.5*(z1+z2));
    //    */
    //GMA Distance between this scattering point and the next scintillator plane
    double dzscatter = 0.5*fabs(TrkClustsData[NewPlane][0].ZPos-TrkClustsData[Plane][0].ZPos);
    double dzscatter2 = pow(dzscatter,2);

    /*  
    UgliGeomHandle ugh(*vldc);
    BField bf(*vldc,-1,0);
    TVector3 uvz;

    uvz(0) = x_k_minus[0];
    uvz(1) = x_k_minus[1];
    uvz(2) = ugh.GetNearestSteelPlnHandle(TrkClustsData[Plane][0].ZPos).GetZ0();

    TVector3 buvz = bf.GetBField(uvz, true);
    buvz[0] *= 0.3;
    buvz[1] *= 0.3;
    buvz[2] *= 0.3;
    */

    double buvz[3] ={0.3*B[0], 0.3*B[1], 0.3*B[2]}; // GMA BX=By=1Tesla, Bz=0 and 0.3 factor comes from
                                     // the conversion factor, pt=0.3Br;


    // we assume the uncertainty in energy loss to follow a Landau distribution,
    // and use the FWHM from Ahlen (RMP, Vol. 52, p. 143, 1980)    

    double eloss = 0.090 * double(DeltaPlane); //GMA energy loss depends on layer thickness and track momentum
                                               // Only energy loss inside iron ?
    
                                               //AAR:won't this be function of angle 
  double sigmaeloss2 = 0.25*eloss*dsdz*qp*qp;
    sigmaeloss2 *= sigmaeloss2;


    // Fill elements of noise matrix
    Q_k_minus[0][0]=dzscatter2*Sigma33Squared;
    Q_k_minus[0][1]=dzscatter2*Sigma34Squared;
    Q_k_minus[0][2]=dzscatter*Sigma33Squared; //GMA is it -ve or +ve ? //GMAA
    Q_k_minus[0][3]=dzscatter*Sigma34Squared;
    Q_k_minus[0][4]=dzscatter*sigmaeloss2*double(DeltaPlane)*buvz[1]*LayerThickness_local*dsdz*(1.+pow(x_k_minus[2],2));

    // cout<<"Q_k_minus "<<dzscatter<<" "<<sigmaeloss2<<" "<<DeltaPlane<<" "<<double(DeltaPlane)<<" "<<buvz[1]<<" "<<LayerThickness_local<<" "<<dsdz<<" "<<(1.+pow(x_k_minus[2],2))<<endl;

    Q_k_minus[1][0]=dzscatter2*Sigma34Squared;
    Q_k_minus[1][1]=dzscatter2*Sigma44Squared;
    Q_k_minus[1][2]=dzscatter*Sigma34Squared;
    Q_k_minus[1][3]=dzscatter*Sigma44Squared;
    Q_k_minus[1][4]=dzscatter*sigmaeloss2*double(DeltaPlane)*buvz[0]*LayerThickness_local*dsdz*(1.+pow(x_k_minus[3],2));

    Q_k_minus[2][0]=dzscatter*Sigma33Squared;
    Q_k_minus[2][1]=dzscatter*Sigma34Squared;
    Q_k_minus[2][2]=Sigma33Squared;
    Q_k_minus[2][3]=Sigma34Squared;
    Q_k_minus[2][4]=sigmaeloss2*double(DeltaPlane)*buvz[1]*LayerThickness_local*dsdz*(1.+pow(x_k_minus[2],2));

    Q_k_minus[3][0]=dzscatter*Sigma34Squared;
    Q_k_minus[3][1]=dzscatter*Sigma44Squared;
    Q_k_minus[3][2]=Sigma34Squared;
    Q_k_minus[3][3]=Sigma44Squared;
    Q_k_minus[3][4]=sigmaeloss2*double(DeltaPlane)*buvz[0]*LayerThickness_local*dsdz*(1.+pow(x_k_minus[3],2));

    Q_k_minus[4][0]=dzscatter*sigmaeloss2*double(DeltaPlane)*buvz[1]*LayerThickness_local*dsdz*(1.+pow(x_k_minus[2],2));
    Q_k_minus[4][1]=dzscatter*sigmaeloss2*double(DeltaPlane)*buvz[0]*LayerThickness_local*dsdz*(1.+pow(x_k_minus[3],2));
    Q_k_minus[4][2]=sigmaeloss2*double(DeltaPlane)*buvz[1]*LayerThickness_local*dsdz*(1.+pow(x_k_minus[2],2));
    Q_k_minus[4][3]=sigmaeloss2*double(DeltaPlane)*buvz[0]*LayerThickness_local*dsdz*(1.+pow(x_k_minus[3],2));
    Q_k_minus[4][4]=sigmaeloss2;
  }
 
  // Display
  if(debug_new) {
      cout<<"------------------------------------------------------------"<<endl;
    cout << "1e6 * Q_k_minus" << endl;
    for(int i=0; i<5; ++i) {
      for(int j=0; j<5; ++j) {
        cout << 1e6*Q_k_minus[i][j] << " ";
      }
      cout << endl;
    }
    cout<<"------------------------------------------------------------"<<endl;
  }

}


void InoTrackFitAlg::ExtrapCovMatrix() {
  // C_k_intermediate = (F_k_minus * C_k_minus * F_k_minus^T) + Q_k_minus
 if(debug_fit)  cout <<" ----------------------InoTrackFitAlg : ExtrapCovMatrix----------------------" << endl;

  for (int i=0; i<5; ++i) {
    for (int j=0; j<5; ++j) {  
      C_k_intermediate[i][j]=0;
      for (int l=0; l<5; ++l) {
        for (int m=0; m<5; ++m) {   
          C_k_intermediate[i][j]+=F_k_minus[i][m]*C_k_minus[m][l]*F_k_minus[j][l];
        }
      }
      C_k_intermediate[i][j]+=Q_k_minus[i][j];
    }
  }
  
  // Diagonal elements should be positive
  double covlim = 1.e-8;

  for(int i=0; i<5; ++i) {
    if(C_k_intermediate[i][i]<covlim) {
      //GMA reopen this      cout <<" InoTrackFitAlg : Negative diagonal element in C_k_intermediate" << endl;
      C_k_intermediate[i][i]=covlim;
    }
  }
  
  // Display
  if(debug_fit) {
      cout<<"------------------------------------------------------------"<<endl;
    cout << "C_k_intermediate" << endl;
    for(int i=0; i<5; ++i) {
      for(int j=0; j<5; ++j) {
        cout << C_k_intermediate[i][j] << " ";
      }
      cout << endl;
    }
    cout<<"------------------------------------------------------------"<<endl;
  }
}

void InoTrackFitAlg::CalcKalmanGain(const int NewPlane)
{
  // K_k = C_k_intermediate * H_k^T * ( V_k + H_k * C_k_intermediate * H_k^T )^-1
  //                                            ==A1_k
  //                                    =    A2_k
  //                                 =  B2_k; 
 if(debug_fit)    cout <<" InoTrackFitAlg : CalcKalmanGain" << endl;

  double A1_k[2][2];
  double B1_k[2][2];
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      A1_k[i][j] = B1_k[i][j] = 0;
    }
  }

  // H_k has only one non-zero element, so we can reduce matrix multiplication required
  int PlaneView = TrkClustsData[NewPlane][0].PlaneView;
  
  if(PlaneView%2==0) {
    A1_k[0][0] =C_k_intermediate[0][0];
    A1_k[0][0]+=TrkClustsData[NewPlane][0].XPosErrSq;  // Add uncertainty in measurement
  }
  //    if ( Denominator==0)  Denominator=1.e-12;
  if (PlaneView >0) {
    A1_k[1][1] = C_k_intermediate[1][1];
    A1_k[1][1]+= TrkClustsData[NewPlane][0].YPosErrSq;  // Add uncertainty in measurement
  }

  if (PlaneView==2) {
    A1_k[0][1] = C_k_intermediate[0][1];
    A1_k[1][0] = C_k_intermediate[1][0];
  }

  double determinant = A1_k[0][0]*A1_k[1][1] - A1_k[0][1]*A1_k[1][0];
  //Inverse matrix;
  //  cout <<" InoTrackFitAlg : V_k " << NewPlane<<" "<<TrkClustsData[NewPlane][0].XPosErrSq <<" " <<TrkClustsData[NewPlane][0].YPosErrSq<<" "<<PlaneView<<" "<<determinant<<" "<<A1_k[0][0]<<" "<<A1_k[0][1]<<" "<<A1_k[1][0]<<" "<<A1_k[1][1]<<endl;
  if (determinant !=0) {
    B1_k[0][0] =  A1_k[1][1]/determinant;
    B1_k[0][1] = -A1_k[0][1]/determinant;
    B1_k[1][0] = -A1_k[1][0]/determinant;
    B1_k[1][1] =  A1_k[0][0]/determinant;
  
    double C1_k[5][2];
    for (int i=0; i<5; ++i) {
      for (int j =0; j<2; j++) {
	K_k[i][j]=0;
	for (int m=0; m<5; ++m) {
	  C1_k[m][j] = 0;
	  for (int n=0; n<2; n++) {
	    //	    C1_k[m][j] +=H_k[n][m]*B1_k[n][j];
	    K_k[i][j] +=C_k_intermediate[i][m]*H_k[n][m]*B1_k[n][j];
	  }
	  //	  K_k[i][j] +=C_k_intermediate[i][m]*C1_k[m][j];
	}
      }
    }

    //    cout <<" InoTrackFitAlg : Kalman Gain: "<<NewPlane<<" "
    //	 << K_k[0][0] << " " << K_k[1][0] << " " << K_k[2][0] << " "
    //	 << K_k[3][0] << " " << K_k[4][0] << " "
    //	 << K_k[0][1] << " " << K_k[1][1] << " " << K_k[2][1] << " "
    //	 << K_k[3][1] << " " << K_k[4][1] << endl;

  } else { 
 //   cout <<" InoTrackFitAlg : V_k + (H_k * C_k_intermediate * H_k_transpose) is zero!" << endl;
  }
}



void InoTrackFitAlg::UpdateStateVector(const int Plane, const int NewPlane, const bool GoForward)
{  
  // x_k = (F_k_minus * x_k_minus) + K_k(m_k - (H_k * F_k_minus * x_k_minus) )
  
  if(TrkFitterDebug>10) {
    cout <<" ----------------------InoTrackFitAlg : UpdateStateVector---------------------- "<<Plane<<" , "<<NewPlane << endl;
  }
  int nswimfail=0;


  // Calculate F_k_minus * x_k_minus, using the Swimmer
  // Also get an accurate measure of dS and Range from the Swimmer
  double StateVector[6];
  double Prediction[6];
  bool GetPrediction=false;

  for(int i=0; i<6; ++i) {StateVector[i]=x_k_minus[i];}

  if(TrkFitterDebug>10) {
    cout <<" InoTrackFitAlg : state: Plane " << Plane<<" NewPlane "
         << NewPlane<<" SV: "<<StateVector[0] << " " 
         << StateVector[1] << " " << StateVector[2] << " " 
         << StateVector[3] << " " << StateVector[4] <<" " << StateVector[5]<<" phi= "
         << atan2(StateVector[3],StateVector[2])<<" theta= "
      //	 << 1./pow(1+pow(StateVector[2],2.)+pow(StateVector[3],2.),0.5)<<" "
         << acos(1./pow(1+pow(StateVector[2],2.)+pow(StateVector[3],2.),0.5))<<" mom= ";
  
    if (StateVector[4] !=0) {
      cout <<1./StateVector[4] <<endl;
    } else {
      cout <<StateVector[4]<<endl;
    }
  }
  // Do the swim
  while(GetPrediction==false && nswimfail<=10) {
    
    if(TrkFitterDebug>10) {
      cout <<" 1122 : state " << Plane<<" "
	   << NewPlane<<" "<<StateVector[0] << " " 
	   << StateVector[1] << " " << StateVector[2] << " " 
	   << StateVector[3] << " " << StateVector[4] <<" " << StateVector[5]<<" phi= "
	   << atan2(StateVector[3],StateVector[2])*180/3.14<<" theta= "
	//	 << 1./pow(1+pow(StateVector[2],2.)+pow(StateVector[3],2.),0.5)<<" "
	   << acos(1./pow(1+pow(StateVector[2],2.)+pow(StateVector[3],2.),0.5))*180/3.14<<" mom= ";
	
      if (StateVector[4] !=0) {
	cout <<1./StateVector[4] <<endl;
      } else {
	cout <<StateVector[4]<<endl;
      }
    }
  
    GetPrediction=Swim(StateVector, Prediction, Plane, NewPlane, GoForward);
    
    //GMA do not allow starting momentum to be very large
    /*
    while (x_k_minus[4]==0 && fabs(Prediction[4])<1.e-2) {
      cout <<"x_4 "<<x_k_minus[4]<<" "<<StateVector[4]<<" "<<Prediction[4]<<endl;
      StateVector[4] = 10.*Prediction[4];
      if (StateVector[4] >1.e-2) StateVector[4] =1.e-2;
      if (StateVector[4] <-1.e-2) StateVector[4] =-1.e-2;
      GetPrediction=Swim(StateVector, Prediction, Plane, NewPlane, GoForward);
    }
    */

    if(TrkFitterDebug>10) {
      cout <<"2233  : State predict " << Plane<<" "
	   << NewPlane<<" SV: "<< Prediction[0] << " " 
	   << Prediction[1] << " " << Prediction[2] << " " 
	   << Prediction[3] << " " << Prediction[4] <<" phi= "
	   << atan2(Prediction[3],Prediction[2])*180/3.14<<" theta= "
	   << acos(1./pow(1+pow(Prediction[2],2.)+pow(Prediction[3],2.),0.5))*180/3.14<<" mom = ";
      if (Prediction[4] !=0) {
	cout <<1./Prediction[4] <<endl;
      } else {
	cout <<Prediction[4]<<endl;
      }
    }

      predicted_k[0] = Prediction[0];
      predicted_k[1] = Prediction[1];
      predicted_k[2] = Prediction[2];
      predicted_k[3] = Prediction[3];
      predicted_k[4] = Prediction[4];
      predicted_k[5] = Prediction[5];

      
    if(GetPrediction==false) {
      StateVector[4]*=0.5; 
      nswimfail++; TotalNSwimFail++;
      // cout <<"  : State predict " << Plane<<" "
      // 	   << NewPlane<<" SV: "<< Prediction[0] << " " 
      // 	   << Prediction[1] << " " << Prediction[2] << " " 
      // 	   << Prediction[3] << " " << Prediction[4] << " " << 1/Prediction[4]
      // 	   <<endl;
      // cout <<" InoTrackFitAlg : UpdateStateVector, Prediction failed - Double momentum and swim again" << endl;
    }
  }

  if(nswimfail>10) {  // Swim shouldn't fail, as it succeeded to get the propagator
    //GMA reopen this    cout <<" InoTrackFitAlg : UpdateStateVector, nswimfail>10, fail track" << endl;
    cout<<"nswimfail "<<nswimfail <<"  : State predict " << Plane<<" "
	<< NewPlane<<" SV: "<< Prediction[0] << " " 
	<< Prediction[1] << " " << Prediction[2] << " " 
	<< Prediction[3] << " " << Prediction[4] << " " << 1/Prediction[4]
	<<endl;
    PassTrack=false;
  }
  
  //  cout <<" InoTrackFitAlg : UpdateStateVector, Check predicted state " << endl;
  CheckValues(Prediction, NewPlane);

  // Calculate H_k * F_k_minus * x_k_minus
  int PlaneView = TrkClustsData[NewPlane][0].PlaneView;
  for (int i=0; i<2; i++) {
    for (int j=0; j<5; j++) {
      H_k[i][j]=0;
    }
  }
  // Calculate A1_k = m_k - H_k * F_k_minus * x_k_minus
  double A1_k[2]={0,0};
  if (PlaneView%2==0) {H_k[0][0]=1; A1_k[0] = TrkClustsData[NewPlane][0].XPos;}
  if (PlaneView   >0) {H_k[1][1]=1; A1_k[1] = TrkClustsData[NewPlane][0].YPos;}
  
  //  if ((Plane >=78 && Plane <=80)||(Plane >=126 && Plane <=128)) {
  if(debug_fit) cout <<"a1_k "<< PlaneView<<" "<<A1_k[0]<<" "<<A1_k[1]<<endl;
  //  }

  for (int i=0; i<2; i++) {
    for (int j=0; j<5; j++) {
      A1_k[i] -= H_k[i][j]*Prediction[j];
      //      if ((Plane >=78 && Plane <=80)||(Plane >=126 && Plane <=128)) {
      if(debug_fit) cout <<"A1_k "<<i<<" "<<j<<" "<<A1_k[i]<<" "<<H_k[i][j]<<" "<<Prediction[j]<<endl;
      //      }
    }
  }

  // Calculate x_k
  x_k[5]=Prediction[5];
  for (int i=0; i<5; ++i) {
    x_k[i]=Prediction[i];
    for (int j=0; j<2; j++) {
      x_k[i]+=K_k[i][j]*A1_k[j];
      //      if ((Plane >=78 && Plane <=80)||(Plane >=126 && Plane <=128)) {
          if(debug_fit) cout <<"x_k "<<i <<" "<<j<<" "<<Prediction[i]<<" "<<x_k[i]<<" "<<K_k[i][j]*A1_k[j]<<endl;
      //      }
    }
  }
  
  //  cout <<" InoTrackFitAlg : UpdateStateVector, Check filtered state " << endl;
   CheckValues(x_k, NewPlane);


  // GMA Care with multiple range corrections - do not want to flip sign
  // (multiple corrections mean sign changes can occur even though absolute value stays same)
  // JAM up to 40 from 4

  double Maxqp = 4.;
  if(fabs(x_k_minus[4])==Maxqp && 
     ( (GoForward==true && ZIncreasesWithTime==true) 
       || (GoForward==false && ZIncreasesWithTime==false) ) ) 
    {
      if(!LastIteration) x_k[4] = (x_k_minus[4]>0 ? Maxqp : -Maxqp);
      //      cout << " resetting in UpdateStateVector " <<  endl;
    }

  //if on last plane in forward swim, disregard sign flip
  if(x_k_minus[4]!=0){
    if ( x_k[4]/x_k_minus[4]<0 && 
	 ( (GoForward==true && ZIncreasesWithTime==true && NewPlane >= EndofRangePlane ) 
	   || (GoForward==false && ZIncreasesWithTime==false && NewPlane <= EndofRangePlane ) )) 
      {
        x_k[4] = -x_k[4];
      }
  }

  x_k[5] = 1; //GMA 09/02/2009  for all upgrade



  // Display
  if(TrkFitterDebug>10) {
    cout <<"  : Filtered:  " << Plane<<"  "
	 << NewPlane<<" ClustDataXY "
	 <<TrkClustsData[NewPlane][0].XPos<<" "
	 <<TrkClustsData[NewPlane][0].YPos<<endl;
    cout<<" SV: "
	<< x_k[0] << " " << x_k[1] << " " << x_k[2] << " "
	<< x_k[3] << " " << x_k[4] <<" ph,th= "
	<< atan2(x_k[3],x_k[2])*180/3.14<<" "
	<< acos(1./pow(1+pow(x_k[2],2.)+pow(x_k[3],2.),0.5))*180/3.14<<" mom= ";
    if (x_k[4] !=0) {
      cout <<1./x_k[4] <<endl;
    } else {
      cout <<x_k[4]<<endl;
    }
  }
}



void InoTrackFitAlg::UpdateCovMatrix()
{
  // C_k = (Identity - (K_k * H_k) ) * C_k_intermediate
  if(debug_fit)
    cout <<"---------------------- InoTrackFitAlg : UpdateCovMatrix----------------------" << endl;
  double A1_k[5][5];

  for (int i=0; i<5; ++i) {
    for (int j=0; j<5; ++j) {  
      C_k[i][j]=0;
      for (int m=0; m<5; ++m) {   
	A1_k[i][m] = 0;
	for (int n = 0; n<2; n++) {
	  A1_k[i][m] += K_k[i][n]*H_k[n][m];
	}
        C_k[i][j]+=(Identity[i][m] - A1_k[i][m]) * C_k_intermediate[m][j];
      }
    }
  }

  // Diagonal elements should be positive
  double covlim = 1.e-8;

  for(int i=0; i<5; ++i) {
    if(C_k[i][i]<covlim) {
      //GMA reopen this      cout <<" InoTrackFitAlg : Negative diagonal element in C_k" << endl;
      C_k[i][i]=covlim;
    }
  }
  
  // Display
  if(debug_fit) {
      cout<<"------------------------------------------------------------"<<endl;
    //cout << "Filtered Covariance matrix" << endl;
    for(int i=0; i<5; ++i) {
      for(int j=0; j<5; ++j) {
       cout << C_k[i][j] << " ";
      }
      cout << endl;  
    }
    cout<<"------------------------------------------------------------"<<endl;
  }


}



void InoTrackFitAlg::MoveArrays()
{
  // Move k to k-1 ready to consider next clust

 if(debug_fit)  cout <<"----------------------InoTrackFitAlg : MoveArrays 1----------------------" << endl;
  for (int i=0; i<5; ++i) {
    for (int j=0; j<5; ++j) { 
      C_k_minus[i][j]=C_k[i][j];
    }
  }
  
  for (int l=0; l<6; ++l) {
    x_k_minus[l]=x_k[l];
  }
 
 if(debug_fit)  cout <<"----------------------InoTrackFitAlg : MoveArrays end 1----------------------" << endl;
}



void InoTrackFitAlg::CheckValues(double* Input, const int NewPlane)
{
  // Make range and gradient corrections
  // Possible source of offset in q/p resolutions
  if(debug_fit)  cout<< " Check value " <<endl;
  // Range check
    
  double Maxqp=4.; double Maxqpfrac=1.2;
  double Range=fTrackCand->GetRange(NewPlane);

  //JAM signal end of range found
  if(fabs(Input[4])>10.0) EndofRange=true; 

   if(Range>0. && (Maxqpfrac*doubleLa/Range)<Maxqp) {Maxqp=(Maxqpfrac*doubleLa/Range);}
     //cout <<" InoTrackFitAlg :  Range " << Range << " Maxqp " << Maxqp << endl;

  if(LastIteration) Maxqp=40;
  if(fabs(Input[4])>Maxqp){
       //cout << " CheckValues: Range check correction " << Input[4] << " " << Maxqp << endl;
    // GMA reopen this    cout <<" InoTrackFitAlg : CheckValues, Range check correction" << endl;
    Input[4]=(Input[4]>0 ? Maxqp : -Maxqp);
  }
  
  // Gradient check
  double Maxgradient=25.;

  if(fabs(Input[2])>Maxgradient) {
    //GMA reopen this    cout <<" InoTrackFitAlg : CheckValues, Gradient correction, U" << endl;
    Input[2]=(Input[2]>0 ? Maxgradient : -Maxgradient);
  }

  if(fabs(Input[3])>Maxgradient) {
    // GMA reopen this    cout <<" InoTrackFitAlg : CheckValues, Gradient correction, V" << endl;
    Input[3]=(Input[3]>0 ? Maxgradient : -Maxgradient);
  }
  
  // Check u and v values are not rubbish
  if(fabs(Input[0])<5000. && fabs(Input[1])<5000.) {PassTrack=true;}
  else { cout << Input[0] << " "<<Input[1] <<endl;
      cout<< " Passtrack 1 "<<" ievt "<<pAnalysis->ievt2<<endl;
      PassTrack=false;}              //AAR: This is set false after the change in swimswimmer.
}

void InoTrackFitAlg::StoreFilteredData(const int NewPlane)
{
  // Store the data required for matching Kalman output data to strips
  if(debug_new) 
    cout <<"----------------------InoTrackFitAlg : StoreFilteredData----------------------"<<NewPlane<<" "<< FilteredData[NewPlane].size() << endl;
  
  for (unsigned ij=0; ij<FilteredData[NewPlane].size(); ij++) {
    if (FilteredData[NewPlane][ij].x_k5==0) {
      FilteredData[NewPlane].erase(FilteredData[NewPlane].begin()+ij); ij-- ;
    }
  }

  for (unsigned ij=0; ij<ExtraPolData[NewPlane].size(); ij++) {
    if (ExtraPolData[NewPlane][ij].x_k5==0) {
      ExtraPolData[NewPlane].erase(ExtraPolData[NewPlane].begin()+ij); ij--;
    }
  }
  
  FiltDataStruct pretemp;

  pretemp.x_k0 = predicted_k[0];
  pretemp.x_k0 = predicted_k[1];
  pretemp.x_k0 = predicted_k[2];
  pretemp.x_k0 = predicted_k[3];
  pretemp.x_k0 = predicted_k[4];
  pretemp.x_k0 = int(predicted_k[5]);
  pretemp.x_k6 = false;

  ExtraPolData[NewPlane].push_back(pretemp);
  
  FiltDataStruct temp;
  // cout<<"Storing "<<NewPlane<<" "<<x_k[0]<<" "<<x_k[1]<<" "<<x_k[2]<<" "<<x_k[3]<<" "<<x_k[4]<<endl;
  temp.x_k0=x_k[0]; temp.x_k1=x_k[1];
  temp.x_k2=x_k[2]; temp.x_k3=x_k[3];
  temp.x_k4=x_k[4]; temp.x_k5=int(x_k[5]);
  temp.x_k6=true;
  //  FilteredData[NewPlane].clear();

  FilteredData[NewPlane].push_back(temp);
  // cout<<".... FilteredData["<<NewPlane<<"].size() "<< FilteredData[NewPlane].size() << endl;
}

void InoTrackFitAlg::StoreFilteredData_sr(const int NewPlane, double* prediction, bool str)
{
  // Store the data required for matching Kalman output data to strips
  if(debug_fit)  
    cout <<"----------------------InoTrackFitAlg : StoreFilteredData_sr----------------------  " <<NewPlane<<" "<<str<<" "<<FilteredData[NewPlane].size()<< endl;
  for (unsigned ij=0; ij<ExtraPolData[NewPlane].size(); ij++) {
    if (ExtraPolData[NewPlane][ij].x_k5==0) {
      ExtraPolData[NewPlane].erase(ExtraPolData[NewPlane].begin()+ij); ij-- ;
    }
  }
  
  for (unsigned ij=0; ij<FilteredData[NewPlane].size(); ij++) {
    if (FilteredData[NewPlane][ij].x_k5==0) {
      FilteredData[NewPlane].erase(FilteredData[NewPlane].begin()+ij); ij-- ;
    }
  }
  
  FiltDataStruct temp;
  
  temp.x_k0=prediction[0]; temp.x_k1=prediction[1];
  temp.x_k2=prediction[2]; temp.x_k3=prediction[3];
  temp.x_k4=prediction[4]; temp.x_k5=0;
  temp.x_k6=str;

  for(int ijx1=0; ijx1<5; ijx1++) {
    for(int jkx1=0; jkx1<5; jkx1++) {
      temp.C_k[ijx1][jkx1] = C_k[ijx1][jkx1];
    }
  }
  // cout<<"Storing "<<NewPlane<<" "<<prediction[0]<<" "<<prediction[1]<<" "<<prediction[2]<<" "<<prediction[3]<<" "<<prediction[4]<<endl;

  //  FilteredData[NewPlane].clear();
  FilteredData[NewPlane].push_back(temp);
  ExtraPolData[NewPlane].push_back(temp);
  if(debug_fit)  cout <<"----------------------InoTrackFitAlg : StoreFilteredData_sr end----------------------" << endl;

}



void InoTrackFitAlg::SetTrackProperties(double* Prediction) { 

  if(nbfield>0) bave /=nbfield;

  // Carry out the assignment of variables to the new fitted track 
  if(TrkFitterDebug>10) {
    cout <<"----------------------InoTrackFitAlg : SetTrackProperties----------------------" <<  nbfield <<" " << bave << endl;
  }  
  
  // Momentum, charge and error on q/p
  
  if(x_k[4]!=0.) {double xx =fabs(1./x_k[4]); fTrackCand->SetMomentumCurve(xx);}
  fTrackCand->SetRangeBiasedQP(x_k4_biased);
  if(x_k[4]>0.) {fTrackCand->SetEMCharge(1.);}
  else if(x_k[4]<0.) {fTrackCand->SetEMCharge(-1.);}
  
  double xx = pow(VtxCov[4],0.5);
  fTrackCand->SetVtxQPError(xx);
    
  // Positions and angles  
  int VtxPlane;
  int EndPlane;
  double dsdz;

  if(ZIncreasesWithTime==true) {VtxPlane=MinPlane; EndPlane=MaxPlane;}
  else {VtxPlane=MaxPlane; EndPlane=MinPlane;}
  
  if(debug_new) {
    for(unsigned int kpc=0; kpc<fTrackCand->ClustsInTrack.size(); kpc++) {
      cout<<"clust "<<kpc<<" = ("<<fTrackCand->ClustsInTrack[kpc]->GetXPos()<<","<<fTrackCand->ClustsInTrack[kpc]->GetYPos()<<","<<fTrackCand->ClustsInTrack[kpc]->GetZPlane()<<");"<<endl;
    }
  }

  // cout<<"VtxPlane "<<VtxPlane<<" "<<FilteredData[VtxPlane].size()<<endl;
  if(FilteredData[VtxPlane].size()>0) {
    // Vtx and end coordinates  //GMA 120809  0 or next one ??????
    fTrackCand->SetVtxU(FilteredData[VtxPlane][0].x_k0);
    fTrackCand->SetVtxV(FilteredData[VtxPlane][0].x_k1);
    fTrackCand->SetVtxZ(ZPosLayer[VtxPlane]); //SlcClustData[VtxPlane][0].csh->GetZPos());
    fTrackCand->SetVtxPlane(VtxPlane);
    
    dsdz=pow(1.+pow(FilteredData[VtxPlane][0].x_k2,2)+pow(FilteredData[VtxPlane][0].x_k3,2),0.5);
    if(ZIncreasesWithTime==false) {dsdz=-dsdz;}
    fTrackCand->SetVtxDirCosU(FilteredData[VtxPlane][0].x_k2/dsdz);
    fTrackCand->SetVtxDirCosV(FilteredData[VtxPlane][0].x_k3/dsdz);
    fTrackCand->SetVtxDirCosZ(1./dsdz);
    // cout<<"mom vtx "<<FilteredData[VtxPlane][0].x_k4<<" "<<FilteredData[EndPlane][0].x_k1<<" "<<FilteredData[EndPlane][0].x_k2<<endl;
    if (FilteredData[VtxPlane][0].x_k4 !=0) {
      fTrackCand->SetMomentumCurve(1./FilteredData[VtxPlane][0].x_k4);
    } else {
      fTrackCand->SetMomentumCurve(-10000000.);
    }

    fTrackCand->SetDirCosU(FilteredData[VtxPlane][0].x_k2/dsdz);
    fTrackCand->SetDirCosV(FilteredData[VtxPlane][0].x_k3/dsdz);
    fTrackCand->SetDirCosZ(1./dsdz);
  } else {
    PassTrack=false;
    cout<<" Vtx Plane empty." << endl;
  }
  // cout<<"PassTrack Vtx "<<PassTrack<<endl;
  // cout<<"EndPlane "<<EndPlane<<" "<<FilteredData[EndPlane].size()<<endl;
  if(FilteredData[EndPlane].size()>0) {
    fTrackCand->SetEndU(FilteredData[EndPlane][0].x_k0);
    fTrackCand->SetEndV(FilteredData[EndPlane][0].x_k1);
    fTrackCand->SetEndZ(ZPosLayer[EndPlane]); //SlcClustData[EndPlane][0].csh->GetZPos());
    fTrackCand->SetEndPlane(EndPlane);
    
    dsdz=pow(1.+pow(FilteredData[EndPlane][0].x_k2,2)+pow(FilteredData[EndPlane][0].x_k3,2),0.5);
    if(ZIncreasesWithTime==false) {dsdz=-dsdz;}
    fTrackCand->SetEndDirCosU(FilteredData[EndPlane][0].x_k2/dsdz);
    fTrackCand->SetEndDirCosV(FilteredData[EndPlane][0].x_k3/dsdz);
    fTrackCand->SetEndDirCosZ(1./dsdz);
    // cout<<"mom end "<<FilteredData[EndPlane][0].x_k4<<" "<<FilteredData[EndPlane][0].x_k0<<" "<<FilteredData[EndPlane][0].x_k1<<" "<<endl;
    if (FilteredData[EndPlane][0].x_k4 !=0) {
      fTrackCand->SetEndMomentumCurve(1./FilteredData[EndPlane][0].x_k4);
    } else {
      fTrackCand->SetEndMomentumCurve(-10000000.);
    }
    fTrackCand->SetEndQP(FilteredData[EndPlane][0].x_k4);
  }  else {
    PassTrack=false;
    cout<<" End Plane empty." << endl;
  }
  // cout<<"PassTrack "<<PassTrack<<endl;
  int newfout = iternew;
  fTrackCand->SetNewFitOut(newfout);
  
  if(newfit_x_k[4] != 0) {
    fTrackCand->SetNewMomentum(1/newfit_x_k[4]);
  } else {
    fTrackCand->SetNewMomentum(-10000000.);
  }
  
  // cout<<"input "<<1./dsdz<<" "<<fTrackCand->GetDirCosZ()<<endl;
  
  // Errors on vtx positions and angles  
  fTrackCand->SetVtxUError(pow(VtxCov[0],0.5));
  fTrackCand->SetVtxVError(pow(VtxCov[1],0.5));
  fTrackCand->SetVtxdUError(pow(VtxCov[2],0.5));
  fTrackCand->SetVtxdVError(pow(VtxCov[3],0.5)); 
  
  // Errors on end positions, angles and q/p
  fTrackCand->SetEndUError(pow(EndCov[0],0.5));
  fTrackCand->SetEndVError(pow(EndCov[1],0.5));
  fTrackCand->SetEnddUError(pow(EndCov[2],0.5));
  fTrackCand->SetEnddVError(pow(EndCov[3],0.5)); 
  fTrackCand->SetEndQPError(pow(EndCov[4],0.5)); 
  
  fTrackCand->SetBave(bave);
  
  // Momentum and thets phi after extrapolation of half length;
  dsdz=pow(1.+pow(Prediction[2],2)+pow(Prediction[3],2),0.5);
  if(ZIncreasesWithTime==false) {dsdz=-dsdz;}
  
  if (Prediction[4] !=0) {
    fTrackCand->SetMomentum(1./Prediction[4]);
  } else {
    fTrackCand->SetMomentum(-10000000.);
  }
  
  fTrackCand->SetTheta(acos(1./dsdz));
  fTrackCand->SetPhi(atan2(Prediction[3]/dsdz,Prediction[2]/dsdz));
  
  // cout<<"Prediction "<<1./Prediction[4]<<" "<<1./dsdz<<" "<<atan2(Prediction[3]/dsdz,Prediction[2]/dsdz)*180/3.14<<endl;

  // if (Prediction[4] !=0) {
  //   double modp_set = fabs(1/Prediction[4]);
  //   if(ZIncreasesWithTime==false) {modp_set=-modp_set;}
  //   double charge = 0.0;
  //   if(Prediction[4]>0.) charge = 1.;
  //   else if(Prediction[4]<0.) charge = -1.;
  
  //   dsdz=pow(1.+pow(Prediction[2],2)+pow(Prediction[3],2),0.5);
  
  //   TVector3 momvec_set(modp_set*(Prediction[2]/dsdz),
  // 			modp_set*(Prediction[3]/dsdz),
  // 			modp_set/dsdz);
  
  //   fTrackCand->SetMomentum(charge*momvec_set.Mag());
  //   fTrackCand->SetTheta(momvec_set.Theta());
  //   fTrackCand->SetPhi(momvec_set.Phi());
  // } else {
  //   dsdz=pow(1.+pow(Prediction[2],2)+pow(Prediction[3],2),0.5);
  //   fTrackCand->SetMomentum(-10000000.);
  //   fTrackCand->SetTheta(20);
  //   fTrackCand->SetPhi(20);
  // }
  
  fTrackCand->SetFCPC(FCorPC);

  // cout<<"SetTrackProp : theta : atan2 : phi "<<180*fTrackCand->GetTheta()/3.14<<" : "<<180*atan2(Prediction[3],Prediction[2])/3.14<<" : "<<180*fTrackCand->GetPhi()/3.14<<endl;
  
  if(fTrackCand->GetFCPC()<0 || fTrackCand->GetFCPC()>255) {
    cout<<"Error Error Error Error Error.... 222 "<<FCorPC<<" "<<fTrackCand->GetFCPC()<<endl;
  }

  // More variables to be set, in order to ensure compatibility
  fTrackCand->SetNTrackStrip(fFinderTrack->ClustsInTrack.size());
  
  fTrackCand->SetNIterate(NIter);
  fTrackCand->SetNSwimFail(TotalNSwimFail);
  
  
  // Obtain "fitting data" for the final track strips
  //  for (unsigned i=0; i<nLayer; ++i) {TrkClustsData[i].clear();} 
  //  GetFitData(MinPlane,MaxPlane);
  
  
  // Set tpos error and Calculate chi2, NDOF
  double Chi2=0; double Chi2Contrib=0; int NDOF=0; double FilteredXPos=0; double FilteredYPos=0;
  
  double momdS=0; double momRange=0;
  
  double sxy=0; // y (distance) = c t + shift
  double sx=0; // y = distance
  double sy=0; //   x = time
  double sn=0;
  double sx2=0;

    
  for (unsigned ijk=0; ijk<fTrackCand->GetEntries(); ijk++) {
    int i = fTrackCand->ClustsInTrack[ijk]->GetZPlane();
    //  for(int i=MinPlane; i<=MaxPlane; ++i) {
    if (i <=int(nLayer)) {
      if(TrkClustsData[i].size()>0) {
	
	if(TrkClustsData[i][0].XPosErrSq>0. || TrkClustsData[i][0].YPosErrSq>0.) {
	  fTrackCand->SetTrackPointXError(i,pow(TrkClustsData[i][0].XPosErrSq,0.5));
	  fTrackCand->SetTrackPointYError(i,pow(TrkClustsData[i][0].YPosErrSq,0.5));

	  
	  momdS += max(0., double(fTrackCand->GetdS(i)));
	  momRange += max(0.,double(fTrackCand->GetRange(i)));
	  
	  //	  if ((!ZIncreasesWithTime && ijk>0) || (ZIncreasesWithTime && ijk <fTrackCand->GetEntries()-1)) {
	  sn +=1;
	  sx += TrkClustsData[i][0].cltime; //Look again for return track
	  sy += momdS;
	  sxy = momdS*TrkClustsData[i][0].cltime;
	  sx2 = (TrkClustsData[i][0].cltime)*(TrkClustsData[i][0].cltime);
	  //	  }
	  
	  Chi2Contrib = 0;
	  
	  for (unsigned jk=0; jk<TrkClustsData[i].size() ;jk++) {
	    if(TrkClustsData[i][jk].PlaneView%2==0) {
	      FilteredXPos=FilteredData[i][jk].x_k0;
	      Chi2Contrib += pow((TrkClustsData[i][jk].XPos-FilteredXPos),2)/TrkClustsData[i][jk].XPosErrSq;
	      NDOF++;
	    }
	    
	    if (TrkClustsData[i][jk].PlaneView >=1) {
	      FilteredYPos=FilteredData[i][jk].x_k1;
	      Chi2Contrib += pow((TrkClustsData[i][jk].YPos-FilteredYPos),2)/TrkClustsData[i][jk].YPosErrSq;
	      NDOF++;
	    }
	  }
	  
	  //	} else if(TrkClustsData[i][0].PlaneView==2) {
	  //	  FilteredXPos=FilteredData[i][0].x_k0;
	  //	  FilteredYPos=FilteredData[i][0].x_k1;
	  //	  Chi2Contrib = pow((TrkClustsData[i][0].XPos-FilteredXPos),2) + 
	  //	    pow((TrkClustsData[i][0].YPos-FilteredYPos),2);
	  //	}
	  // TrkClustsData[i][0].TPosError;
	  fTrackCand->SetPlaneChi2(i,Chi2Contrib);        
	  
	  Chi2+=Chi2Contrib;
	  
          //     NDOF++;
	}
      }
    } else {
      if (TrkClustsData[i+shiftLa].size()>0) {
	momdS += fTrackCand->Get2dS(i+shiftLa);
	momRange += fTrackCand->Get2Range(i+shiftLa);
	sn +=1;
	sx += TrkClustsData[i+shiftLa][0].cltime; //Look again for return track
	sy += momdS;
	sxy += momdS*TrkClustsData[i+shiftLa][0].cltime;
	sx2 += (TrkClustsData[i][0+shiftLa].cltime)*(TrkClustsData[i][0].cltime);
        if(pAnalysis->isXtermOut==1){
	  cout <<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
	  cout <<"ijk "<<ijk<<" "<<i<<" "<<momdS<<" "<<momRange<<" "<<endl;
	  cout <<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;	
	}//isXterm
      }
    }
  }

  
  double velocity = 0;
  if (sn >0 && (sx2*sn - sx*sx) !=0) {
    velocity = 10*(sxy*sn - sx*sy)/(sx2*sn - sx*sx); // x 10^8m/s
  }
  if(pAnalysis->isXtermOut==1){
    cout<<endl;  
    cout <<"chisq : Chi2 NDOF momdS momRange  velocity"<< endl;
    cout <<"chisq "<<Chi2 <<" "<<NDOF-5<<" "<<momdS<<" "<<momRange<<" "<<velocity<< endl;
    cout<<"---------------------------------------------------------------------------"<<endl;
  }//isXterm
  fTrackCand->SetChi2(Chi2);
  fTrackCand->SetNDOF(NDOF-5); // Number of constraints set to 5
  fTrackCand->SetMomentumdS(momdS);
  fTrackCand->SetMomentumRange(momRange);
  fTrackCand->Setcval(velocity);

  //  fTrackCand->SetMomentum(momRange);  

  // Assign U, V and q/p values
  // cout<<"MinPlane "<<MinPlane<<" "<<MaxPlane<<endl;
  for(int ij=MinPlane; ij<=MaxPlane; ++ij) {                   //asm_170311
    if(TrkFitterDebug>10) {cout<<"PlaneInfostrore     "<<ij<<" "<<FilteredData[ij].size()<<endl;}
    if(FilteredData[ij].size()>0) {
      // cout<<"par "<<FilteredData[ij][0].x_k0<<" "<<FilteredData[ij][0].x_k1<<" "<<FilteredData[ij][0].x_k4<<endl;
      fTrackCand->SetU(ij,FilteredData[ij][0].x_k0);
      fTrackCand->SetV(ij,FilteredData[ij][0].x_k1);
      fTrackCand->SetPlaneQP(ij,FilteredData[ij][0].x_k4);
      fTrackCand->xfitpos1.push_back(FilteredData[ij][0].x_k0);
      fTrackCand->yfitpos1.push_back(FilteredData[ij][0].x_k1);
      fTrackCand->zfitpos1.push_back(ZPosLayer[ij]);
      fTrackCand->zfitlay1.push_back(ij);
      fTrackCand->extrapolx0.push_back(ExtraPolData[ij][0].x_k0);
      fTrackCand->extrapolx1.push_back(ExtraPolData[ij][0].x_k1);
      if(ExtraPolData[ij][0].x_k4!=0) {fTrackCand->extrapolmom.push_back(1./ExtraPolData[ij][0].x_k4);}
      else {fTrackCand->extrapolmom.push_back(-1000000);}
      fTrackCand->filteredx2.push_back(FilteredData[ij][0].x_k2);
      fTrackCand->filteredx3.push_back(FilteredData[ij][0].x_k3);
      fTrackCand->filteredx4.push_back(FilteredData[ij][0].x_k4);
      double dsdz1=pow(1.+pow(FilteredData[ij][0].x_k2,2)+pow(FilteredData[ij][0].x_k3,2),0.5);
      if(ZIncreasesWithTime==false) {dsdz1=-dsdz1;}
      if(FilteredData[ij][0].x_k4!=0) {fTrackCand->filteredmom.push_back(1./FilteredData[ij][0].x_k4);}
      else {fTrackCand->filteredmom.push_back(-1000000);}
      fTrackCand->filteredthe.push_back(acos(1./dsdz1));
      fTrackCand->filteredphi.push_back(atan2(FilteredData[ij][0].x_k2/dsdz1,FilteredData[ij][0].x_k3/dsdz1));
      if(TrkFitterDebug>10) {
	cout<<"PlaneInfostrore pln,x,y "<<ij<<" "<<FilteredData[ij][0].x_k0<<" "<<FilteredData[ij][0].x_k1<<" q/p,p "<<FilteredData[ij][0].x_k4<<" "<<1./FilteredData[ij][0].x_k4<<" the,phi "<<180*acos(1./dsdz1)/3.14<<" "<<atan2(FilteredData[ij][0].x_k2/dsdz1,FilteredData[ij][0].x_k3/dsdz1)<<endl;
      }
    }
  }
  int wq1size = fTrackCand->filteredmom.size();
  for(unsigned wq1=1; wq1<wq1size-1 && wq1size>2; wq1++) {
    double tmpv2 = pow(0.5*(fTrackCand->filteredmom[wq1-1]+fTrackCand->filteredmom[wq1+1])-fTrackCand->filteredmom[wq1],2);
    fTrackCand->momvecdiff1.push_back(tmpv2);
  }
  
  vector<Hep3Vector> finderhitpos;
  vector<int> zlayin;
  double strtslope[2];
  double strtintercpt[2];
  double strtchisq[2];
  int strtnhits[2];
  double strtexpecpos[2];
  for(unsigned aq1=0; aq1<fTrackCand->ClustsInTrack.size(); aq1++) {
    Hep3Vector tmp3vect;
    tmp3vect.setX(fTrackCand->ClustsInTrack[aq1]->GetXPos());
    tmp3vect.setY(fTrackCand->ClustsInTrack[aq1]->GetYPos());
    tmp3vect.setZ(fTrackCand->ClustsInTrack[aq1]->GetZPos());
    finderhitpos.push_back(tmp3vect);
    int tmpplnz = fTrackCand->ClustsInTrack[aq1]->GetZPlane();
    zlayin.push_back(tmpplnz);
    if(ExtraPolData[tmpplnz].size()) {
      double tmpv2 = pow(pow(ExtraPolData[tmpplnz][0].x_k0-tmp3vect.x(),2)+pow(ExtraPolData[tmpplnz][0].x_k1-tmp3vect.y(),2),0.5);
      fTrackCand->radialdiff1.push_back(tmpv2);
    }
    if(TrkFitterDebug>10) {
      cout<<"aq1 "<<aq1<<" "<<fTrackCand->ClustsInTrack[aq1]->GetZPlane()<<" "<<tmp3vect<<endl;
    }
  }
  
  
  StraightLineFit(finderhitpos,zlayin,fTrackCand->GetEndPlane(),strtslope,strtintercpt,strtchisq,strtnhits,strtexpecpos);
  fTrackCand->SetStraightLineSlopeX(strtslope[0]);
  fTrackCand->SetStraightLineSlopeY(strtslope[1]);
  fTrackCand->SetStraightLineInterceptX(strtintercpt[0]);
  fTrackCand->SetStraightLineInterceptY(strtintercpt[1]);
  fTrackCand->SetStraightLineChi2X(strtchisq[0]);
  fTrackCand->SetStraightLineChi2Y(strtchisq[1]);
  fTrackCand->SetStraightLineNhitsX(strtnhits[0]);
  fTrackCand->SetStraightLineNhitsY(strtnhits[1]);
  fTrackCand->SetStraightLineXExpec(strtexpecpos[0]);
  fTrackCand->SetStraightLineYExpec(strtexpecpos[1]);

  double simplecurv,simpleradii,simplexexpec;
  double simplex0[2],simplechi2[3],simplexavg[4];
  int simplenhits;
  simple_track_fit(finderhitpos,fTrackCand->GetEndPlane(),simplecurv,simpleradii,simplex0,simplechi2,simplexavg,simplenhits,simplexexpec);

  fTrackCand->SetSimpleCurv(simplecurv);
  fTrackCand->SetSimpleRadii(simpleradii);
  fTrackCand->SetSimpleX0(simplex0[0]);
  fTrackCand->SetSimpleZ0(simplex0[1]);
  fTrackCand->SetSimpleChi2Pos(simplechi2[0]);
  fTrackCand->SetSimpleChi2Neg(simplechi2[1]);
  fTrackCand->SetSimpleAvgXPos(simplexavg[0]);
  fTrackCand->SetSimpleAvgXNeg(simplexavg[1]);
  fTrackCand->SetSimpleChi2Cndn(simplechi2[2]);
  fTrackCand->SetSimpleAvgXCndn(simplexavg[2]);
  fTrackCand->SetSimpleAvgXMeas(simplexavg[3]);
  fTrackCand->SetSimpleNhits(simplenhits);
  fTrackCand->SetSimpleXExpec(simplexexpec);
  
  // for(int dw1=0; dw1<2; dw1++) {
  //   cout<<"dw1 "<<dw1<<" "<<strtslope[dw1]<<" "<<strtintercpt[dw1]<<" "<<strtchisq[dw1]<<" "<<strtnhits[dw1]<<endl;
  // }
  
  if(timecheck) {
    double tmpslope=-1000;
    double tmpintercept = -10000;
    double tmpxtexp = -10000;
    bool updown = DirectionFromFitterHits(fTrackCand,EndPlane,tmpslope,tmpintercept,tmpxtexp);
    //yexp[ij]=tmpslope*xval[ij]+intersect;
    fTrackCand->SetNewTimeEndPlaneExp(tmpxtexp);
    fTrackCand->SetNewTimeSlope(tmpintercept);
    fTrackCand->SetNewTimeIntercept(tmpslope);
  }
  
  InoHit_Manager *pHitCollection = InoHit_Manager::APointer;
  
  for(int iw1=EndPlane-2; iw1>-1; iw1--) {
    for(unsigned qw1=0; qw1<pHitCollection->InoHit_list.size(); qw1++) {
      if(pHitCollection->InoHit_list[qw1]->GetZPlane()==iw1 && fabs(pHitCollection->InoHit_list[qw1]->GetTime()-fTrackCand->GetNewTimeEndPlaneExp())<15.0) {
	fTrackCand->HitsNotInTrack.push_back(pHitCollection->InoHit_list[qw1]);
      }
    }
  }

  int tmcnt1 = 0;
  int tmcnt2 = 0;
  for(unsigned qw1=0; qw1<pHitCollection->InoHit_list.size(); qw1++) {
    // cout<<"inloop "<<qw1<<" "<<pHitCollection->InoHit_list[qw1]->GetZPlane()<<" "<<EndPlane<<" "<<pHitCollection->InoHit_list[qw1]->GetTime()<<" "<<fTrackCand->GetNewTimeEndPlaneExp()<<endl;
    if(pHitCollection->InoHit_list[qw1]->GetZPlane()==EndPlane && abs(pHitCollection->InoHit_list[qw1]->GetTime()-fTrackCand->GetNewTimeEndPlaneExp())<15.0) {
      tmcnt1++;
    }
    if(pHitCollection->InoHit_list[qw1]->GetZPlane()==EndPlane-1 && abs(pHitCollection->InoHit_list[qw1]->GetTime()-fTrackCand->GetNewTimeEndPlaneExp())<15.0) {
      tmcnt2++;
    }
  }
  // cout<<tmcnt1<<" "<<tmcnt2<<endl;
  fTrackCand->SetNhitsEndPlane(tmcnt1);
  fTrackCand->SetNhitsEndPlaneM1(tmcnt2);
  
  // Fill time and range maps
  SetT() ;// &cth);
  //  SetRangeAnddS();//cth);
  
  // Set all time related variables
  TimingFit(); //cth);

  
  //  Calibrate(&cth);
  // cout<<"End.."<<endl;

}

void InoTrackFitAlg::TimingFit() //CandFitTrackCamHandle &cth)
{
  if(debug_fit)
    cout <<"----------------------InoTrackFitAlg : TimingFit----------------------" << endl;
  
  // Initialisations
  double s; double t; double q; int n=0; 
  double MinUncertainty = 0.; double MinCT=-3000.;

  // Time of first strip in track
  StripListTime=9.e10;

  // Create an offset such that dS=0 at the MinPlane
  double dSOffset=0.; double Sign=-1.; double dS[doubleLa]; // GMA need to put from db
  if(ZIncreasesWithTime==true) {dSOffset=fTrackCand->GetdS(MinPlane); Sign=1.;}             //asm_170311

  // Store data needed in arrays. Pulse is in PEs.
  double Qp[doubleLa];  double Qm[doubleLa];
  double CTp[doubleLa]; double CTm[doubleLa];
  int Skipp[doubleLa];  int Skipm[doubleLa];
  double C=3.e8;   

  double ErrorParam[3];
  ErrorParam[0]=0.; ErrorParam[1]=0.; ErrorParam[2]=0.;

  // Zero the arrays
  for(unsigned int i=0; i<nLayer; ++i) { 
    dS[i]=0.; Qp[i]=0.; Qm[i]=0.; CTp[i]=0.; 
    CTm[i]=0.; Skipp[i]=0; Skipm[i]=0;
  }
  
  // Organise timing for the Far Detector
  // Parameters for PE vs time fit residual
  //GMA all thee number need to be updated for INO
  MinUncertainty=0.56;
  ErrorParam[0]=0.56; ErrorParam[1]=0.50; ErrorParam[2]=-0.34;

  //  cout <<" fInoTrackFitAlg : TimingFit" << endl;
  // Loop over all planes
  for(int i=MinPlane; i<=MaxPlane; ++i) {
    
    if(InitTrkClustData[i].size()>0) {
      dS[i]=Sign*(dSOffset-fTrackCand->GetdS(i));
      
      CTp[i]=C*fTrackCand->GetT(i); //,StripEnd::kPositive);
      //      CTm[i]=C*fTrackCand->GetT(i,StripEnd::kNegative);
      
      if(CTp[i]>MinCT && CTp[i]<StripListTime) {StripListTime=CTp[i];}
      //        if(CTm[i]>MinCT && CTm[i]<StripListTime) {StripListTime=CTm[i];}
      
      for(unsigned int j=0; j<InitTrkClustData[i].size(); ++j) {
	Qp[i]+=InitTrkClustData[i][j].csh->GetPulse(); //StripEnd::kPositive);
	//          Qm[i]+=InitTrkClustData[i][j].csh->GetPulse(StripEnd::kNegative);
      }
    }
  }
  
  // Subtract StripList time
  if(StripListTime<8.e30) {
    for(int i=MinPlane; i<=MaxPlane; ++i) {
      if(InitTrkClustData[i].size()>0) {
	CTp[i]-=StripListTime;
	//          CTm[i]-=StripListTime;
      }
    }
  } else {
    StripListTime=0.;
  }
  //  cout <<" eInoTrackFitAlg : TimingFit" << endl;  
  // Carry out a simple straight line fit for T vs dS
  // Sqt: sum of charge*time, Sqss: sum of charge*dS*dS, etc.
  double Sqs=0; double Sqt=0; double Sqss=0; double Sqst=0; double Sqtt=0; double Sq=0;
  double TimeSlope=-999; double TimeOffset=-999; double RMS=-999;
  double CTCut = 0.; bool CalculateChi2=true;
  
  // On first iteration, carry out simple fit. Remove outlying points on subsequent passes.
  for(int itr=0; itr<3; ++itr) {
 
    for(int i=MinPlane; i<=MaxPlane; ++i) { 
      // Only consider planes where we found our final strips
      if(InitTrkClustData[i].size()>0) { 
        
        // For positive strip ends      
        s=dS[i]; q=Qp[i]; t=CTp[i];
        
        if(q>0. && t>MinCT && Skipp[i]==0) {
          if(itr==0) {Sq+=q; Sqs+=q*s; Sqt+=q*t; Sqss+=q*s*s; Sqst+=q*s*t; Sqtt+=q*t*t; n++;}
          
          else if(fabs(t-TimeOffset-(s*TimeSlope)) > CTCut) {
            Sqs-=q*s; Sqt-=q*t; Sqss-=q*s*s; Sqst-=q*s*t; Sqtt-=q*t*t; Sq-=q; n--; Skipp[i]=1;
          }
        }
	
//		
//        // For negative strip ends
//        q=Qm[i]; t=CTm[i];
//        if(q>0. && t>MinCT && Skipm[i]==0) {
//          if(itr==0) {Sq+=q; Sqs+=q*s; Sqt+=q*t; Sqss+=q*s*s; Sqst+=q*s*t; Sqtt+=q*t*t; n++;}
//          
//          else if(fabs(t-TimeOffset-(s*TimeSlope)) > CTCut) {
//            Sqs-=q*s; Sqt-=q*t; Sqss-=q*s*s; Sqst-=q*s*t; Sqtt-=q*t*t; Sq-=q; n--; Skipm[i]=1;
//          }
//        }
//	
      }
    }
    
    // Calculate parameters
    if( (Sq*Sqss-Sqs*Sqs)!=0. && Sq!=0. ) {
      TimeSlope  = (Sq*Sqst-Sqs*Sqt)/(Sq*Sqss-Sqs*Sqs);
      TimeOffset = (Sqt*Sqss-Sqs*Sqst)/(Sq*Sqss-Sqs*Sqs);
      if( ((Sqtt/Sq)-((Sqt/Sq)*(Sqt/Sq)))>0. ) {
        RMS = pow((Sqtt/Sq)-((Sqt/Sq)*(Sqt/Sq)),0.5);
        CTCut = 3.+RMS;
      }
      else {CTCut = 3.5;}
    }
    else  {CalculateChi2=false; break;}
  }

  //    cout <<" dInoTrackFitAlg : TimingFit" << endl;

  // Set timing properties for the fitted track
  if(n!=0 && CalculateChi2==true) {

    // Offset, slope and vtx/end times
    fTrackCand->SetTimeOffset((TimeOffset+StripListTime)/C);
    fTrackCand->SetTimeSlope(TimeSlope/C);

    if(ZIncreasesWithTime==true) {
      fTrackCand->SetVtxT((TimeOffset+StripListTime)/C); 
      fTrackCand->SetEndT((TimeOffset+StripListTime)/C+(dS[MaxPlane]*TimeSlope/C));
    }
    else {
      fTrackCand->SetEndT((TimeOffset+StripListTime)/C); 
      fTrackCand->SetVtxT((TimeOffset+StripListTime)/C+(dS[MaxPlane]*TimeSlope/C));
    }
    //    cout <<" cInoTrackFitAlg : TimingFit" << endl;
    // Chi2
    double Uncertainty; double Residual2; double Chi2=0; 
 
    for(int i=MinPlane; i<=MaxPlane; ++i) { 
      
      if(InitTrkClustData[i].size()>0) { 
        // For positive strip ends
        s=dS[i]; q=Qp[i]; t=CTp[i];
        if(q>0. && t>MinCT && Skipp[i]==0) {
          Residual2=pow(t-TimeOffset-(s*TimeSlope),2);
	  
          // From a rough parameterisation of uncertainty (in CT) vs number of PEs
          if (q<20) {Uncertainty = ErrorParam[0]+exp(ErrorParam[1]+ErrorParam[2]*q);}
          else {Uncertainty=MinUncertainty;}
	  
          if(Uncertainty!=0.) {Chi2+=Residual2/pow(Uncertainty,2);}
        }
//	
//        // For negative strip ends
//        q=Qm[i]; t=CTm[i];
//        if(q>0. && t>MinCT && Skipm[i]==0) {
//          Residual2=pow(t-TimeOffset-(s*TimeSlope),2);
//
//          // From a rough parameterisation of uncertainty (in CT) vs number of PEs
//          if (q<20) {Uncertainty = ErrorParam[0]+exp(ErrorParam[1]+ErrorParam[2]*q);}
//          else {Uncertainty=MinUncertainty;}
//          
//          if(Uncertainty!=0.) {Chi2+=Residual2/pow(Uncertainty,2);}
//        }
//
      }
      
    }
    // Set these properties
    fTrackCand->SetTimeFitChi2(Chi2);
    fTrackCand->SetNTimeFitDigit(n);
  }
  
  // Now carry out fits with gradients constrained to be +/- c
  double CTIntercept[2]; double Csigma[2]; double Ctrunc[2]; 
  double ChiSqPositive=-999; double ChiSqNegative=-999; 
  int ChiSqNdfPos=-999; int ChiSqNdfNeg=-999;  
  double Swtt[2]; double Swt[2]; double Sw[2]; int npts[2]={0,0};
  //  cout <<" bInoTrackFitAlg : TimingFit" << endl;
  if(Sq!=0.) {
    CTIntercept[0]=Sqt/Sq; Csigma[0]=-99999.9; Ctrunc[0]=-99999.9;
    CTIntercept[1]=Sqt/Sq; Csigma[1]=-99999.9; Ctrunc[1]=-99999.9;
    
    for(int itr=0; itr<2; ++itr)
      {
        Swtt[0]=0.; Swt[0]=0.; Sw[0]=0.; npts[0]=0;
        Swtt[1]=0.; Swt[1]=0.; Sw[1]=0.; npts[1]=0;
        
        for(unsigned int i=0; i<nLayer; ++i)
          {
            // For positive strip ends
            if(Qp[i]>0. && CTp[i]>MinCT) {
              q=Qp[i]; 
              
              t=CTp[i]-dS[i]+CTIntercept[0];
              if(Ctrunc[0]<0. || fabs(t)<Ctrunc[0]) {Swtt[0]+=q*t*t; Swt[0]+=q*t; Sw[0]+=q; ++npts[0];}
              
              t=CTp[i]+dS[i]+CTIntercept[1];
              if(Ctrunc[1]<0. || fabs(t)<Ctrunc[1]) {Swtt[1]+=q*t*t; Swt[1]+=q*t; Sw[1]+=q; ++npts[1];}
            }
//
//            // For negative strip ends
//            if(Qm[i]>0. && CTm[i]>MinCT) {
//              q=Qm[i]; 
//              
//              t=CTm[i]-dS[i]+CTIntercept[0];
//              if(Ctrunc[0]<0. || fabs(t)<Ctrunc[0]) {Swtt[0]+=q*t*t; Swt[0]+=q*t; Sw[0]+=q; ++npts[0];}
//              
//              t=CTm[i]+dS[i]+CTIntercept[1];
//              if(Ctrunc[1]<0. || fabs(t)<Ctrunc[1]) {Swtt[1]+=q*t*t; Swt[1]+=q*t; Sw[1]+=q; ++npts[1];}
//            }
//	   
          }
        
        // Results for fit with gradient +C
        if(npts[0]>1 && Sw[0]!=0.) {
          CTIntercept[0]=CTIntercept[0]-Swt[0]/Sw[0]; Csigma[0]=0.; 
          if((Swtt[0]/Sw[0])-(Swt[0]/Sw[0])*(Swt[0]/Sw[0])>0.) {Csigma[0]=pow((Swtt[0]/Sw[0])-(Swt[0]/Sw[0])*(Swt[0]/Sw[0]),0.5);}
          ChiSqPositive=Csigma[0]; ChiSqNdfPos=npts[0]-1; 
          Ctrunc[0]=Csigma[0]+3.;
        }
        
        // Results for fit with gradient -C
        if(npts[1]>1 && Sw[1]!=0.) {
          CTIntercept[1]=CTIntercept[1]-Swt[1]/Sw[1]; Csigma[1]=0.;
          if((Swtt[1]/Sw[1])-(Swt[1]/Sw[1])*(Swt[1]/Sw[1])>0.) {Csigma[1]=pow((Swtt[1]/Sw[1])-(Swt[1]/Sw[1])*(Swt[1]/Sw[1]),0.5);}
          ChiSqNegative=Csigma[1]; ChiSqNdfNeg=npts[1]-1; 
          Ctrunc[1]=Csigma[1]+3.;
        }
        
      }
  }
  // Set these properties
  fTrackCand->SetTimeForwardFitRMS(ChiSqPositive);
  fTrackCand->SetTimeForwardFitNDOF(ChiSqNdfPos);
  fTrackCand->SetTimeBackwardFitRMS(ChiSqNegative);
  fTrackCand->SetTimeBackwardFitNDOF(ChiSqNdfNeg);

  //cout <<" aInoTrackFitAlg : TimingFit" << ChiSqPositive <<" "<<ChiSqNdfPos<<" "<<ChiSqNegative<<" "<<ChiSqNdfNeg<< endl;
  if(debug_fit)  cout <<"----------------------InoTrackFitAlg : TimingFit----------------------" << endl;
}


void InoTrackFitAlg::SetRangeAnddS( ) //CandFitTrackCamHandle& cth)
{

  // Set range and dS as calculated by the swimmer
//  cout <<" InoTrackFitAlg : SetRangeAnddS from swimmer values " << endl;

  //  int ZDir; 
  //  int VtxPlane=5000; int EndPlane=-20; int Increment=1;
  int VtxPlane=5000; int EndPlane=0; int Increment=1;
  double dS; double dRange; double dP;

  // Start at the end of the track and calculate the required additions to range

  // find ending Z position (defined as Z position where muon has 100 MeV of residual energy.  This corresponds to 1/2 inch of Fe.

  // NOTE: Average dP for 1" iron is 95 MeV.
  /* 
  if(ZIncreasesWithTime==true) {ZDir=1; EndPlane=MaxPlane; VtxPlane=MinPlane; Increment=-1;}
  else {ZDir=-1; EndPlane=MinPlane; VtxPlane=MaxPlane; Increment=1;}

  PlexPlaneId plnid(detector,EndPlane,false);
  PlexPlaneId endplnid(detector,EndPlane,false);

  double Zscint = ZPosLayer[EndPlane]; //SlcClustData[EndPlane][0].csh->GetZPos();
  double Znextscint = Zscint;

  UgliScintPlnHandle scintpln = ugh.GetScintPlnHandle(plnid); 
  double Zend = Zscint + double(ZDir)*scintpln.GetHalfThickness();

  PlexPlaneId nextscint =  endplnid.GetAdjoinScint(ZDir);
  UgliScintPlnHandle nextscintpln = ugh.GetScintPlnHandle(nextscint);
  if(nextscintpln.IsValid() && nextscint.GetPlaneView()!=PlaneView::kUnknown){
    Znextscint = nextscintpln.GetZ0();
  }
  else{
    nextscint = endplnid;
  }

  plnid = plnid.GetAdjoinSteel(ZDir);
  if(plnid.IsValid()){
    UgliSteelPlnHandle steelpln = ugh.GetSteelPlnHandle(plnid);
    Zend = steelpln.GetZ0() - double(ZDir)*steelpln.GetHalfThickness();
  }

  // add two planes of steel for the ND spectrometer
  if(detector==Detector::kNear && EndPlane>=121) {
    for(int i=0;i<2;i++){
      if(plnid.GetAdjoinSteel(ZDir).IsValid()){
        PlexPlaneId plnid_after = plnid.GetAdjoinSteel(ZDir);
        if(plnid_after.IsValid()) {
          plnid = plnid_after;
          UgliSteelPlnHandle steelpln = ugh.GetSteelPlnHandle(plnid);
          Zend = steelpln.GetZ0() - double(ZDir)*steelpln.GetHalfThickness();
        }
      }
    }
  }
  */

  //GMA This is to arbitray value just to have compilation
  //    Need to put the RPC/IRON postion properly
  double Zend = LayerThickness*(EndPlane-1);

  //  double angle = 0;
  //  double ct = cos(angle);
  //  double st = sin(angle);
  // Determine whether track stops in coil
  //  float u_end = FilteredData[EndPlane][0].x_k0;
  //  float v_end = FilteredData[EndPlane][0].x_k1;
  //  float du_end = FilteredData[EndPlane][0].x_k2;
  //  float dv_end = FilteredData[EndPlane][0].x_k3;
  //  float delz = 0.0852; // Znextscint-Zscint;
  //  float u_extrap = u_end +delz*du_end;
  //  float v_extrap = v_end +delz*dv_end;
  //  float x_extrap = (ct*u_extrap-st*v_extrap);
  //  float y_extrap = (st*u_extrap+ct*v_extrap);

  /*
  PlaneCoverage::PlaneCoverage_t kPC = PlaneCoverage::kComplete;
  if(detector==Detector::kNear) kPC=PlaneCoverage::kNearFull;

  bool isInOutline = fPL.IsInside(x_extrap,y_extrap,nextscint.GetPlaneView(),kPC,false);
  bool isInCoil = isInOutline && !fPL.IsInside(x_extrap,y_extrap,nextscint.GetPlaneView(),kPC,true);
  */

  double S = 0; double Range = 0; double Prange = 0.0;
  double StateVector[6]; double Output[6];
  double chargesign = -1;
  bool GoForward = true; bool done=true; bool swimOK=true;

  // if in coil find midpoint and swim towards last clust from there
    /*
  if(isInCoil){
    float zCoil = Znextscint;
    float u_extrapC = u_extrap;
    float v_extrapC = v_extrap;
    float x_extrapC = x_extrap;
    float y_extrapC = y_extrap;
    while(isInCoil){
      zCoil -= 1.0*Munits::cm*ZDir;
      float delzC = zCoil - Zscint;
      u_extrapC = u_end +delzC*du_end;
      v_extrapC = v_end +delzC*dv_end;
      x_extrapC = (ct*u_extrapC-st*v_extrapC);
      y_extrapC = (st*u_extrapC+ct*v_extrapC);
      //      isInCoil = !fPL.IsInside(x_extrapC,y_extrapC,nextscint.GetPlaneView(),kPC,true);        
    }
    float zMinCoil = zCoil;
    if(zMinCoil<Zscint && ZDir==1)  zMinCoil=Zscint;
    if(zMinCoil>Zscint && ZDir==-1) zMinCoil=Zscint;

    zCoil = Znextscint;

    isInCoil = true;
    while(isInCoil){
      zCoil += 1.0*Munits::cm*ZDir;
      float delzC = zCoil - Zscint;
      u_extrapC = u_end +delzC*du_end;
      v_extrapC = v_end +delzC*dv_end;
      x_extrapC = (ct*u_extrapC-st*v_extrapC);
      y_extrapC = (st*u_extrapC+ct*v_extrapC);
      isInCoil = !fPL.IsInside(x_extrapC,y_extrapC,nextscint.GetPlaneView(),kPC,true);        
    }


    float zMaxCoil = zCoil;
    float zmin; float zmax;
    ugh.GetZExtent(zmin,zmax);
    if(zMaxCoil>zmax && ZDir==1)  zMaxCoil=zmax;
    if(zMaxCoil<zmin && ZDir==-1) zMaxCoil=zmin;

    // now swim from mid-coil back to endplane
    float zMidCoil = 0.5*(zMinCoil + zMaxCoil);
    float delzC  = zMidCoil -Zscint;
    u_extrapC = u_end +delzC*du_end;
    v_extrapC = v_end +delzC*dv_end;
    x_extrapC = 0.707*(u_extrapC-v_extrapC);
    y_extrapC = 0.707*(u_extrapC+v_extrapC);
    
    StateVector[0] = u_extrapC; Output[0]=StateVector[0];
    StateVector[1] = v_extrapC; Output[1]=StateVector[1];
    StateVector[2] = FilteredData[EndPlane][0].x_k2; Output[2]=StateVector[2];
    StateVector[3] = FilteredData[EndPlane][0].x_k3; Output[3]=StateVector[3];
    chargesign = -1;
    if(FilteredData[EndPlane][0].x_k4!=0.) {chargesign =  FilteredData[EndPlane][0].x_k4/fabs(FilteredData[EndPlane][0].x_k4);}
    
    GoForward = !ZIncreasesWithTime;
    StateVector[4] = 10.*chargesign; Output[4]=StateVector[4]; 
    
    StateVector[5] = FilteredData[EndPlane][0].x_k5; Output[5]=StateVector[5];

    double dsdz = pow((1. + pow(StateVector[2],2) + pow(StateVector[3],2)),0.5);    
    // set fallback to nominal energy loss in case coil swim fails
    Prange = 0.095*dsdz;
    if(detector==Detector::kNear && EndPlane>121) Prange = 0.2*dsdz;  
    Prange += 0.5*dsdz*0.1*fabs(zMaxCoil-zMinCoil)*2.357*1.97;

    swimOK = Swim(StateVector, Output, zMidCoil, EndPlane , GoForward, &dS, &dRange, &dP);
   
    if(swimOK ){
      S = dS; Range = dRange; Prange = fabs(dP);
      fTrackCand->SetdS(EndPlane,S); 
      fTrackCand->SetRange(EndPlane,Range);
    }
    if(!swimOK) {Output[4] = chargesign/Prange;}
  } else
    */

  {
    // // normal case - track does not end in coil
    //    if((Zend<Zscint && ZDir==1) || (Zend>Zscint && ZDir==-1)) {
    //      cout <<" InoTrackFitAlg :  Zend on wrong side of last scint! " << endl;
    //      Zend=Zscint;
    //    }
    
    // now swim to Zend
    StateVector[0]=FilteredData[EndPlane][0].x_k0; Output[0]=StateVector[0];
    StateVector[1]=FilteredData[EndPlane][0].x_k1; Output[1]=StateVector[1];
    StateVector[2]=FilteredData[EndPlane][0].x_k2; Output[2]=StateVector[2];
    StateVector[3]=FilteredData[EndPlane][0].x_k3; Output[3]=StateVector[3];
    StateVector[4]=FilteredData[EndPlane][0].x_k4; Output[4]=StateVector[4];
    StateVector[5]=FilteredData[EndPlane][0].x_k5; Output[5]=StateVector[5];

    chargesign = -1;
    if(StateVector[4]!=0.) {chargesign = StateVector[4]/fabs(StateVector[4]);}
    
    GoForward = ZIncreasesWithTime;
    done = Swim(StateVector, Output, EndPlane, Zend, GoForward,  &dS, &dRange, &dP);
    
    GoForward = !ZIncreasesWithTime;
    double dsdz = pow((1. + pow(StateVector[2],2) + pow(StateVector[3],2)),0.5); 
    S = 0; Range = 10.0*dsdz; Prange = 0.095*dsdz;
    swimOK = false;
    if(done){
      for(int j=0;j<6;j++) {StateVector[j]=Output[j];}
      
      // now swim from Zend to EndPlane
      StateVector[4] = chargesign * 10.52;  // start @ P = 100 MeV (Eloss in 1/2 " Iron)
      swimOK = Swim(StateVector, Output, Zend, EndPlane , GoForward, &dS, &dRange, &dP);
      if(swimOK){
        S += dS; Range += dRange; Prange += fabs(dP);
        fTrackCand->SetdS(EndPlane,S); 
        fTrackCand->SetRange(EndPlane,Range);
      }
    }
    if(!swimOK) {Output[4] = chargesign/Prange;}
  }
  
  int thisplane = EndPlane;
  // now swim back to vertex
  bool firstplane=true;
  for(int i=EndPlane+Increment; Increment*i<=Increment*VtxPlane; i+=Increment) {
    if(FilteredData[i].size()>0) {
      double delU = FilteredData[i][0].x_k0 - StateVector[0] ;
      double delV = FilteredData[i][0].x_k1 - StateVector[1] ;
      double dSperPlane=0.;
      if(thisplane!=i) {dSperPlane = pow(delU*delU + delV*delV,0.5)/double(abs(thisplane-i));}
      
      // only update state vector if change in U/V is reasonable.
      if(dSperPlane < 1.5) { // *TMinuit::m) {
        StateVector[0]=FilteredData[i][0].x_k0;
        StateVector[1]=FilteredData[i][0].x_k1;
        StateVector[2]=FilteredData[i][0].x_k2;
        StateVector[3]=FilteredData[i][0].x_k3;
	
        chargesign=-1;
        if(FilteredData[i][0].x_k4!=0.) {chargesign = FilteredData[i][0].x_k4/fabs(FilteredData[i][0].x_k4);}
	
	StateVector[5]=FilteredData[i][0].x_k5;
      }
      
      StateVector[4] = chargesign * fabs(Output[4]);
      done = Swim(StateVector, Output, thisplane, i , GoForward, &dS, &dRange, &dP);
      if(done){
        S+=dS; Range+=dRange; Prange+=fabs(dP);
        fTrackCand->SetdS(i,S); fTrackCand->SetRange(i,Range);
        firstplane=false;
      }
      else {
        cout <<" InoTrackFitAlg :  swim fail " << endl;
      }
      thisplane=i;
    }
  }

  /*
  //  PlexPlaneId vtxplnid(detector,VtxPlane,false);
  //  PlexPlaneId plnid_before = vtxplnid.GetAdjoinSteel(-ZDir);
  
  if(plnid_before.IsValid()) {
    plnid = plnid_before;
    UgliSteelPlnHandle steelpln = ugh.GetSteelPlnHandle(plnid);
    double Zstart = steelpln.GetZ0();
    StateVector[0]=FilteredData[VtxPlane][0].x_k0;
    StateVector[1]=FilteredData[VtxPlane][0].x_k1;
    StateVector[2]=FilteredData[VtxPlane][0].x_k2;
    StateVector[3]=FilteredData[VtxPlane][0].x_k3;
    StateVector[4]=Output[4];
    StateVector[5]=FilteredData[VtxPlane][0].x_k5;

    Swim(StateVector, Output, VtxPlane, Zstart, GoForward, &dS,&dRange,&dP);
    S+=dS; Range+=dRange; Prange+=fabs(dP);

    track->SetRange(VtxPlane,Range);
    track->SetdS(VtxPlane,S);
  }
  */
  
  // if Prange < 21 GeV, use this value.  Otherwise, use finder track energy, which is somewhat less prone to gross errors.
  
  // apply fudge factor for nominal steel thickness in ND geometry.
  
  //  track->SetMomentumRange(Prange*ecorr);
  //  CandTrackHandle* findertrack = track->GetFinderTrack();
  //  if(((detector==Detector::kFar && Prange>21.) || (detector==Detector::kNear && Prange>12.)) && findertrack) {track->SetMomentumRange(findertrack->GetMomentum());}
  
}

/*
void InoTrackFitAlg::SetPropertiesFromFinderTrack(CandFitTrackCamHandle &cthx)
{
  // This method is only called if the fit fails. We set properties from finder track.
  // This clearly does not include fitted properties such as q/p or QPVtxError.
 cout <<" InoTrackFitAlg : SetPropertiesFromFinderTrack" << endl;
  cthx.SetDirCosU(track->GetDirCosU());
  cthx.SetDirCosV(track->GetDirCosV());
  cthx.SetDirCosZ(track->GetDirCosZ());
  cthx.SetVtxU(track->GetVtxU());
  cthx.SetVtxV(track->GetVtxV());
  cthx.SetVtxZ(track->GetVtxZ());
  cthx.SetVtxT(track->GetVtxT());
  cthx.SetVtxPlane(track->GetVtxPlane());

  cthx.SetEndDirCosU(track->GetEndDirCosU());
  cthx.SetEndDirCosV(track->GetEndDirCosV());
  cthx.SetEndDirCosZ(track->GetEndDirCosZ());
  cthx.SetEndU(track->GetEndU());
  cthx.SetEndV(track->GetEndV());
  cthx.SetEndZ(track->GetEndZ());
  cthx.SetEndT(track->GetEndT());
  cthx.SetEndPlane(track->GetEndPlane());

  cthx.SetMomentumRange(track->GetMomentum());
  cthx.SetMomentum(track->GetMomentum());

  cthx.SetTimeSlope(track->GetTimeSlope());
  cthx.SetTimeOffset(track->GetTimeOffset());
  cthx.SetTimeFitChi2(track->GetTimeFitChi2());
  cthx.SetNTimeFitDigit(track->GetNTimeFitDigit());
  cthx.SetTimeForwardFitRMS(track->GetTimeForwardFitRMS());
  cthx.SetTimeForwardFitNDOF(track->GetTimeForwardFitNDOF());
  cthx.SetTimeBackwardFitRMS(track->GetTimeBackwardFitRMS());
  cthx.SetTimeBackwardFitNDOF(track->GetTimeBackwardFitNDOF());


  // Set quantities required at each plane in finder track
  int direction=1;
  if(ZIncreasesWithTime==false) {direction=-1;}

  for(int i=cthx.GetVtxPlane(); i*direction<=cthx.GetEndPlane()*direction; i+=direction){
    if(track->IsTPosValid(i)) {
      cthx.SetTrackPointError(i,track->GetTrackPointError(i));
      cthx.SetdS(i,track->GetdS(i));
      cthx.SetRange(i,track->GetRange(i));
      cthx.SetU(i,track->GetU(i));
      cthx.SetV(i,track->GetV(i));
    }
  }

  CalculateTrace(); //cthx);  
  SetT(); //&cthx);

  Calibrate(); //&cthx);

}
*/

void InoTrackFitAlg::Trace(const char * /* c */) const
{;
}


void InoTrackFitAlg::SetT()
{

// we take a weighted average of clusts in the same
// plane.  the proper way to do this would be to consider clusts coming
// from the same PMT and use the earliest time.

  //cout <<" inoTrackFitAlg : starting SetT" << endl;
  for (unsigned i=0; i<fFinderTrack->ClustsInTrack.size(); i++) {
    

  }
  
}

/*
void AlgTrack::CalculateTrace( )
  
  double m[2]; double c[2]; 
  TIter stripItr(fTrackCand->GetDaughterIterator());
  CandStripHandle *firststrip = 0;
  firststrip = dynamic_cast<CandStripHandle*>(stripItr());
  if (!firststrip) {
    MSG("RecoBaseAlgTrack",Msg::kDebug) << "no strips, returning" << endl;
  }
  if (!firststrip) return; // if no strips do nothing
  //RWH: is this really necesary? can't we get the VldContext
  //     from the FTrackCandHandle itself?

  const VldContext *vldc = firststrip->GetVldContext();
 
  UgliGeomHandle ugh(*vldc); 
  PlaneOutline pl;
  float detzmin=0;
  float detzmax=999;
  ugh.GetZExtent(detzmin,detzmax,-1);

  int dir=1;
  if(fTrackCand->GetVtxZ()>fTrackCand->GetEndZ())dir=-1;

  // For vertex - loop over planes upstream from vertex until exit detector
  // trace = linear extrapolation to track entry point
  bool IsOutside = false; bool IsInside = true;
  double trace=0;
  double tracez=0;
  double z = fTrackCand->GetVtxZ();
  double u=0; double v=0; double x=0; double y=0; double r=0;
  float xedge,yedge,dist;
  int nActiveUpstream = 0;
  Detector::Detector_t detector = vldc->GetDetector();
  if(fTrackCand->GetVtxDirCosZ()!=0.) {
    m[0]=fTrackCand->GetVtxDirCosU()/fTrackCand->GetVtxDirCosZ();
    m[1]=fTrackCand->GetVtxDirCosV()/fTrackCand->GetVtxDirCosZ();
    c[0]=fTrackCand->GetVtxU()-(m[0]*fTrackCand->GetVtxZ());
    c[1]=fTrackCand->GetVtxV()-(m[1]*fTrackCand->GetVtxZ());
    PlexPlaneId plnid(detector,fTrackCand->GetVtxPlane(),false);
    while(!IsOutside){
      plnid = plnid.GetAdjoinScint(-dir);

      if(plnid.IsValid()){
        UgliScintPlnHandle scintpln = ugh.GetScintPlnHandle(plnid);
        z = scintpln.GetZ0();
        u = c[0] + m[0]*z;
        v = c[1] + m[1]*z;
        x = 0.707*(u-v);
        y = 0.707*(u+v);
        r = pow(x*x+y*y,0.5);
        pl.DistanceToEdge(x, y,
                          plnid.GetPlaneView(),
                          plnid.GetPlaneCoverage(),
                          dist, xedge, yedge);
        
       IsInside = pl.IsInside( x, y,
                               plnid.GetPlaneView(),
                               plnid.GetPlaneCoverage(),true);
       IsInside |= r<0.5;
       if(IsInside){
         nActiveUpstream++;
         trace = pow(pow(u-fTrackCand->GetVtxU(),2)+pow(v-fTrackCand->GetVtxV(),2)+pow(z-fTrackCand->GetVtxZ(),2),0.5);
         tracez = fabs(z - fTrackCand->GetVtxZ());
       }
      }
      IsOutside = !IsInside | !plnid.IsValid();
    }
    if(plnid.IsValid()) trace += dist;
  }
  fTrackCand->SetVtxTrace(trace);
  fTrackCand->SetVtxTraceZ(tracez);
  fTrackCand->SetVtxnActiveUpstream(nActiveUpstream);

  // set distance from track vertext to nearest edge on that plane
  z = fTrackCand->GetVtxZ();
  u = fTrackCand->GetVtxU();
  v = fTrackCand->GetVtxV();
  x = 0.707*(u-v);
  y = 0.707*(u+v);
  r = pow(x*x+y*y,0.5);
  PlexPlaneId plnid(detector,fTrackCand->GetVtxPlane(),false);
  dist = 0.;
  IsInside=true;
  if(plnid.IsValid()){
    pl.DistanceToEdge(x, y,
                      plnid.GetPlaneView(),
                      plnid.GetPlaneCoverage(),
                      dist, xedge, yedge);
    IsInside = pl.IsInside( x, y,
                            plnid.GetPlaneView(),
                            plnid.GetPlaneCoverage(),true);
    
  }
  if(!IsInside && r<0.5) dist=0;
  fTrackCand->SetVtxDistToEdge(dist);

  // For vertex - loop over planes upstream from vertex until exit detector
  // trace = linear extrapolation to track entry point
  IsOutside = false; IsInside=true;
  trace=0; tracez=0;
  z = fTrackCand->GetEndZ();
  u=0; v=0; x=0; y=0; r=0;
  int nActiveDownstream = 0;
  if(fTrackCand->GetEndDirCosZ()!=0.) {
    m[0]=fTrackCand->GetEndDirCosU()/fTrackCand->GetEndDirCosZ();
    m[1]=fTrackCand->GetEndDirCosV()/fTrackCand->GetEndDirCosZ();
    c[0]=fTrackCand->GetEndU()-(m[0]*fTrackCand->GetEndZ());
    c[1]=fTrackCand->GetEndV()-(m[1]*fTrackCand->GetEndZ());
   
    PlexPlaneId plnid(detector,fTrackCand->GetEndPlane(),false);

    while(!IsOutside){
      plnid = plnid.GetAdjoinScint(dir);
      if(plnid.IsValid()){
        UgliScintPlnHandle scintpln = ugh.GetScintPlnHandle(plnid);
        z = scintpln.GetZ0();
        u = c[0] + m[0]*z;
        v = c[1] + m[1]*z;
        x = 0.707*(u-v);
        y = 0.707*(u+v);
        r = pow(x*x+y*y,0.5);
        pl.DistanceToEdge(x, y,
                          plnid.GetPlaneView(),
                          plnid.GetPlaneCoverage(),
                          dist, xedge, yedge);
        
        IsInside = pl.IsInside( x, y,
                                plnid.GetPlaneView(),
                                plnid.GetPlaneCoverage(),true);
        IsInside |= r<0.5;
        if(IsInside){
          nActiveDownstream++;
          trace = pow(pow(u-fTrackCand->GetEndU(),2)+pow(v-fTrackCand->GetEndV(),2)+pow(z-fTrackCand->GetEndZ(),2),0.5);
          tracez = fabs(z - fTrackCand->GetEndZ());
        }
      }
      IsOutside = !IsInside | !plnid.IsValid();
    }
    if(plnid.IsValid()) trace +=dist;
  }
  fTrackCand->SetEndTrace(trace);
  fTrackCand->SetEndTraceZ(tracez);
  fTrackCand->SetEndnActiveDownstream(nActiveDownstream);

  // set distance from track end to nearest edge on that plane
  z = fTrackCand->GetEndZ();
  u = fTrackCand->GetEndU();
  v = fTrackCand->GetEndV();
  x = 0.707*(u-v);
  y = 0.707*(u+v);
  r = pow(x*x+y*y,0.5);
       
  PlexPlaneId plnidend(detector,fTrackCand->GetEndPlane(),false);
  dist = 0.;
  IsInside=true;
  if(plnidend.IsValid()){
    pl.DistanceToEdge(x, y,
                      plnidend.GetPlaneView(),
                      plnidend.GetPlaneCoverage(),
                      dist, xedge, yedge);
    IsInside = pl.IsInside( x, y,
                            plnidend.GetPlaneView(),
                            plnidend.GetPlaneCoverage(),true);

  }
  // set distance to edge =0 if inside coil.  This forces tracks that end in
  // coil to be !contained. 
  if(!IsInside && r<0.5) dist=0; 
  
  fTrackCand->SetEndDistToEdge(dist);
}
 -0.0192623 0.00426438 -0.0664026 0.0257764 0
/afs/cern.ch/user/m/majumder/scratch0/anal/CMSSW_2_0_0_pre1/src/Test/InclMuonPdfAnalyser/test/MyFirstFamosFile.root

*/

//GMA 070809 : modified code looks for  Plane+-1 layer, but not NewPlane. We have to modify this

bool InoTrackFitAlg::GetCombiPropagator_new(const int Plane, const int NewPlane, const bool GoForward)
{
  // Combination propagator, essentially same as SR propagator, but
  // generation of matrix reduces calls to swimmer by 80%
  cout <<" InoTrackFitAlg : GetCombiPropagator" << endl;

  for (int i=0; i<5; ++i) {
    for (int j=0; j<5; ++j) {
      F_k_minus[i][j]=0;  
    }
  }

  F_k_minus[0][0]=1; F_k_minus[1][1]=1; 
  F_k_minus[2][2]=1; F_k_minus[3][3]=1;
  
  DeltaZ=fabs(TrkClustsData[NewPlane][0].ZPos-TrkClustsData[Plane][0].ZPos);
  DeltaPlane=abs(NewPlane-Plane);

  // Swimmer section for last column
  double PState[6];  double NState[6];  double StateVector[6];
  double Increment=0.01;
  bool SwimInc=false; bool SwimDec=false;
  int nswimfail=0;
  
  if(GoForward==true) {F_k_minus[0][2]=DeltaZ; F_k_minus[1][3]=DeltaZ;}
  else if(GoForward==false) {F_k_minus[0][2]=-DeltaZ; F_k_minus[1][3]=-DeltaZ;}
  

  // Give swimmer fixed number of opportunities for successful swim
  while((SwimInc==false || SwimDec==false) && nswimfail<=10) {

    Increment=0.05*fabs(x_k_minus[4]); 
    if(Increment<.01) {Increment=.01;}
    
    //GMAA
    //    if (x_k_minus[4]==0) Increment=0.2;


    for(int j=0; j<6; ++j) {StateVector[j]=x_k_minus[j];}  

    int nextp=0;
    // Increment then swim
    StateVector[4]+=Increment;
    SwimInc=Swim_new(StateVector, NState, Plane, nextp, GoForward);
    

    if (NewPlane !=nextp) {
      cout<<"CombipropagatorF error "<< Plane<<" "<<nextp<<" "<<NewPlane<<endl;
      cout<<" "<<StateVector[0]<<" "  <<StateVector[1]<<" "<<StateVector[2]<<" "<<StateVector[3]<<" "<<StateVector[4]<<endl;
      cout<<" "<<NState[0]<<" "  <<NState[1]<<" "<<NState[2]<<" "<<NState[3]<<" "<<NState[4]<<endl;
    }
    StateVector[4]=x_k_minus[4];
    
    // Decrement then swim
    StateVector[4]-=Increment;
    SwimDec=Swim_new(StateVector, PState, Plane, nextp, GoForward);
    if (!(NewPlane !=nextp)) {
      cout<<"CombipropagatorB error "<< Plane<<" "<<nextp<<" "<<NewPlane<<endl;
      cout<<" "<<StateVector[0]<<" "  <<StateVector[1]<<" "<<StateVector[2]<<" "<<StateVector[3]<<" "<<StateVector[4]<<endl;
      cout<<" "<<PState[0]<<" "  <<PState[1]<<" "<<PState[2]<<" "<<PState[3]<<" "<<PState[4]<<endl;
    }
    

    // If swim failed, double momentum and swim again
    if(SwimInc==false || SwimDec==false) {
      // GMA reopen this      cout <<" InoTrackFitAlg : GetCombiPropagator, Swim failed - Double momentum and swim again" << endl;
      x_k_minus[4]*=0.5;
      nswimfail++; TotalNSwimFail++;
      break;                           //AAR: note this break statement
    }
    
    // Form last row of propagator matrix.  Need to transpose to get proper Kalman F_k_minus
    else {
      if(Increment!=0.) {
        for(int j=0; j<5; ++j) {
          F_k_minus[j][4] = (NState[j]-PState[j]) / (2*Increment);
        }
      }
      else {F_k_minus[4][4]=1;}
    }
    
  } // End while statement
  
  if(nswimfail>10) {cout <<" InoTrackFitAlg : GetCombiPropagator, nswimfail>10, fail track" << endl; return false;}


  // Display
 // if(!debug_fit) {
      cout<<"------------------------------------------------------------"<<endl;
    cout << "Combi F_k_minus" << endl;
    for(int i=0; i<5; ++i) {  
      for(int j=0; j<5; ++j) {
        cout << F_k_minus[i][j] << " ";
      }
      cout << endl;
    }
    cout<<"------------------------------------------------------------"<<endl;
//  }
  
  return true;
}

void InoTrackFitAlg::UpdateStateVector_new(const int Plane, const int NewPlane, double* Prediction, const bool GoForward)
{  
  // x_k = (F_k_minus * x_k_minus) + K_k(m_k - (H_k * F_k_minus * x_k_minus) )
   //cout <<" InoTrackFitAlg : UpdateStateVector" << endl;

  int nswimfail=0;


  // Calculate F_k_minus * x_k_minus, using the Swimmer
  // Also get an accurate measure of dS and Range from the Swimmer
  double StateVector[6];
  double dummyPrediction[6];
  bool GetPrediction=false;

  for(int i=0; i<6; ++i) {StateVector[i]=x_k_minus[i];}
  
  // Do the swim
  while(GetPrediction==false && nswimfail<=10) {
    
    int tmpplane = NewPlane;
    GetPrediction=Swim_new(StateVector, dummyPrediction, Plane, tmpplane, GoForward);
    
    //GMA do not allow starting momentum to be very large
    /*
    while (x_k_minus[4]==0 && fabs(Prediction[4])<1.e-2) {
      cout <<"x_4 "<<x_k_minus[4]<<" "<<StateVector[4]<<" "<<Prediction[4]<<endl;
      StateVector[4] = 10.*Prediction[4];
      if (StateVector[4] >1.e-2) StateVector[4] =1.e-2;
      if (StateVector[4] <-1.e-2) StateVector[4] =-1.e-2;
      GetPrediction=Swim_new(StateVector, dummyPrediction, Plane, NewPlane, GoForward);
    }
    */

    if(GetPrediction==false) {
      StateVector[4]*=0.5; 
      nswimfail++; TotalNSwimFail++;
      //GMA reopen this      cout <<" InoTrackFitAlg : UpdateStateVector, Prediction failed - Double momentum and swim again" << endl;
    }
  }

  if(nswimfail>10) {  // Swim shouldn't fail, as it succeeded to get the propagator
    //GMA reopen this    cout <<" InoTrackFitAlg : UpdateStateVector, nswimfail>10, fail track" << endl;
   cout<< " Passtrack 2 "<<" ievt "<<pAnalysis->ievt2<<endl;
       PassTrack=false;
  }
  
  //  cout <<" InoTrackFitAlg : UpdateStateVector, Check predicted state " << endl;
  CheckValues(Prediction, NewPlane);

  // Calculate H_k * F_k_minus * x_k_minus
  int PlaneView = TrkClustsData[NewPlane][0].PlaneView;
  for (int i=0; i<2; i++) {
    for (int j=0; j<5; j++) {
      H_k[i][j]=0;
    }
  }
  // Calculate A1_k = m_k - H_k * F_k_minus * x_k_minus
  double A1_k[2]={0,0};
  if (PlaneView%2==0) {H_k[0][0]=1; A1_k[0] = TrkClustsData[NewPlane][0].XPos;}
  if (PlaneView   >0) {H_k[1][1]=1; A1_k[1] = TrkClustsData[NewPlane][0].YPos;}
  
  //  if ((Plane >=78 && Plane <=80)||(Plane >=126 && Plane <=128)) {
  //  cout <<"a1_k "<< PlaneView<<" "<<A1_k[0]<<" "<<A1_k[1]<<endl;
  //  }

  
  //GMA 090809 Check this shift

  double tmpPrediction[5];
  for (int j=0; j<5; j++) { tmpPrediction[j] =Prediction[j];}
  tmpPrediction[0] +=shiftvector.x();
  tmpPrediction[1] +=shiftvector.y();

  for (int i=0; i<2; i++) {
    for (int j=0; j<5; j++) {
      A1_k[i] -= H_k[i][j]*Prediction[j];



      //      if ((Plane >=78 && Plane <=80)||(Plane >=126 && Plane <=128)) {
      //      cout <<"A1_k "<<i<<" "<<j<<" "<<A1_k[i]<<" "<<H_k[i][j]<<" "<<Prediction[j]<<endl;
      //      }
    }
  }

  // Calculate x_k
  x_k[5]=Prediction[5];
  for (int i=0; i<5; ++i) {
    x_k[i]=Prediction[i];
    for (int j=0; j<2; j++) {
      x_k[i]+=K_k[i][j]*A1_k[j];
      //      if ((Plane >=78 && Plane <=80)||(Plane >=126 && Plane <=128)) {
      //      cout <<"x_k "<<i <<" "<<j<<" "<<Prediction[i]<<" "<<x_k[i]<<" "<<K_k[i][j]*A1_k[j]<<endl;
      //      }
    }
  }
  
  //  cout <<" InoTrackFitAlg : UpdateStateVector, Check filtered state " << endl;
  CheckValues(x_k, NewPlane);


  // GMA Care with multiple range corrections - do not want to flip sign
  // (multiple corrections mean sign changes can occur even though absolute value stays same)
  // JAM up to 40 from 4

  double Maxqp = 4.;
  if(fabs(x_k_minus[4])==Maxqp && 
     ( (GoForward==true && ZIncreasesWithTime==true) 
       || (GoForward==false && ZIncreasesWithTime==false) ) ) 
    {
      if(!LastIteration) x_k[4] = (x_k_minus[4]>0 ? Maxqp : -Maxqp);
      //      cout << " resetting in UpdateStateVector " <<  endl;
    }

  //if on last plane in forward swim, disregard sign flip
  if(x_k_minus[4]!=0){
    if ( x_k[4]/x_k_minus[4]<0 && 
	 ( (GoForward==true && ZIncreasesWithTime==true && NewPlane >= EndofRangePlane ) 
	   || (GoForward==false && ZIncreasesWithTime==false && NewPlane <= EndofRangePlane ) )) 
      {
        x_k[4] = -x_k[4];
      }
  }

  x_k[5] = 1; //GMA 09/02/2009  for all upgrade



   //Display
    cout <<" InoTrackFitAlg : State: Filtered  " << Plane<<" "
         << NewPlane<<" "
         <<TrkClustsData[NewPlane][0].XPos<<" "
         <<TrkClustsData[NewPlane][0].YPos<<" "
         << x_k[0] << " " << x_k[1] << " " << x_k[2] << " "
         << x_k[3] << " " << x_k[4] <<" ";
    if (x_k[4] !=0) {
      cout <<1./x_k[4] <<endl;
    } else {
      cout <<x_k[4]<<endl;
    }

}


bool InoTrackFitAlg::Swim_new(double* StateVector, double* Output, const int Plane, 
			  int& NewPlane, const bool GoForward_loc, double* dS, double* Range, double* dE)
{

  if(debug_new) cout<<"Swim_new...."<<" "<<Plane<<" "<<NewPlane<<endl;
  double loc_layerthickness = 0.096; // in meter
  
  if( (ZIncreasesWithTime==true && GoForward_loc==false) || 
      (GoForward_loc==false && ZIncreasesWithTime==false) ){
    if(Plane>0) {
      loc_layerthickness = ZPosLayer[Plane] - ZPosLayer[Plane-1];
    }
  } else if ( (ZIncreasesWithTime==false && GoForward_loc==true) ||
	      (GoForward_loc==true && ZIncreasesWithTime==false) ) {
    if(Plane<nLayer-1) {
      loc_layerthickness = ZPosLayer[Plane+1] - ZPosLayer[Plane];
    }
  }
  
  // Initialisations
  // customize for bfield scaling.
  //  BField * bf = new BField(*vldc,-1,0);
  //  SwimSwimmer* myswimmer = new SwimSwimmer(fabs(LayerThickness*(Plane-NewPlane)), 0.5*LayerThickness); //*vldc,bf);
  SwimSwimmer* myswimmer = new SwimSwimmer(Plane, loc_layerthickness, 0.5*LayerThickness); //*vldc,bf);
  int cplane = (ZIncreasesWithTime)? Plane+1 : Plane;
  myswimmer->SetBPlane(cplane);
  //  if(UseGeoSwimmer) GeoSwimmer::Instance()->Initialize(*vldc);
 
  //  double invSqrt2 = pow(1./2.,0.5);
  double charge = 0.;
  bool done = false;

  if(fabs(StateVector[4])>1.e-10) {
    double modp = fabs(1./StateVector[4]);
    if(debug_new) cout<<"In Swim_new First Option.. "<<endl;
    // Fix, to account for fact the cosmic muons could move in direction of negative z
    if(ZIncreasesWithTime==false) {modp=-modp;}

    double dsdz = pow((1.+pow(StateVector[2],2)+pow(StateVector[3],2)),0.5);

    double angle = 0;
    double ct=cos(angle);
    double st=sin(angle);

    //    double dxdz = invSqrt2*(StateVector[2]-StateVector[3]);
    //    double dydz = invSqrt2*(StateVector[2]+StateVector[3]);
    double dxdz = ct*StateVector[2]-st*StateVector[3];
    double dydz = st*StateVector[2]+ct*StateVector[3];

    // Set up current muon details
    if(StateVector[4]>0.) charge = 1.;
    else if(StateVector[4]<0.) charge = -1.;

    //    cout <<"charge2 === "<<endl;
    //    cout <<"ct "<< ct<<" "<<st<<" "<<StateVector[0]<<" "<<StateVector[1]<<" "<<ZPosLayer[Plane]<<endl; //SlcClustData[Plane][0].csh->GetZPos()<<endl;
    //    cout <<"modp "<< modp<<" "<< dxdz<<" "<<dsdz<<" "<<dydz<<endl;
    TVector3 position(ct*StateVector[0]-st*StateVector[1],
                      st*StateVector[0]+ct*StateVector[1],
                      ZPosLayer[Plane]); //SlcClustData[Plane][0].csh->GetZPos());

    TVector3 momentum(modp*(dxdz/dsdz),
                      modp*(dydz/dsdz),
                      modp/dsdz);

    //    TVector3 bfield = bf->GetBField(position);
    //TVector3 bfield(1.,1.,0.); //GMA-magnetic field   //AAR : change the field here   //AAR:commented out
    //    TVector3 bfield(1.5,0.,0.); //GMA-magnetic field
    //bave += TMath::Sqrt(bfield[0]*bfield[0]+bfield[1]*bfield[1]+bfield[2]*bfield[2]); //AAR : commented out
    //    bave += pow(bfield[0]*bfield[0]+bfield[1]*bfield[1]+bfield[2]*bfield[2],0.5);
    //nbfield++;  //AAR : commented out

    SwimParticle muon(position,momentum);

    muon.SetCharge(charge);
    //    cout <<"charge === "<<charge<<" "<<momentum.X()<<" "<<momentum.Y()<<" "<<momentum.Z()<<" "<<position.X()<<" "<<position.Y()<<" "<<position.Z()<<" "<<dxdz<<" "<<dydz<<" "<<dsdz<<" st "<<StateVector[0]<<" "<<StateVector[1]<<" "<<StateVector[2]<<" "<<StateVector[3]<<" "<<StateVector[4]<<" "<<muon.GetMomentum().Z()<<endl;

    int nextplane=0;
    //GMA    SwimZCondition zc(ZPosLayer[NewPlane]); //SlcClustData[NewPlane][0].csh->GetZPos());
    // Do the swim, accounting for direction of motion w.r.t time too

    if(debug_new) cout<<"Before Swim_new::SwimSwimmer::Swim(... "<<" "<<done<<" "<<nextplane<<" "<<muon.GetCharge()*muon.GetMomentumModulus()<<endl;

    if( (GoForward_loc==true && ZIncreasesWithTime==true)  || (GoForward_loc==false && ZIncreasesWithTime==false) ) {
      if(UseGeoSwimmer) {
	//        done = GeoSwimmer::Instance()->SwimForward(muon,ZPosLayer[NewPlane]); //SlcClustData[NewPlane][0].csh->GetZPos());
      } else {
        done = myswimmer->SwimForward(muon, nextplane,t_bave);
      }
    }
    else if( (GoForward_loc==true && ZIncreasesWithTime==false)  || (GoForward_loc==false && ZIncreasesWithTime==true) ) {
      if(UseGeoSwimmer) {
	//        done = GeoSwimmer::Instance()->SwimBackward(muon,ZPosLayer[NewPlane]); //SlcClustData[NewPlane][0].csh->GetZPos());
      } else {
        done = myswimmer->SwimBackward(muon, nextplane,t_bave);
      }
    }

    if(debug_new) cout<<"After Swim_new::SwimSwimmer::Swim(... "<<" "<<done<<" "<<nextplane<<" "<<muon.GetCharge()*muon.GetMomentumModulus()<<endl;
    
    bave += t_bave;

    nbfield++;  //AAR : commented out
    shiftvector = myswimmer->getCrossingShift();
    if(done==true) {
      double angle = 0;
      double ct=cos(angle);
      double st=sin(angle);
      if(muon.GetDirection().Z()!=0. && muon.GetMomentumModulus()!=0.) {
        Output[0]=(st*muon.GetPosition().Y()+ct*muon.GetPosition().X());
        Output[1]=(ct*muon.GetPosition().Y()-st*muon.GetPosition().X());
        Output[2]=(st*(muon.GetDirection().Y()/muon.GetDirection().Z())+ct*(muon.GetDirection().X()/muon.GetDirection().Z()));
        Output[3]=(ct*(muon.GetDirection().Y()/muon.GetDirection().Z())-st*(muon.GetDirection().X()/muon.GetDirection().Z()));
        Output[4]=muon.GetCharge()/muon.GetMomentumModulus();

	Output[5]= StateVector[5];
        // Get range and dS from the Swimmer
        if(dS) {*dS=muon.GetS();} 
	if(Range) {*Range=muon.GetRange();} 
	if(dE){*dE=muon.GetMomentumModulus()-momentum.Mag();} 
	
	NewPlane = nextplane;
	//GMA put this more elegantly
	fTrackCand->fdS[NewPlane] =muon.GetS();
	fTrackCand->fRange[NewPlane] =muon.GetRange();

      } else {
	done=false;
      }
    } else {
      if(debug_new) cout<<"done.. not true" <<endl;
    }
  } else {

    // If infinite momentum, use straight line extrapolation
    if(debug_new) cout<<"In Swim_new Second Option.. "<<endl;

    
    //    double delz = LayerThickness;
    //    cout <<"delz "<< delz<<endl;
    //    if (SlcClustData[NewPlane].size()>0 && SlcClustData[Plane].size()>0) {
    //      delz = ZPosLayer[NewPlane]-ZPosLayer[Plane]; //(SlcClustData[NewPlane][0].csh->GetZPos()-SlcClustData[Plane][0].csh->GetZPos());
    //    }

    //    cout <<"delz "<< delz<<endl;

    double costhe = 1./pow(1. + Output[2]*Output[2] + Output[3]*Output[3], 0.5);
    if(ZIncreasesWithTime==false) {costhe = -costhe;}
    if (costhe >0) {NewPlane = Plane+1; } else {NewPlane = Plane-1;}

    double delz = (NewPlane-Plane)*LayerThickness;
    
    Output[0]=StateVector[0] + StateVector[2]*delz;
    Output[1]=StateVector[1] + StateVector[3]*delz;
    Output[2]=StateVector[2];
    Output[3]=StateVector[3];
    Output[4]=StateVector[4];
    Output[5]=StateVector[5];




    done=true;
  }
  //if(StateVector[4]/Output[4] < 0){
  if(debug_new) {
    cout <<" StateVector "<< StateVector[0]<<" "<<StateVector[1]<<" "<<StateVector[2]<<" "<<StateVector[3]<<" "<<StateVector[4]<<" "<<1/StateVector[4]<<endl;
    cout <<" Output "<< Output[0]<<" "<<Output[1]<<" "<<Output[2]<<" "<<Output[3]<<" "<<Output[4]<<" "<<1/Output[4]<<" "<<done<<endl;
  }
  //}

  
  delete myswimmer;
  //  delete bf;
  return done;
}



void InoTrackFitAlg::GoForwards_new(bool first) { 
  // Carry out the Kalman fit along the track in the direction of increasing z
  if(TrkFitterDebug>0) {
    cout<< endl;
    cout << "InoTrackFitAlg : GoForwards_new, carry out fit in positive z direction --- " <<first << endl;
  }

  // JAM in 2nd iteration, stop when end of range is reached.

  double StateVector[6]; double Prediction[6];

  Int_t StartPlane = MinPlane;// Int_t EndPlane=MaxPlane;
  if(!ZIncreasesWithTime){
    StartPlane = EndofRangePlane;
  }
  else EndofRangePlane = MaxPlane;
  double x_k_old[6];

  unsigned ntrk = fTrackCand->GetEntries();
  int iend = (first) ? 0 : 1;
  //a  cout <<"ntrksize "<< ntrk<<endl;
  for (unsigned ijk=0; ijk<ntrk; ijk++) {
    //    int i=fTrackCand->ClustsInTrack[ijk]->GetZPlane();
    //a    cout<<"ijkfd_new  "<< ijk<<" "<< i<<" "<< TrkClustsData[i].size()<<" "<<fTrackCand->ClustsInTrack[ijk]->GetZPlane()<<" "<<fTrackCand->ClustsInTrack[ijk]->GetXPos()<<" "<< fTrackCand->ClustsInTrack[ijk]->GetYPos()<<" "<<fTrackCand->ClustsInTrack[ijk]->GetZPos()<<" "<<int(fTrackCand->ClustsInTrack[ijk]->GetStraight())<<endl; 
    //    for (unsigned int j=0; j<TrkClustsData[i].size(); j++) {
    //      cout <<TrkClustsData[i][j].XPos<<" "<<TrkClustsData[i][j].YPos<<" "<<TrkClustsData[i][j].ZPos<<" "<< TrkClustsData[i][j].PlaneView<<endl;      
    //    }
  }
  
  for (unsigned ijk=0; ijk<ntrk-iend; ijk++) {
    //  for (unsigned ijk=0; ijk<ntrk; ijk++) {
    //  for (unsigned ijk=0; ijk<ntrk-1; ijk++) {
    if (!fTrackCand->ClustsInTrack[ijk]->GetStraight()) continue;
    int i=fTrackCand->ClustsInTrack[ijk]->GetZPlane();

    //    cout<<"ijk  "<< ijk<<" "<< i<<" "<< TrkClustsData[i].size()<<endl; 
    //    for (unsigned int j=0; j<TrkClustsData[i].size(); j++) {
    //      cout <<"ikcl "<< ijk <<" "<<i<<" "<<TrkClustsData[i][j].XPos<<" "<<TrkClustsData[i][j].YPos<<" "<<TrkClustsData[i][j].ZPos<<" "<< TrkClustsData[i][j].PlaneView<<endl;      
    //    }
    //  for (int i=StartPlane; i<=EndPlane; ++i) {
    EndofRange = false;
    if (TrkClustsData[i].size()>0) {
      if (PassTrack) {
        // Find Next Plane
        int NewPlane=-99;
	//	int k=(i+1);
        //	if (!first && ijk !=ntrk-1) {
	//	  k=fTrackCand->ClustsInTrack[ijk+1]->GetZPlane();
	//	}
	
	for (int ij=0; ij<6; ij++) {x_k_old[ij] = x_k_minus[ij];}
	//a	cout <<"ijlcs "<< ijk <<" "<<i<<" "<< x_k_minus[0]<<" "<<" "<< x_k_minus[1]<<" "<< x_k_minus[2]<<" "<< x_k_minus[3]<<" "<< x_k_minus[4]<<" "<< x_k_minus[5]<<" "<<int(first)<< endl;
	if (!first && ijk !=ntrk-1 && fTrackCand->ClustsInTrack[ijk+1]->GetStraight()) {
	  //	if (!first && ijk !=ntrk-1) {
	  NewPlane=fTrackCand->ClustsInTrack[ijk+1]->GetZPlane();
	} else {
          if (fTrackCand->GetEntries()>=MINLAYER) {
	    int plane = i;
	    int loopmx = ((ijk<ntrk-1) ) ? max(MaxPlane-i,5) : 4;
	    //	    cout <<"loopmx "<< loopmx<<" "<<i<<" "<<MaxPlane<<endl;
	    double dsBefore=0;
	    double drangeBefore=0;
	    
	    for (int nloop=0; nloop <loopmx; nloop++) {
	      if (nloop==0) {
		for(int ij=0; ij<6; ++ij) {StateVector[ij]=x_k_minus[ij];}
	      } else {
		for(int ij=0; ij<6; ++ij) {StateVector[ij]= Prediction[ij];}
	      }
	      int nextplane = -99;
	      double ds=0;
	      double drange=0;
	      if(debug_new) cout <<"before "<< i<<" "<<plane<<" "<<nextplane <<" "<<1/StateVector[4] <<" "<<1/Prediction[4]<<endl;
	      bool GetPrediction=Swim_new(StateVector, Prediction, plane, nextplane, true, &ds, &drange);
	      if(debug_new) 
		cout <<"planecom "<< i<<" "<<plane<<" "<<nextplane<<" "<<int(GetPrediction) <<" "<<1/StateVector[4] <<" "<<1/Prediction[4]<<endl;
	      
	      if (!GetPrediction || nextplane<0 || nextplane>=int(nLayer)) break;
	      if (GetPrediction) {
		if(loopmx<10) {
		  double dirGeom[3] = {0.0,0.0,1.0};
		  double posGeom[3] = {100*Prediction[0],100*Prediction[1],100*ZPosLayer[nextplane]};
		  icalGeometry->InitTrack(posGeom, dirGeom);
		  TGeoMaterial* localmat= icalGeometry->GetCurrentVolume()->GetMaterial();
		  // int matchk = int(strstr(localmat->GetName(),"rpcgas"));
		  if(!(strstr(localmat->GetName(),"rpcgas"))) {
		    loopmx += 10;
		  }
		}
		if (nextplane < plane) {
		  if (TrkClustsData[nextplane+shiftLa].size()>0) {
		    // cout<<"Storing Filtered Data absr 1"<<endl;
		    StoreFilteredData_sr(nextplane+shiftLa, Prediction, false);
		    fTrackCand->f2dS[nextplane+shiftLa] =ds+dsBefore;
		    fTrackCand->f2Range[nextplane+shiftLa] =drange+drangeBefore;
		    dsBefore = drangeBefore = 0;
		  } else {

		    // cout<<"Storing Filtered Data absr 2"<<endl;
		    StoreFilteredData_sr(nextplane+shiftLa, Prediction, false);
		    dsBefore +=ds;
		    drangeBefore +=drange;
		  }
		} else {
		  plane = nextplane;
		  if (TrkClustsData[plane].size() >0) { NewPlane= plane; break;}
		  // cout<<"Storing Filtered Data absr 3"<<endl;
		  StoreFilteredData_sr(nextplane, Prediction, true);
		}

	      }

	    }
	  }
	}

	for (int ij=0; ij<6; ij++) {x_k_minus[ij] = x_k_old[ij];}

	//	cout<<"goforpref "<<i<<" "<<Prediction[4]<<endl;
        if (NewPlane!=-99) {
          // Define measurement function
	  //	  cout <<"forTrkClustsData[NewPlane] "<<NewPlane<<" "<< TrkClustsData[NewPlane].size()<<endl;
	  int PlaneView = TrkClustsData[NewPlane][0].PlaneView;
	  //GMA for Clusts this condition is fine, but for cluster, this is not correct
	  // Carry out the Kalman fit
	  for (int ij=0; ij<2; ij++) {
	    for (int jk=0; jk<5; jk++) {
	      H_k[ij][jk]=0;
	    }
	  }

	  if (PlaneView%2==0) {H_k[0][0]=1;}
	  if (PlaneView   >0) {H_k[1][1]=1;}
	  //cout <<"InoTrackFitAlg.GoForwards : WARNING : PlaneView for clusts is not matching with 0/1/2"<<endl;

	  //	  bool direction = false; if (NewPlane>i) direction= true;
	  //	  bool CombiPropagatorOk=GetCombiPropagator(i,NewPlane,direction);
	  // bool CombiPropagatorOk=GetCombiPropagator(i,NewPlane,true);
	  double ds_gpl = 0.0;
	  double drange = 0.0;
	  bool CombiPropagatorOk=GetCombiPropagator(i,NewPlane,true,&ds_gpl, &drange);
	  
          if(CombiPropagatorOk ) {
            
	    GetNoiseMatrix(i,NewPlane);
            ExtrapCovMatrix();
            CalcKalmanGain(NewPlane);

	    //            UpdateStateVector(i,NewPlane,direction); //true);
	    if(debug_new) cout<<"Go Forwards .... "<<" "<<first<<endl;
	    UpdateStateVector(i,NewPlane,true);
            UpdateCovMatrix();
            MoveArrays();
            if(SaveData) {
	      // cout<<"Storing Filtered Data ab8"<<endl;
	      StoreFilteredData(NewPlane);
	      fTrackCand->fdS[NewPlane] = ds_gpl;
	      fTrackCand->fRange[NewPlane] = drange;
	    }
	    //	    cout <<" InoTrackFitAlg : GoForwards, Filtered q/p " << x_k[4] << endl;

          }
          else {
	    cout <<"combifail "<< i<<" "<<NewPlane<<endl;

	    PassTrack=false;
	  }
        }
        
      }
      //else {cout <<" InoTrackFitAlg : GoForwards, Outside of detector - track failed" <<" "<<ijk<<" "<<i<< endl;}
    }
    //JAM end of range found
    if(EndofRange && LastIteration && ZIncreasesWithTime){
      EndofRangePlane=i;
      break;
    }
  }
  // Store entries from covariance matrix for use in setting track properties
  if(LastIteration) { // 07/02/2009 NIter==2) {
    if(ZIncreasesWithTime==true) {
      EndCov[0]=C_k[0][0]; EndCov[1]=C_k[1][1];
      EndCov[2]=C_k[2][2]; EndCov[3]=C_k[3][3];
      EndCov[4]=C_k[4][4];
    }
    else {
      VtxCov[0]=C_k[0][0]; VtxCov[1]=C_k[1][1];
      VtxCov[2]=C_k[2][2]; VtxCov[3]=C_k[3][3];
      VtxCov[4]=C_k[4][4];
    }
  }
  if(TrkFitterDebug>0) {
    cout <<" InoTrackFitAlg : GoForwards complete "<<first<<endl;
  }
}

void InoTrackFitAlg::GoBackwards_new(bool first) { 
  // Carry out the Kalman fit along the track in the direction of decreasing z
  if(TrkFitterDebug>0) {
    cout<<endl;
    cout <<" InoTrackFitAlg : GoBackwards_new, carry out fit in negative z direction --- " <<first<< endl;
  }
  double StateVector[6]; double Prediction[6];
  Int_t StartPlane = MaxPlane;// Int_t EndPlane=MinPlane;
  if(ZIncreasesWithTime){
    StartPlane = EndofRangePlane;
  }
  else EndofRangePlane = MinPlane;
  double x_k_old[6];
 
  int ntrk = fTrackCand->GetEntries();
  int iend = (first) ? 0 : 1;
  for (int ijk=ntrk-1; ijk>=iend; ijk--) {
    //  for (int ijk=ntrk-1; ijk>=0; ijk--) {
    //  for (int ijk=ntrk-1; ijk>0; ijk--) {
    
    if (!fTrackCand->ClustsInTrack[ijk]->GetStraight()) continue;
    int i=fTrackCand->ClustsInTrack[ijk]->GetZPlane();
    //  for (int i=StartPlane; i>=EndPlane; --i) {  
    if (TrkClustsData[i].size()>0) {
      if (PassTrack) {
	
        //Find Prev Plane
        int NewPlane=-99;
	//        int k=(i-1);

	for (int ij=0; ij<6; ij++) {x_k_old[ij] = x_k_minus[ij];}
	
	//a	cout <<"backfirst "<<int(first)<<" "<<ijk<<" "<<ntrk-1<<endl;
	if (!first && ijk !=0 && fTrackCand->ClustsInTrack[ijk-1]->GetStraight()) {
	  //	if (!first && ijk !=0 ) { 
	  NewPlane=fTrackCand->ClustsInTrack[ijk-1]->GetZPlane();
	} else {
          if (fTrackCand->GetEntries()>=MINLAYER) {
	    int plane = i;
	    //	    int loopmx = ((ijk<ntrk-1) && (ijk>0)) ? i - MinPlane : 3;
	    int loopmx = ((ijk>0)) ? max(i-MinPlane, 5) : 4;
	    double dsBefore=0;
	    double drangeBefore=0;
	    for (int nloop=0; nloop <loopmx; nloop++) {
	      if (nloop==0) {
		for(int ij=0; ij<6; ++ij) {StateVector[ij]=x_k_minus[ij];}
	      } else {
		for(int ij=0; ij<6; ++ij) {StateVector[ij]= Prediction[ij];}
	      }
	      double ds=0;
	      double drange=0;
	      int nextplane = -99;
	      bool GetPrediction=Swim_new(StateVector, Prediction, plane, nextplane, false, &ds, &drange);
	      // if(debug_new) 
	      // cout <<"planecom "<< i<<" "<<plane<<" "<<nextplane<<" "<<int(GetPrediction) <<" "<<1/StateVector[4] <<" "<<1/Prediction[4]<<endl;
	      
	      if (!GetPrediction || nextplane<0 || nextplane>=int(nLayer)) break;
	      if (GetPrediction) {
		if(loopmx<10) { 
		  double dirGeom[3] = {0.0,0.0,1.0};
		  double posGeom[3] = {100*Prediction[0],100*Prediction[1],100*ZPosLayer[nextplane]};
		  icalGeometry->InitTrack(posGeom, dirGeom);
		  TGeoMaterial* localmat= icalGeometry->GetCurrentVolume()->GetMaterial();
		  if(!(strstr(localmat->GetName(),"rpcgas"))) {
		    loopmx += 10;
		  }
		}
		if (nextplane > plane) {
		  if (TrkClustsData[nextplane+shiftLa].size()>0) {
		    // cout<<"Storing Filtered Data absr 4"<<endl;
		    StoreFilteredData_sr(nextplane+shiftLa, Prediction, false);
		    fTrackCand->f2dS[nextplane+shiftLa] =ds+dsBefore;
		    fTrackCand->f2Range[nextplane+shiftLa] =drange+drangeBefore;
		    dsBefore = drangeBefore = 0;
		  } else {
		    // cout<<"Storing Filtered Data absr 5"<<endl;
                    StoreFilteredData_sr(nextplane+shiftLa, Prediction, false);
		    dsBefore +=ds;
		    drangeBefore +=drange;
		  }
		} else {
		  plane = nextplane;
		  if (TrkClustsData[plane].size() >0) {
		    NewPlane= plane;
		    break;
		  }
		  // cout<<"Storing Filtered Data absr 6"<<endl;
		  StoreFilteredData_sr(nextplane, Prediction, true);
		}

	      }
	    }
	  }
	}
        for (int ij=0; ij<6; ij++) {x_k_minus[ij] = x_k_old[ij];}
        if (NewPlane!=-99) {
	  // Define measurement function
	  int PlaneView = TrkClustsData[NewPlane][0].PlaneView;
	  
	  //GMA for Clusts this condition is fine, but for cluster, this is not correct
	  // Carry out the Kalman fit
	  for (int ij=0; ij<2; ij++) {
	    for (int jk=0; jk<5; jk++) {
	      H_k[ij][jk]=0;
	    }
	  }
	  if (PlaneView%2==0) {H_k[0][0]=1;}
	  if (PlaneView   >0) {H_k[1][1]=1;}
	  
	  // cout <<"InoTrackFitAlg.GoBackwards : WARNING : PlaneView for clusts is not matching with 0/1/2"<<endl;
	  // bool direction = false; if (NewPlane>i) direction = true;
	  //	  bool CombiPropagatorOk=GetCombiPropagator(i,NewPlane,direction);
	  double ds_gpl = 0.0;
	  double drange = 0.0;
	  bool CombiPropagatorOk=GetCombiPropagator(i,NewPlane,false,&ds_gpl, &drange);
          if(CombiPropagatorOk ) {
            GetNoiseMatrix(i,NewPlane);
            ExtrapCovMatrix();
            CalcKalmanGain(NewPlane);
	    //            UpdateStateVector(i,NewPlane,direction); // false);
	    // if(debug_new)
	    UpdateStateVector(i,NewPlane,false);
	    UpdateCovMatrix();
            MoveArrays();
	    if(SaveData) {
	      // cout<<"Storing Filtered Data ab9"<<endl;
	      StoreFilteredData(NewPlane);
	      fTrackCand->fdS[NewPlane] = ds_gpl;
	      fTrackCand->fRange[NewPlane] = drange;
	    }
	    //a	    cout <<" InoTrackFitAlg : GoBackwards, Filtered q/p " << x_k[4] << endl;
	    //	    cout <<"xarrayBO "<<i<<" "<<NewPlane<<" "<<x_k_minus[0]<<" "<<x_k_minus[1]<<" "<<x_k_minus[2]<<" "<<x_k_minus[3]<<" "<<x_k_minus[4]<<" "<<x_k_minus[5]<<endl;
          }
          else {cout<<" PassTrack 5" <<" ievt "<<pAnalysis->ievt2<< endl;PassTrack=false;}
        }
      }
      //else {cout <<" InoTrackFitAlg : GoBackwards, Outside detector - track failed" << endl;}
    }
    //JAM end of range found
    if(EndofRange && LastIteration && !ZIncreasesWithTime){ //VALGRIND
      EndofRangePlane = i;
      break;
    }
  }

  // Store entries from covariance matrix for use in setting track properties
  if(LastIteration) { // 07/02/09 NIter==2) {
    if(ZIncreasesWithTime==true) {
      VtxCov[0]=C_k[0][0]; VtxCov[1]=C_k[1][1];
      VtxCov[2]=C_k[2][2]; VtxCov[3]=C_k[3][3];
      VtxCov[4]=C_k[4][4];
    }
    else {
      EndCov[0]=C_k[0][0]; EndCov[1]=C_k[1][1];
      EndCov[2]=C_k[2][2]; EndCov[3]=C_k[3][3];
      EndCov[4]=C_k[4][4];
    }
  }
  if(TrkFitterDebug>0) {
    cout <<" InoTrackFitAlg : GoBackwards complete "<<first<<endl;
  }
}

bool InoTrackFitAlg::DirectionFromFinderHits( InoTrack *trk) {
//  cout<< "  InoTrackFitAlg::DirectionFromFinderHits " <<endl;
  unsigned nhits = trk->GetEntries();
  double sxy=0; // y (distance) = c t + shift
  double sx=0; // y = distance
  double sy=0; //   x = time
  double sn=0;
  double sx2=0;
  double dist=0;
  for (unsigned ij=0; ij<nhits; ij++) {
    if (ij>0) {
      dist += pow( pow(trk->ClustsInTrack[ij]->GetXPos() - trk->ClustsInTrack[ij-1]->GetXPos(), 2.) +
		   pow(trk->ClustsInTrack[ij]->GetYPos() - trk->ClustsInTrack[ij-1]->GetYPos(), 2.) +
		   pow(trk->ClustsInTrack[ij]->GetZPos() - trk->ClustsInTrack[ij-1]->GetZPos(), 2.), 0.5);
    }
    double time = trk->ClustsInTrack[ij]->GetTime();
    sn +=1;
    sx += time; //Look again for return track
    sy  += dist;
    sxy += dist*time;
    sx2 += time*time;   
  }
  double velocity = 0;
  if (sn >0 && (sx2*sn - sx*sx) !=0) {
    velocity = 10*(sxy*sn - sx*sy)/(sx2*sn - sx*sx); // x 10^8m/s
  }
//a  cout <<"velocity "<< velocity <<" "<<dist<<" "<<nhits<<endl;
  return (velocity>0) ? true : false;

}

bool InoTrackFitAlg::DirectionFromFitterHits(InoTrackCand *trk, int epln, double& xslope, double& xintercept, double& xexp) {
//  cout<< "  InoTrackFitAlg::DirectionFromFinderHits " <<endl;
  unsigned nhits = trk->GetEntries();
  double sxy=0; // y (distance) = c t + shift
  double sx=0; // y = distance
  double sy=0; //   x = time
  double sn=0;
  double sx2=0;
  double dist=0;
  vector <double> xval_loc;
  vector <double> yval_loc;
  double xxloc;
  double yyloc;
  double yexp_loc = -10000;
  for (unsigned ij=0; ij<nhits; ij++) {
    if (ij>0) {
      dist += pow( pow(trk->ClustsInTrack[ij]->GetXPos() - trk->ClustsInTrack[ij-1]->GetXPos(), 2.) +
		   pow(trk->ClustsInTrack[ij]->GetYPos() - trk->ClustsInTrack[ij-1]->GetYPos(), 2.) +
		   pow(trk->ClustsInTrack[ij]->GetZPos() - trk->ClustsInTrack[ij-1]->GetZPos(), 2.), 0.5);
    }
    double time = trk->ClustsInTrack[ij]->GetTime();
    xval_loc.push_back(dist);
    yval_loc.push_back(time);
    
    if(trk->ClustsInTrack[ij]->GetZPlane()==epln) {
      xxloc = dist;
      yyloc = time;
    } else {
      sn +=1;
      sx += dist; //Look again for return track
      sy  += time;
      sxy += dist*time;
      sx2 += dist*dist;   
    }
  }
  double velocity = 0;
  if (sn >0 && (sx2*sn - sx*sx) !=0) {
    xslope = (sxy*sn - sx*sy)/(sx2*sn - sx*sx); // ns/m
    xintercept = sy/sn - xslope*sx/sn;
    // cout<<"nhits "<<nhits<<endl;
    yexp_loc = xslope*xxloc + xintercept;
    xexp = yexp_loc;

    for(unsigned kl=0; kl<nhits; kl++) {
      int tninla = trk->ClustsInTrack[kl]->GetZPlane();
      double y1_loc = xslope*xval_loc[kl] + xintercept;
      // cout<<"kl "<<kl<<" "<<trk->ClustsInTrack[kl]->HitsInCluster.size()<<endl;
      for(unsigned jk=0; jk<trk->ClustsInTrack[kl]->HitsInCluster.size(); jk++) {
	// if(trk->ClustsInTrack[kl]->HitsInCluster.size()==1) {
	// if(jk==0) {
	pAnalysis->hdifftime2[tninla]->Fill(trk->ClustsInTrack[kl]->HitsInCluster[jk]->GetXTimeCorr() - trk->ClustsInTrack[kl]->HitsInCluster[jk]->GetYTimeCorr());
	// }
	pAnalysis->hxtime_ext[tninla]->Fill(trk->ClustsInTrack[kl]->HitsInCluster[jk]->GetXTimeCorr() - y1_loc);
	pAnalysis->hytime_ext[tninla]->Fill(trk->ClustsInTrack[kl]->HitsInCluster[jk]->GetYTimeCorr() - y1_loc);
	pAnalysis->h_hit_time_ext[tninla]->Fill(trk->ClustsInTrack[kl]->HitsInCluster[jk]->GetTime() - y1_loc);
	
      }
    }
  } else {
    xslope = -10000;
    xintercept = -10000;
    xexp = -10000;
  }
  
  //a  cout <<"velocity "<< velocity <<" "<<dist<<" "<<nhits<<endl;
  return (xslope>0) ? true : false;
  
}

int InoTrackFitAlg::CheckFCPCUpOrDn(double *aax_k, bool DirExtraPol, int MaxMinPlane, bool GoDir) {
  //bit 1 : LocalPos[1];
  //bit 2 : CheckMat[1]
  //bit 3 : LocalPos[0]
  //bit 4 : CheckMat[0]
  
  if(debug_fcpc) cout<<"int InoTrackFitAlg::CheckFCPCUpOrDn(double *"<<aax_k[0]<<","<<aax_k[1]<<", bool "<<DirExtraPol<<", int "<<MaxMinPlane<<", bool "<<GoDir<<") {...."<<endl;
  double startStateVector[5];

  const int nEmpyExtr=2;
  const int nMxShowerlayer=6;
  int shwLyrMx = (DirExtraPol == ZIncreasesWithTime) ? 2 : nMxShowerlayer;
  unsigned int CheckMat[nEmpyExtr] = {0};
  unsigned int LocalPos[nEmpyExtr] = {0};
  
  double dir[3]={0.0,0.0,1.0};
  int nextplane =0;
  int startplane = -1;
  double Prediction[6]	= {0.0};

  double tmpaax_k[5];
  double tmpaC_k_intermediate[5][5];
  
  for (int ij=0; ij<5; ij++) {
    tmpaax_k[ij] = aax_k[ij];
    for (int jk=0; jk<5; jk++) {
      tmpaC_k_intermediate[ij][jk] = C_k_intermediate[ij][jk];
    }
  }

  //DirExtraPol = true for upward extrapolation aka Plane Increases
  //DirExtraPol = false for downward extrapolation aka Plane Decreases

  int iext=0;
  int iempty=0;
  while(iext<shwLyrMx && iempty<nEmpyExtr) {
    double dsExt = 0.0;
    double drangeExt = 0.0;
    //nextplane = -10;
    if(DirExtraPol) {
      startplane = MaxMinPlane + iext;
      nextplane = startplane + 1;
      if (nextplane >= int(nLayer)) { break;} //{nextplane = nLayer-1;}
      if(debug_fcpc) cout<<"Extrapolating upwards aka plane increases... "<<startplane<<" "<<nextplane<<endl;
    } else {
      startplane = MaxMinPlane - iext;
      nextplane = startplane - 1;
      if (nextplane <0) { break;} // {nextplane = 0;}
      if(debug_fcpc) cout<<"Extrapolating downwards aka plane decreases... "<<startplane<<" "<<nextplane<<endl;
    }
    iext++;
    // if(debug_fcpc) cout<<"debug_fcpc "<<debug_fcpc<<endl;
    
    bool GetPrediction=Swim(aax_k, Prediction, startplane, nextplane, GoDir);//, &dsExt, &drangeExt);
    // bool GetPrediction= PredictedStateCov(aax_k, startplane, nextplane, GoDir, Prediction, 0, &dsExt, &drangeExt);
    if(debug_fcpc) cout<<"SwimExtrapolate("<<aax_k[0]<<", "<<aax_k[1]<<", "<<Prediction[0]<<", "<<Prediction[1]<<", "<<startplane<<", "<<nextplane<<", "<<GoDir<<", &dsExt, &drangeExt); = "<<GetPrediction<<endl;
    if (GetPrediction) {
      double pos[3] = {100*Prediction[0],100*Prediction[1],100*ZPosLayer[nextplane]};
      double dsdz_fcpc_ini = sqrt(1.0 + aax_k[2]*aax_k[2] + aax_k[3]*aax_k[3]);
      if(ZIncreasesWithTime==false) {dsdz_fcpc_ini = - dsdz_fcpc_ini;}
      double theta_fcpc_ini = 180*acos(1./dsdz_fcpc_ini)/3.1415;
      double phi_fcpc_ini = 180*atan2(aax_k[3],aax_k[2])/3.1415;
      double momentum_fcpc_ini = 0;
      if(aax_k[4]!=0) {
	momentum_fcpc_ini = 1./aax_k[4];
      }

      double dsdz_fcpc_fin = sqrt(1.0 + Prediction[2]*Prediction[2] + Prediction[3]*Prediction[3]);
      if(ZIncreasesWithTime==false) {dsdz_fcpc_fin = - dsdz_fcpc_fin;}
      double theta_fcpc_fin = 180*acos(1./dsdz_fcpc_fin)/3.1415;
      double phi_fcpc_fin = 180*atan2(Prediction[3],Prediction[2])/3.1415;
      double momentum_fcpc_fin = 0;
      if(Prediction[4]!=0) {
	momentum_fcpc_fin = 1./Prediction[4];
      }
      
      if(debug_fcpc)  {
	cout<<"Initial Position Info: momentum = "<<momentum_fcpc_ini<<", cos(theta) = "<<1./dsdz_fcpc_ini<<" theta = "<<theta_fcpc_ini<<" phi = "<<phi_fcpc_ini<<" x = "<<100*aax_k[0]<<" y = "<<100*aax_k[1]<<" z = "<<100*ZPosLayer[startplane]<<" Plane No. "<<startplane<<endl;
	cout<<"Final Position Info: momentum = "<<momentum_fcpc_fin<<", cos(theta) = "<<1./dsdz_fcpc_fin<<" theta = "<<theta_fcpc_fin<<" phi = "<<phi_fcpc_fin<<" x = "<<pos[0]<<" y = "<<pos[1]<<" z = "<<pos[2]<<" Plane No. "<<nextplane<<endl;
      }
      
      icalGeometry->InitTrack(pos, dir);
      double localpos[3];
      icalGeometry->MasterToLocal(pos, localpos);
      TGeoMaterial* localmat= icalGeometry->GetCurrentVolume()->GetMaterial();
      // cout<<"localmat "<<localmat->GetName()<<endl;
      CheckMat[iempty] = 0;
      iempty++;
      if (strstr(localmat->GetName(),"rpcgas")) {
	CheckMat[iempty-1] = 1;
	LocalPos[iempty-1] = 1;
	if(debug_fcpc) {
	  cout<<".............................................."<<endl;
	  cout<<"Fully Contained..."<<DirExtraPol<<"..."<<GoDir<<"...NextPlane..."<<nextplane<<endl;
	  cout <<"pos "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<endl;
	  cout <<"pos "<<localpos[0]<<" "<<localpos[1]<<" "<<localpos[2]<<endl;
	  cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
	}
	
      } // if (strstr(localmat->GetName(),"rpcgas")) {
      for (int ix=0; ix<5; ix++) {
	aax_k[ix] = Prediction[ix];
      }
    } else { // else for if (GetPrediction) {
      double pos[3] = {100*Prediction[0],100*Prediction[1],100*ZPosLayer[nextplane]};
      double dsdz_fcpc_ini = sqrt(1.0 + aax_k[2]*aax_k[2] + aax_k[3]*aax_k[3]);
      if(ZIncreasesWithTime==false) {dsdz_fcpc_ini = - dsdz_fcpc_ini;}
      double theta_fcpc_ini = 180*acos(1./dsdz_fcpc_ini)/3.1415;
      double phi_fcpc_ini = 180*atan2(aax_k[3],aax_k[2])/3.1415;
      double momentum_fcpc_ini = 0;
      if(aax_k[4]!=0) {
	momentum_fcpc_ini = 1./aax_k[4];
      }

      double dsdz_fcpc_fin = sqrt(1.0 + Prediction[2]*Prediction[2] + Prediction[3]*Prediction[3]);
      if(ZIncreasesWithTime==false) {dsdz_fcpc_fin = - dsdz_fcpc_fin;}
      double theta_fcpc_fin = 180*acos(1./dsdz_fcpc_fin)/3.1415;
      double phi_fcpc_fin = 180*atan2(Prediction[3],Prediction[2])/3.1415;
      double momentum_fcpc_fin = 0;
      if(Prediction[4]!=0) {
	momentum_fcpc_fin = 1./Prediction[4];
      }
      if(debug_fcpc)  {
	cout<<"ELSE CONDITION::::"<<endl;
	cout<<"Initial Position Info: momentum = "<<momentum_fcpc_ini<<", cos(theta) = "<<1./dsdz_fcpc_ini<<" theta = "<<theta_fcpc_ini<<" phi = "<<phi_fcpc_ini<<" x = "<<100*aax_k[0]<<" y = "<<100*aax_k[1]<<" z = "<<100*ZPosLayer[startplane]<<" Plane No. "<<startplane<<endl;
	cout<<"Final Position Info: momentum = "<<momentum_fcpc_fin<<", cos(theta) = "<<1./dsdz_fcpc_fin<<" theta = "<<theta_fcpc_fin<<" phi = "<<phi_fcpc_fin<<" x = "<<pos[0]<<" y = "<<pos[1]<<" z = "<<pos[2]<<" Plane No. "<<nextplane<<endl;
      }

      if(abs(momentum_fcpc_ini)<0.1) {
	CheckMat[0] = 1;
	LocalPos[0] = 1;
	CheckMat[1] = 1;
	LocalPos[1] = 1;
      } else {
	if(CheckMat[0]) {
	  CheckMat[0] = 0;
	  LocalPos[0] = 0;
	}
	CheckMat[1] = 0;
	LocalPos[1] = 0;
      }
      goto ENDEXTRAPOL;
    } // if (GetPrediction) {
  } // while(iext<nMxShowerlayer && iempty<nEmpyExtr) {
  
 ENDEXTRAPOL: 
  
  for (int ij=0; ij<5; ij++) {
    aax_k[ij] = tmpaax_k[ij];
    for (int jk=0; jk<5; jk++) {
      C_k_intermediate[ij][jk] = tmpaC_k_intermediate[ij][jk];
    }
  }


  unsigned int alltags = 0;
  alltags +=8*CheckMat[0]+4*LocalPos[0]+2*CheckMat[1]+LocalPos[1];
  if(debug_fcpc)  cout <<"alltags "<<alltags<<endl;
  if(debug_fcpc)  cout<<"...}"<<endl;
  return alltags;
}

void InoTrackFitAlg::StraightLineFit(vector<Hep3Vector>& finder, vector<int>& loczlay, int eplnstrt, double* locslope, double* locintercpt, double* locchi2, int* locnhit, double* locexpecpos) {

  const int size = finder.size();
  double szx=0, sz=0, sx=0, snx=0, sz2=0, szy=0, sy=0, sny=0;
  double xerr2 = pow(StripXWidth/pow(12.,0.5),2);
  double yerr2 = pow(StripYWidth/pow(12.,0.5),2);

  // cout<<xerr2<<" "<<yerr2<<endl;
  
  for (int ij=0; ij<size; ij++) {
    szx +=finder[ij].z()*finder[ij].x()/xerr2;
    szy +=finder[ij].z()*finder[ij].y()/yerr2;
    sz +=finder[ij].z()/xerr2;
    sz2 +=finder[ij].z()*finder[ij].z()/xerr2;
    sx +=finder[ij].x()/xerr2;
    sy +=finder[ij].y()/yerr2;
    snx +=1/xerr2;
    sny +=1/yerr2;
  }
  
  if (snx >0 && sz2*snx - sz*sz !=0) { 
    double slope = locslope[0] = (szx*snx - sz*sx)/(sz2*snx - sz*sz);
    double intersect = locintercpt[0] = sx/snx - slope*sz/snx;  // <x> = a + b <z>
    double determ = (snx*sz2 - sz*sz);
    int nused = 0;
    double chiSquare = 0.0;

    for(int ij=0; ij<size; ij++) {
      double xexpec = slope*finder[ij].z() + intersect;
      chiSquare +=  pow(finder[ij].x()-xexpec,2.)/xerr2;
      // cout<<"chis2 ijx "<<ij<<" "<<finder[ij].x()<<" "<<xexpec<<" "<<pow(finder[ij].x()-xexpec,2.)<<" "<<pow(finder[ij].x()-xexpec,2.)/xerr2<<endl;
      pAnalysis->hxpos_ext[loczlay[ij]]->Fill((finder[ij].x()-xexpec)/StripXWidth);
      nused++;
    }
    locchi2[0] = chiSquare;
    locnhit[0] = nused;
    //    cout<<"Straight "<< slope <<" "<<intersect <<" "<<errcst<<" "<<errcov<<" "<<errlin <<endl;

    if(eplnstrt-1<0 || eplnstrt>9) {
      locexpecpos[0] = -100000.0;
    } else {
      locexpecpos[0] = locslope[0]*ZPosLayer[eplnstrt-1] + locintercpt[0];
    }
    
  } else {
    locchi2[0] = locslope[0] = locintercpt[0] = -1; locnhit[0] = -1; locexpecpos[0] = -100000.0;
  }

  if (sny >0 && sz2*sny - sz*sz !=0) { 
    double slope = locslope[1] = (szy*sny - sz*sy)/(sz2*sny - sz*sz);
    double intersect = locintercpt[1] = sy/sny - slope*sz/sny;  // <y> = a + b <z>
    double determ = (sny*sz2 - sz*sz);

    int nused = 0;
    double chiSquare = 0.0;
    
    for(int ij=0; ij<size; ij++) {
      double yexpec = slope*finder[ij].z() + intersect;
      chiSquare +=  pow(finder[ij].y()-yexpec,2.)/yerr2;
      // cout<<"chis2 ijy "<<ij<<" "<<finder[ij].y()<<" "<<yexpec<<" "<<pow(finder[ij].y()-yexpec,2.)<<" "<<pow(finder[ij].y()-yexpec,2.)/yerr2<<endl;
      pAnalysis->hypos_ext[loczlay[ij]]->Fill((finder[ij].y()-yexpec)/StripYWidth);
      nused++;
    }
    locchi2[1] = chiSquare;
    locnhit[1] = nused;

    if(eplnstrt-1<0 || eplnstrt>9) {
      locexpecpos[1] = -100000.0;
    } else {
      locexpecpos[1] = locslope[1]*ZPosLayer[eplnstrt-1] + locintercpt[1];
    }
    
  } else {
    locchi2[1] = locslope[1] = locintercpt[1] = -1; locnhit[1] = -1; locexpecpos[1] = -100000.0;
  }
  
}

void InoTrackFitAlg::simple_track_fit(vector<Hep3Vector>& finder, int elpnsimple, double& curve, double& radii, double* x0, double* locchisq, double* xavg, int& locnhits, double& locexpecpos) {
  double xerr2 = pow(StripXWidth/pow(12.,0.5),2);
  double yerr2 = pow(StripYWidth/pow(12.,0.5),2);
  errxy2 = xerr2;
  int size = finder.size();
  
  double x1 = 0.5*(finder[0].x()+finder[1].x());
  double y1 = 0.5*(finder[0].z()+finder[1].z());
  double x2 = 0.5*(finder[size-1].x()+finder[size-2].x());
  double y2 = 0.5*(finder[size-1].z()+finder[size-2].z());

  int sign =(sfdph(atan2(y1,x1), atan2(y2,x2))>0.0) ? -1 : +1;
  double slope = -x1/y1; //Perpendicular to the initial direction

  double xx = (x2*x2 - x1*x1 + y2*y2 - y1*y1)/(2*(slope*(y2-y1) + x2-x1));
  double yy = slope*xx;
  double radius = min(228.0, sqrt( (xx-x1)*(xx-x1) + (yy-y1)*(yy-y1)))/1.1; 

  if (sqrt(xx*xx+yy*yy) > radius) {
    xx = xx*radius/sqrt(xx*xx+yy*yy);
    yy = yy*radius/sqrt(xx*xx+yy*yy);
  }

  curve = radius; //momentum using simple curvature

  const int nsgpr=3;
  double fitres[nsgpr];
  double parerr[nsgpr];
  
  TMinuit *gMinuit = new TMinuit(nsgpr);
  gMinuit->SetPrintLevel(-1);

  TString hname[nsgpr] = {"x0", "y0", "radius"};
  double strt[nsgpr] = {xx, yy, radius};
  double alow[nsgpr] = {xx-100.0, yy-100.0, 0.25*radius};
  double ahig[nsgpr] = {xx+100.0, yy+100.0, 14*radius};
  double step[nsgpr] = {0.001, 0.001, 0.01};
  
  gMinuit->SetFCN(fcnsg); //Give input function of TMinuit
  
  double arglist[10];
  int ierflg = 0;
  arglist[0] =  1 ;
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
  
  for (int jk=0; jk<nsgpr; jk++) {
    gMinuit->mnparm(jk, hname[jk], strt[jk], step[jk], alow[jk], ahig[jk],ierflg);
  }

  arglist[0] = 0;

  nsize1 = size;
  for(int ij=0; ij<size; ij++) {
    xvalin[ij] = finder[ij].x();
    yvalin[ij] = finder[ij].z();
  }
  
  gMinuit->mnexcm("MINIMIZE", arglist, 0, ierflg);
  
  TString chnam;
  double parv,err,xlo,xup, plerr, mierr, eparab, gcc;
  int iuit;
  
  for (int ij=0; ij<nsgpr; ij++) {
    gMinuit->mnpout(ij, chnam, parv, err, xlo, xup, iuit);
    gMinuit->mnerrs(ij, plerr, mierr, eparab, gcc);
    fitres[ij] = parv;
    parerr[ij] = err;
  }
  
  x0[0] = fitres[0]; 
  x0[1] = fitres[1]; 
  radii = fitres[2]; 

  double x22[1];
  double fval_r;
  Int_t errflga;
  Int_t npar1 = 3;
  fcnsg(npar1,x22,fval_r,fitres,errflga);
  double loc_chisq_new = fval_r;
  double loc_chisq_all = 0;
  double loc_chisq_pos = 0;
  double loc_chisq_neg = 0;
  double fitxavg_pos=0;
  double fitxavg_neg=0;
  double xavg_meas=0;
  int nused = 0;
  for(int ij=0; ij<size; ij++) {
    double loc_val = pow((fitres[2]*fitres[2])-((yvalin[ij]-fitres[1])*(yvalin[ij]-fitres[1])),0.5);
    double xexpec = fitres[0] + loc_val;
    double xexpec2 = fitres[0] - loc_val;
    loc_chisq_pos += pow(xvalin[ij]-xexpec,2.)/xerr2;
    fitxavg_pos += xexpec;
    loc_chisq_neg += pow(xvalin[ij]-xexpec2,2.)/xerr2;
    fitxavg_neg += xexpec2;
    loc_chisq_all += min(pow(xvalin[ij]-xexpec,2.),pow(xvalin[ij]-xexpec2,2.))/xerr2;
    xavg_meas += xvalin[ij];
    nused++;
  }

  locchisq[0] = loc_chisq_pos;
  xavg[0] = fitxavg_pos/nused;
  locchisq[1] = loc_chisq_neg;
  xavg[1] = fitxavg_neg/nused;
  if(loc_chisq_pos>loc_chisq_neg) {
    if(elpnsimple-1<0 && elpnsimple>9) {
      locexpecpos = -100000.0;
    } else {
      locexpecpos = fitres[0] - pow((fitres[2]*fitres[2])-((ZPosLayer[elpnsimple-1]-fitres[1])*(ZPosLayer[elpnsimple-1]-fitres[1])),0.5);
    }
    // locchisq[2] = loc_chisq_neg;
    xavg[2] = fitxavg_neg/nused;
  } else {
    if(elpnsimple-1<0 && elpnsimple>9) {
      locexpecpos = -100000.0;
    } else {
      locexpecpos = fitres[0] + pow((fitres[2]*fitres[2])-((ZPosLayer[elpnsimple-1]-fitres[1])*(ZPosLayer[elpnsimple-1]-fitres[1])),0.5);
    }
    // locchisq[2] = loc_chisq_pos;
    xavg[2] = fitxavg_pos/nused;
  }
  locchisq[2] = loc_chisq_all;
  xavg[3] = xavg_meas;
  locnhits = nused;

  if (gMinuit) { delete gMinuit; gMinuit=0;}
}
