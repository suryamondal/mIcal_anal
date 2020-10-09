#include "InoNewFitAlg.h"

InoNewFitAlg::InoNewFitAlg() {
  int VtxPlane = -100;
  for(int ij=0; ij<nvectormx; ij++) {
    input_x_k[ij] = -10000.;
    for(int kl=0; kl<nvectormx; kl++) {
      input_x_k_err[ij][kl] = -10000.;
    }
    finalStateVector[ij] = -10000.;
    finalSVerror[ij] = -10000.;
    for(int jk=0; jk<2*nlayermx; jk++) {
      AMatrix[jk][ij] = -10000.;
    }
  }
  ZIncreasesWithTime = false;
  for(int ij=0; ij<2*nlayermx; ij++) {
    datapts[ij] = -100000.;
  }
  

  for(int ij=0; ij<2*nlayermx; ij++) {
    exppts[ij] = -100000.;
    for(int jk=0; jk<2*nlayermx; jk++) {
      CovMatrix[ij][jk] = -100000.;
    }
  }
  fgood = 0;
}

InoNewFitAlg::~InoNewFitAlg() {

}

InoNewFitAlg::InoNewFitAlg(double* psvVtx, double* mpts,double* mptserr, double* mptsz, int* occu, int nmiss, int vtxp, bool TrkDir) {
  ZIncreasesWithTime = TrkDir;
  int VtxPlane = vtxp;
  for(int ij=0; ij<nvectormx; ij++) {
    input_x_k[ij] = psvVtx[ij];
    finalStateVector[ij] = -10000;
    finalSVerror[ij] = -10000;
  }
  int lm =0;
  for(int ij=0; ij<5; ij++) {
    for(int jk=0; jk<2*nlayermx; jk++) {
      AMatrix[jk][ij] = 0.0;
    }
  }

  for(int ij=0; ij<nmiss; ij++) {
    occulyr.push_back(occu[ij]);
    // cout<<occu[ij]<<",";
  }
  // cout<<endl;
  // cout<<"occulyr "<<occulyr.size();
  for(int ij=0; ij<2*nlayermx; ij++) {
    posin[ij] = mpts[ij];
    posinerr2[ij] = mptserr[ij];
  }
  for(int ij=0; ij<nlayermx; ij++) {
    posinz[ij] = mptsz[ij];
  }
  

  for(int ij=0; ij<2*nlayermx; ij++) {
    exppts[ij] = 0.0;
    for(int jk=0; jk<2*nlayermx; jk++) {
      CovMatrix[ij][jk] = 0.0;
    }
  }

  Ical0DetectorParameterDef* paradef = Ical0DetectorParameterDef::AnPointer;
  LayerThickness = (1/m)*2*(paradef->GetParlay(2)+paradef->GetParirlay(2));
  nLayer      = paradef->GetnLayer();
  ShiftInX = (1/m)*paradef->GetRPCShift(0);
  ShiftInY = (1/m)*paradef->GetRPCShift(1);
  ShiftInZ = (1/m)*paradef->GetRPCShift(2);
  for (int ijk=0; ijk<nLayer; ijk++) {
    ZPosLayer[ijk] = (1/m)*(-paradef->GetParino(2) + 2*(paradef->GetParhcoil(2)+paradef->GetParcoilsupport(2)) + 2*(ijk+1)*paradef->GetParirlay(2) + (2*ijk+1)*(paradef->GetParlay(2))) + ShiftInZ; //AAR need to correct this   8/12/2011 for layer thickness 5.6cm this is fine
    //    cout<<"ijk "<< ijk<<" "<<ZPosLayer[ijk]<<endl;
  }

  debug_fit = false;
}

bool InoNewFitAlg::DoIterations() {

  MultiSimAnalysisDigi* pAnalysis = MultiSimAnalysisDigi::AnPointer;

  for(int nint = 0; nint<20; nint++) {
    TMatrixD Pi_TMatrix(5,1);

    for(int ij=0; ij<5; ij++) {
      if(nint==0) {Pi_TMatrix[ij][0] = input_x_k[ij];}
      else {Pi_TMatrix[ij][0] = vtx_x_k[ij];}
    }
    if(nint==0) {
      // Pi_TMatrix[4][0] = -0.2;
      // cout<<"Primary input matrix "<<endl;
      // Pi_TMatrix.Print();
    }
    for(int ij=0; ij<6; ij++) {
      if(nint==0) {vtx_x_k[ij] = input_x_k[ij];}
    };
    if(ZIncreasesWithTime) {GoFrd = true;}
    else {GoFrd = false;}

    for(int pij=0; pij<2; pij++) {
      for(int oij=0; oij<5; oij++) {
	AMatrix[pij][oij] = 0.0;
      }
    }
    AMatrix[0][0] = 1;
    AMatrix[1][1] = 1;
    exppts[0] = input_x_k[0];
    exppts[1] = input_x_k[1];
    for(int mno=0; mno<2*nlayermx; mno++) {
      CovMatrix[mno][mno] = posinerr2[mno];
    }
    int ilay = 0;
    // CovMatrix[0][0] = 1;
    // CovMatrix[1][1] = 1;
    for(int nextplane=1; nextplane<10; nextplane++) {
      //   int nextplane = ilay+1;
      int DeltaPlane=abs(nextplane-ilay);
      double DeltaZ=fabs(DeltaPlane*LayerThickness);
        
      double PState[6];  
      double NState[6];  
      double StateVector[6];
      double Prediction[6];
      double Increment=0.01;
      bool SwimInc=false; bool SwimDec=false;
      for(int ij=0; ij<6; ij++) {StateVector[ij] = vtx_x_k[ij];}
      bool GetPrediction = Swim(StateVector, Prediction, ilay, nextplane, true);
      for(int ij=0; ij<6; ij++) {x_k_minus[ij] = StateVector[ij];}
      // cout<<"GetPrediction "<<GetPrediction<<" "<<nextplane<<endl;
      exppts[2*nextplane] = Prediction[0];
      exppts[2*nextplane+1] = Prediction[1];
      for(unsigned int ghi=0; ghi<occulyr.size();ghi++) {
	if(nextplane==occulyr[ghi]) {
	  posin[2*nextplane] = Prediction[0];
	  posin[2*nextplane+1] = Prediction[1];
	  CovMatrix[2*nextplane][2*nextplane] = 1;
	  CovMatrix[2*nextplane+1][2*nextplane+1] = 1;
	}
      }
      //      if(GetPrediction) {
      	GetNoiseMatrix(ilay,nextplane);
	CovMatrix[2*nextplane][2*nextplane] += Q_k_minus[0][0];
	CovMatrix[2*nextplane][2*nextplane+1] += Q_k_minus[0][1];
	CovMatrix[2*nextplane+1][2*nextplane] += Q_k_minus[1][0];
	CovMatrix[2*nextplane+1][2*nextplane+1] += Q_k_minus[1][1];
	for(int djk=0; djk<5; djk++) {
	  Increment = 0.05*fabs(vtx_x_k[djk]);
	  if(Increment<.01) {Increment=.01;}
	  for(int kl=0; kl<6; kl++) {StateVector[kl] = vtx_x_k[kl];}
	  StateVector[djk] += Increment;
	  SwimInc=Swim(StateVector, PState, ilay, nextplane, GoFrd);
	  StateVector[djk] = vtx_x_k[djk];
	  StateVector[djk] -= Increment;
	  SwimDec=Swim(StateVector, NState, ilay, nextplane, GoFrd);
	  AMatrix[2*nextplane][djk] = (PState[0]-NState[0])/(2*Increment);
	  AMatrix[2*nextplane+1][djk] = (PState[1]-NState[1])/(2*Increment);
	  // cout<<"AMatrix "<<" "<<nextplane<<" "<<djk<<" "<<AMatrix[2*nextplane][djk]<<" "<<AMatrix[2*nextplane+1][djk]<<endl;
	} // for(int djk=0; djk<5; djk++) 
	//  } // if(GetPrediction) {
    
      // for(int kl=0; kl<6; kl++) {vtx_x_k[kl] = Prediction[kl];}
    } // for(int ilay=0; ilay<10; ilay++) 

    TMatrixD A_TMatrix(20,5);
    TMatrixD V_TMatrix(20,20);
    TMatrixD P0_TMatrix = Pi_TMatrix;
    TMatrixD D_TMatrix(20,1);
    TMatrixD F_TMatrix(20,1);
  
    for(int row=0; row<20; row++) {
      for(int acol=0; acol<5; acol++) {
	A_TMatrix[row][acol] = 100*AMatrix[row][acol];
      }
      for(int vcol=0; vcol<20; vcol++) {
	V_TMatrix[row][vcol] = 10000*CovMatrix[row][vcol];
      }
      D_TMatrix[row][0] = 100*posin[row];
      F_TMatrix[row][0] = 100*exppts[row];
    }

    // cout<<"A_TMatrix "<<endl;
    // A_TMatrix.Print();
    // cout<<"V_TMatrix "<<endl;
    // V_TMatrix.Print();
    // cout<<"D_TMatrix "<<endl;
    // D_TMatrix.Print();
    // cout<<"F_TMatrix "<<endl;
    // F_TMatrix.Print();

    Double_t det1;
    TMatrixD VI_TMatrix = V_TMatrix;
    VI_TMatrix.Invert(&det1);
    // cout<<"nint1 = "<<nint<<" "<<det1<<endl;
    if(det1 !=0) {
      // cout<<"VI_TMatrix "<<endl;
      // VI_TMatrix.Print();
      fgood++;
    } else {
      V_TMatrix.Print();
      cout<<"nint1 = "<<nint<<" "<<det1<<endl;
      fgood = 0;
      return false;
    }
    // A_TMatrix.Print();
    TMatrixD AT_TMatrix(TMatrixD::kTransposed,A_TMatrix);
    // cout<<"AT_TMatrix "<<endl;
    // AT_TMatrix.Print();
    TMatrixD Product1M = AT_TMatrix*VI_TMatrix;
    // cout<<"AT*VI "<<endl;
    // Product1M.Print();
    TMatrixD Product2M = Product1M*A_TMatrix;
    TMatrixD tmpProd2M = Product2M;
    // cout<<"AT*VI*A "<<endl;
    // tmpProd2M.Print();
    // Product2M.Print();
    Double_t det2;
    Product2M.Invert(&det2);
    // cout<<"nint2 = "<<nint<<" "<<det2<<endl;
    if(det2 !=0) {
      // cout<<"(AT*VI*A).inverse"<<endl;
      // Product2M.Print();
      fgood++;
    } else {
      A_TMatrix.Print();
      tmpProd2M.Print();
      cout<<"nint2 = "<<nint<<" "<<det2<<endl;
      fgood = 0;
      return false;
    }
    TMatrixD Product3M = Product2M*AT_TMatrix;
    TMatrixD Product4M = Product3M*VI_TMatrix;
    // cout<<"(AT*VI*A).inverse()*AT*VI"<<endl;
    // Product4M.Print();
    TMatrixD Product5M = D_TMatrix - F_TMatrix;
    // cout<<"m-f matrix"<<endl;
    // Product5M.Print();
    TMatrixD Product6M = Product4M*Product5M;
    // cout<<"full product matrix"<<endl;
    // Product6M.Print();
    TMatrixD FinalProduct = P0_TMatrix + Product6M;
    // cout<<"P matrix after the iteration: "<<endl;
    // pAnalysis->hkalmanx[1][nint]->Fill(FinalProduct[0][0]-(0.01*D_TMatrix[0][0]));
    // pAnalysis->hkalmany[1][nint]->Fill(FinalProduct[1][0]-(0.01*D_TMatrix[1][0]));

    // FinalProduct.Print();
    for(int fij=0; fij<5; fij++) {
      vtx_x_k[fij] = FinalProduct[fij][0];
    }
    // cout<<"nint3 = "<<nint<<" p = "<<1/vtx_x_k[4]<<endl;
  }
  for(int iij=0; iij<5; iij++) {
    finalStateVector[iij] = vtx_x_k[iij];
  }

  return true;
}
bool InoNewFitAlg::Swim(double* StateVector, double* Output, const int Plane, int& NewPlane, const bool GoForward) {

  if(debug_fit)   cout <<" ----------------------InoTrackFitAlg : Swim ----------------------"<< Plane<<" "<<NewPlane<<endl;
  
  // Initialisations
  // customize for bfield scaling.
  //  BField * bf = new BField(*vldc,-1,0);
  SwimSwimmer* myswimmer = new SwimSwimmer(fabs(LayerThickness*(Plane-NewPlane)), 0.5*LayerThickness); //*vldc,bf);
  //if(debug_fit)
  // cout <<" InoTrackFitAlg : Swim "<< Plane<<" "<<NewPlane<<endl;
  // cout <<" InoTrackFitAlg : Swim StateVector "<<StateVector[0]<<" "<<StateVector[1]<<" "<<StateVector[2]<<" "<<StateVector[3]<<" "<<StateVector[4]<<" "<<endl;
  double charge = 0.;
  bool done = false;
    
  if(fabs(StateVector[4])>1.e-10) {

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

    Int_t nbfield;
    Double_t bave;
    Double_t t_bave;
    
    //    cout <<"charge === "<<charge<<" "<<momentum.X()<<" "<<momentum.Y()<<" "<<momentum.Z()<<" "<<position.X()<<" "<<position.Y()<<" "<<position.Z()<<" "<<dxdz<<" "<<dydz<<" "<<dsdz<<" st "<<StateVector[0]<<" "<<StateVector[1]<<" "<<StateVector[2]<<" "<<StateVector[3]<<" "<<StateVector[4]<<" "<<muon.GetMomentum().Z()<<endl;
    // cout<<"SwimFunc: Plane "<<Plane<<" "<<StateVector[0]<<" "<<StateVector[1]<<" "<<muon.GetMomentumModulus()<<endl;
    // Do the swim, accounting for direction of motion w.r.t time too
    if( (GoForward==true && ZIncreasesWithTime==true)  || (GoForward==false && ZIncreasesWithTime==false) ) {
      done = myswimmer->SwimForward(muon,t_bave);
    } else if( (GoForward==true && ZIncreasesWithTime==false)  || (GoForward==false && ZIncreasesWithTime==true) ) {
      done = myswimmer->SwimBackward(muon,t_bave);
    }
 
     
    bave += t_bave;  //AAR:  added
    nbfield++;       //AAR:  added
    if(done==true) {

      double angle11 = 0;
      double ct=cos(angle11);
      double st=sin(angle11);
      if(muon.GetDirection().Z()!=0. && muon.GetMomentumModulus()!=0.) {
        Output[0]=(st*muon.GetPosition().Y()+ct*muon.GetPosition().X());
        Output[1]=(ct*muon.GetPosition().Y()-st*muon.GetPosition().X());
        Output[2]=(st*(muon.GetDirection().Y()/muon.GetDirection().Z())+ct*(muon.GetDirection().X()/muon.GetDirection().Z()));
        Output[3]=(ct*(muon.GetDirection().Y()/muon.GetDirection().Z())-st*(muon.GetDirection().X()/muon.GetDirection().Z()));
        Output[4]=muon.GetCharge()/muon.GetMomentumModulus();
	// cout<<"SwimFunc: NewPlane "<<NewPlane<<" "<<Output[0]<<" "<<Output[1]<<" "<<muon.GetMomentumModulus()<<endl;

	Output[5]= StateVector[5];
        // Get range and dS from the Swimmer
	
      } else {done=false;}
    }
  
  } else {

    // If infinite momentum, use straight line extrapolation
    
    
    double delz = fabs((NewPlane-Plane)*LayerThickness);
    //    cout <<"delz "<< delz<<endl;

    Output[0]=StateVector[0] + StateVector[2]*delz;
    Output[1]=StateVector[1] + StateVector[3]*delz;
    Output[2]=StateVector[2];
    Output[3]=StateVector[3];
    Output[4]=StateVector[4];
    Output[5]=StateVector[5];

    done=true;
  }
  //cout <<" Input S1 "<< StateVector[0]<<" "<<StateVector[1]<<" "<<StateVector[2]<<" "<<StateVector[3]<<" "<<StateVector[4]<<endl;
  //cout <<" Output S1 "<< Output[0]<<" "<<Output[1]<<" "<<Output[2]<<" "<<Output[3]<<" "<<Output[4]<<endl;
  delete myswimmer;
  //  delete bf;
  return done;
}

void InoNewFitAlg::GetNoiseMatrix(const int Plane, const int NewPlane) {

  FieldPropagator* pFieldMap = FieldPropagator::FdPointer;
  double B[6];
  double x[3];
  double Bx,By;
  int DeltaPlane=abs(NewPlane-Plane);
  double DeltaZ=fabs(DeltaPlane*LayerThickness);
  x[0]= (posin[2*NewPlane]+posin[2*Plane])*0.5;
  x[1]= (posin[2*NewPlane+1]+posin[2*Plane+1])*0.5;
  x[2]= (posinz[NewPlane]+posin[Plane])*0.5;
  int cplane = (NewPlane>Plane)? Plane+1 : Plane-1;
 pFieldMap->ElectroMagneticField(x, Bx,By,cplane);
  B[0] =Bx*1000;B[1]= By*1000; B[2]=0;
  // This method is essentially the same as in SR fitter
  //   cout <<" InoTrackFitAlg : GetNoiseMatrix" << endl;

  for (int p=0; p<2; ++p) {
    for (int q=0; q<2; ++q) {
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
    double dzscatter = 0.5*fabs(posinz[NewPlane]-posin[Plane]);
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

    // double buvz[3] ={0.3*B[0], 0.3*B[1], 0.3*B[2]}; // GMA BX=By=1Tesla, Bz=0 and 0.3 factor comes from
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
    Q_k_minus[1][0]=dzscatter2*Sigma34Squared;
    Q_k_minus[1][1]=dzscatter2*Sigma44Squared;
    
    // Display
    if(debug_fit) {
      cout<<"------------------------------------------------------------"<<endl;
      cout << "1e6 * Q_k_minus" << endl;
      for(int i=0; i<2; ++i) {
	for(int j=0; j<2; ++j) {
	  cout << 1e6*Q_k_minus[i][j] << " ";
	}
	cout << endl;
      }
      cout<<"------------------------------------------------------------"<<endl;
    }

  }

}
