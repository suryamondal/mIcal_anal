#ifndef Fcnsg_h
#define Fcnsg_h

double xvalin[200];
double yvalin[200];
double errxy2;
int nsize1;
void fcnsg(Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t flag) {
  double fval=0;

  for (int ij=0; ij<nsize1; ij++) {
    fval += abs((xvalin[ij]-par[0])*(xvalin[ij]-par[0]) 
	       +(yvalin[ij]-par[1])*(yvalin[ij]-par[1])
		- par[2]*par[2]); //use proper cov matrix for X/Y
    //    cout<<" fval = "<<ij<<" "<<fit_input3v[ij]<<" "<<fval<<endl;
  }
  f = fval/errxy2; //Use some sort of uncorrelated error
  //  cout<<"par "<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<f<<endl;
  //but use covariant error matrix of X-Y, which is also depend on layer number
  //  cout<<" fvalxxxxxxxxxxxxx = "<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<f<<endl;
}

double sfdph(double phi1, double phi2) {
  double dphi = phi2 - phi1;
  //  cout <<"dphi "<<dphi<<endl;
  if (dphi > M_PI) { dphi = 2*M_PI - dphi;}
  if (dphi < -M_PI){ dphi = 2*M_PI + dphi;}
  //  cout <<"dphi2 "<<dphi<<endl;
  return dphi;
}

#endif
