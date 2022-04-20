#include "MultiSimAnalysisDigi.hh"

MultiSimAnalysisDigi *MultiSimAnalysisDigi::AnPointer;

MultiSimAnalysisDigi::MultiSimAnalysisDigi() {
  AnPointer = this;
  paradef = DetectorParameterDef::AnPointer;
  CardFile = ParameterMessenger::AnPointer;
  DetectorType = CardFile->GetDetectorType();

  numberInLA = paradef->GetnLayer();
  isInOut = CardFile->GetInputOutput();
  isVisOut = CardFile->GetVisualOutput();
  isXtermOut = CardFile->GetXTermOutput();
  collatedIn = CardFile->GetCollatedIn();
  nmxhit	= 10000;
  ievt		= 0;

  ievt_wt=0;
  intxn_id=0;
  
  CorrTimeError = CardFile->GetCorrTimeSmr();
  UnCorrTimeError = CardFile->GetUnCorrTimeSmr();
  TimeToDigiConv = CardFile->GetTimeToDigiConv();
  // cout<<"Hellllllllllll"<<TimeToDigiConv<<endl;
  SignalSpeed = CardFile->GetSignalSpeed();

  H	=0;
  Hp	=0;
  
  pRootFile	=0;
  inputRootFile=0;
  pVisFile = 0;
  pEventTree	=0;
  inputEventTree=0;
  visTree = 0;
}

MultiSimAnalysisDigi::~MultiSimAnalysisDigi() {
  cout<<"Deleting MultiSimAnalysisDigi class pointer ..."<<endl;
}

void MultiSimAnalysisDigi::OpenInputRootFiles(char* infile) {
  
  char inRootFile[400];
  
  if(isInOut==0) {
    sprintf(inRootFile,"%s.root",infile);
    inputRootFile = new TFile(inRootFile, "read");
    cout<< "Data is read from simulation output rootfile : "<< inRootFile <<endl;
    
    inputEventTree = (TTree*)inputRootFile->Get("T3");
    inputEventTree->SetBranchAddress("irun",&irun);
    inputEventTree->SetBranchAddress("ievt",&ievt);
    inputEventTree->SetBranchAddress("ngent",&ngent);
    inputEventTree->SetBranchAddress("pidin",pidin);
    inputEventTree->SetBranchAddress("momin",momin);
    inputEventTree->SetBranchAddress("thein",thein);
    inputEventTree->SetBranchAddress("phiin",phiin);
    inputEventTree->SetBranchAddress("posxin",posxin);
    inputEventTree->SetBranchAddress("posyin",posyin);
    inputEventTree->SetBranchAddress("poszin",poszin); 
    inputEventTree->SetBranchAddress("nsimht", &nsimht);
    inputEventTree->SetBranchAddress("detid", detid);
    inputEventTree->SetBranchAddress("simpdgid", simpdgid);
    inputEventTree->SetBranchAddress("simtime", simtime);
    inputEventTree->SetBranchAddress("simenr", simenr);
    inputEventTree->SetBranchAddress("simvx", simvx);
    inputEventTree->SetBranchAddress("simvy", simvy);
    inputEventTree->SetBranchAddress("simvz", simvz);
    inputEventTree->SetBranchAddress("simpx", simpx);
    inputEventTree->SetBranchAddress("simpy", simpy);
    inputEventTree->SetBranchAddress("simpz", simpz);
    inputEventTree->SetBranchAddress("simlocvx", simlocvx);
    inputEventTree->SetBranchAddress("simlocvy", simlocvy);
    inputEventTree->SetBranchAddress("ngenerated", &ngenerated);
    inputEventTree->SetBranchAddress("naperture", &naperture);
  } else if(isInOut==1) {
    sprintf(inRootFile,"%s.root",infile);
    inputRootFile = new TFile(inRootFile, "read");
    cout<< "Data is read from digitization output file : "<< inRootFile  <<endl;
   
    inputEventTree = (TTree*)inputRootFile->Get("T2");
    inputEventTree->SetBranchAddress("irun",&irun);
    inputEventTree->SetBranchAddress("ievt",&ievt);
    inputEventTree->SetBranchAddress("ngent",&ngent);
    inputEventTree->SetBranchAddress("pidin",pidin);
    inputEventTree->SetBranchAddress("momin",momin);
    inputEventTree->SetBranchAddress("thein",thein);
    inputEventTree->SetBranchAddress("phiin",phiin);
    inputEventTree->SetBranchAddress("posxin",posxin);
    inputEventTree->SetBranchAddress("posyin",posyin);
    inputEventTree->SetBranchAddress("poszin",poszin); 
    inputEventTree->SetBranchAddress("ndigiht", &ndigiht);
    inputEventTree->SetBranchAddress("trigx", &trigx);
    inputEventTree->SetBranchAddress("trigy", &trigy);
    inputEventTree->SetBranchAddress("ngenerated", &ngenerated);
    inputEventTree->SetBranchAddress("naperture", &naperture);
    inputEventTree->SetBranchAddress("triggeracceptance", &triggeracceptance);
    inputEventTree->SetBranchAddress("stripid", stripid);
    inputEventTree->SetBranchAddress("digipdgid", digipdgid);
    inputEventTree->SetBranchAddress("digitime", digitime);
    inputEventTree->SetBranchAddress("digitruetime", digitruetime);
    inputEventTree->SetBranchAddress("digienr", digienr);
    inputEventTree->SetBranchAddress("digivx", digivx);
    inputEventTree->SetBranchAddress("digivy", digivy);
    inputEventTree->SetBranchAddress("digivz", digivz);
    inputEventTree->SetBranchAddress("digipx", digipx);
    inputEventTree->SetBranchAddress("digipy", digipy);
    inputEventTree->SetBranchAddress("digipz", digipz);
    // inputEventTree->SetBranchAddress("diginoise", diginoise);
  } else if(isInOut==2) { // datareading
    sprintf(inRootFile,"%s",infile);
    inputRootFile = new TFile(inRootFile, "read");
    cout<< "Data is read from data file : "<< inRootFile  <<endl;
    
    inputEventTree = (TTree*)inputRootFile->Get("evetree");
    data_event = new RPCEve(inputEventTree);
    // data_event = new evetree(inputEventTree);
    data_event->Loop();
    
  }
    
}

void MultiSimAnalysisDigi::OpenOutputRootFiles(char* outfile) {

  char outRootFile[400];
  sprintf(outRootFile,"%s.root",outfile);
  if(isInOut==0) {
    pRootFile = new TFile(outRootFile, "RECREATE");
    if (!pRootFile) {
      cout << "Error opening output root file !" << endl;
      exit(-1);
    } else {
      cout<< "Output stored in root file: "<< outRootFile <<endl;
    }
    pEventTree = new TTree("T2", "INODIGI");
    pEventTree->Branch("irun",&irun,"irun/i"); //VALGRIND
    pEventTree->Branch("ievt",&ievt2,"ievt/i");
    pEventTree->Branch("ngent",&ngent,"ngent/i");
    pEventTree->Branch("pidin",pidin,"pidin[ngent]/I");
    pEventTree->Branch("ievt_wt",&ievt_wt,"ievt_wt/F");
    pEventTree->Branch("intxn_id",&intxn_id,"intxn_id/I");
    pEventTree->Branch("momin",momin,"momin[ngent]/F");
    pEventTree->Branch("thein",thein,"thein[ngent]/F");
    pEventTree->Branch("phiin",phiin,"phiin[ngent]/F");
    pEventTree->Branch("posxin",posxin,"posxin[ngent]/F");
    pEventTree->Branch("posyin",posyin,"posyin[ngent]/F");
    pEventTree->Branch("poszin",poszin,"poszin[ngent]/F");
    pEventTree->Branch("ndigiht", &ndigiht, "ndigiht/i");
    pEventTree->Branch("trigx", &trigx, "trigx/I");
    pEventTree->Branch("trigy", &trigy, "trigy/I");
    pEventTree->Branch("ngenerated",&ngenerated,"ngenerated/i");
    pEventTree->Branch("naperture",&naperture,"naperture/i");
    pEventTree->Branch("triggeracceptance",&triggeracceptance,"triggeracceptance/i");
    pEventTree->Branch("stripid", stripid, "stripid[ndigiht]/i");
    pEventTree->Branch("digipdgid", digipdgid, "digipdgid[ndigiht]/I");
    pEventTree->Branch("digitime", digitime, "digitime[ndigiht]/I");
    pEventTree->Branch("digitruetime", digitruetime, "digitruetime[ndigiht]/I");
    pEventTree->Branch("digienr", digienr, "digienr[ndigiht]/F");
    pEventTree->Branch("digivx", digivx, "digivx[ndigiht]/F"); 
    pEventTree->Branch("digivy", digivy, "digivy[ndigiht]/F");
    pEventTree->Branch("digivz", digivz, "digivz[ndigiht]/F");
    pEventTree->Branch("digipx", digipx, "digipx[ndigiht]/F");
    pEventTree->Branch("digipy", digipy, "digipy[ndigiht]/F");
    pEventTree->Branch("digipz", digipz, "digipz[ndigiht]/F");
    
    
    if (!pEventTree) {
      cout << "Error allocating Tree !" << endl;
      exit(-1);
    } else {
      cout << "Output stored in root tree: T2; INODIGI" << endl;
    }
  } else if(isInOut<3) {
    pRootFile = new TFile(outRootFile, "RECREATE");
    if (!pRootFile) {
      cout << "Error opening output root file !" << endl;
      exit(-1);
    } else {
      cout<< "Output stored in root file: "<< outRootFile <<endl;
    }

    if(isInOut==1) {
      pEventTree = new TTree("T1", "INORECO");
    } else if(isInOut==2) {
      pEventTree = new TTree("T4", "INODATA");
    }
    pEventTree->Branch("irun",&irun,"irun/i"); //VALGRIND
    pEventTree->Branch("ievt",&ievt2,"ievt/i");
    pEventTree->Branch("ngent",&ngent,"ngent/i");
    pEventTree->Branch("pidin",pidin,"pidin[ngent]/I");
    pEventTree->Branch("ievt_wt",&ievt_wt,"ievt_wt/F");
    pEventTree->Branch("intxn_id",&intxn_id,"intxn_id/I");
    pEventTree->Branch("momin",momin,"momin[ngent]/F");
    pEventTree->Branch("thein",thein,"thein[ngent]/F");
    pEventTree->Branch("phiin",phiin,"phiin[ngent]/F");
    pEventTree->Branch("posxin",posxin,"posxin[ngent]/F");
    pEventTree->Branch("posyin",posyin,"posyin[ngent]/F");
    pEventTree->Branch("poszin",poszin,"poszin[ngent]/F");
    pEventTree->Branch("ngenerated",&ngenerated,"ngenerated/i");
    pEventTree->Branch("naperture",&naperture,"naperture/i");
    pEventTree->Branch("triggeracceptance",&triggeracceptance,"triggeracceptance/i");
    pEventTree->Branch("hw_trig",&hw_trig,"hw_trig/I");
    // pEventTree->Branch("hw_trigy",&hw_trigy,"hw_trigy/I");
    pEventTree->Branch("sw_trigx",&sw_trigx,"sw_trigx/I");
    pEventTree->Branch("sw_trigy",&sw_trigy,"sw_trigy/I");
    pEventTree->Branch("ntrkt",&ntrkt,"ntrkt/i");
    pEventTree->Branch("itype",itype,"itype[ntrkt]/I");
    pEventTree->Branch("nLayer",&nLayer,"nLayer/I");
    pEventTree->Branch("nhits",nhits,"nhits[ntrkt]/I");
    pEventTree->Branch("nhits_finder",nhits_finder,"nhits_finder[ntrkt]/I");
    pEventTree->Branch("chisq",chisq,"chisq[ntrkt]/F");
    pEventTree->Branch("cvalue",cvalue,"cvalue[ntrkt]/F");
    pEventTree->Branch("trkmm",trkmm,"trkmm[ntrkt]/F");
    pEventTree->Branch("trkth",trkth,"trkth[ntrkt]/F");
    pEventTree->Branch("trkph",trkph,"trkph[ntrkt]/F");
    pEventTree->Branch("momvx",momvx,"momvx[ntrkt]/F");
    pEventTree->Branch("thevx",thevx,"thevx[ntrkt]/F");
    pEventTree->Branch("phivx",phivx,"phivx[ntrkt]/F");
    pEventTree->Branch("posxvx",posxvx,"posxvx[ntrkt]/F");
    pEventTree->Branch("posyvx",posyvx,"posyvx[ntrkt]/F");
    pEventTree->Branch("poszvx",poszvx,"poszvx[ntrkt]/F");
    pEventTree->Branch("momend",momend,"momend[ntrkt]/F");
    pEventTree->Branch("theend",theend,"theend[ntrkt]/F");
    pEventTree->Branch("phiend",phiend,"phiend[ntrkt]/F");
    pEventTree->Branch("posxend",posxend,"posxend[ntrkt]/F");
    pEventTree->Branch("posyend",posyend,"posyend[ntrkt]/F");
    pEventTree->Branch("strpxend",strpxend,"strpxend[ntrkt]/I");
    pEventTree->Branch("strpyend",strpyend,"strpyend[ntrkt]/I");
    pEventTree->Branch("poszend",poszend,"poszend[ntrkt]/F");
    pEventTree->Branch("momds"   ,momds   ,"momds[ntrkt]/F");
    pEventTree->Branch("momrg"   ,momrg   ,"momrg[ntrkt]/F");
    pEventTree->Branch("vtxzplane",vtxzplane,"vtxzplane[ntrkt]/I");
    pEventTree->Branch("endzplane",endzplane,"endzplane[ntrkt]/I");
    pEventTree->Branch("nvisht", &nvisht, "nvisht/i");
    pEventTree->Branch("fitposxx",fitposxx,"fitposxx[nvisht]/F");
    pEventTree->Branch("fitposyy",fitposyy,"fitposyy[nvisht]/F");
    pEventTree->Branch("fitposzz",fitposzz,"fitposzz[nvisht]/F");
    pEventTree->Branch("fitlayzz",fitlayzz,"fitlayzz[nvisht]/F");
    // pEventTree->Branch("fitlayx2",fitlayx2,"fitlayx2[nvisht]/F");
    // pEventTree->Branch("fitlayx3",fitlayx3,"fitlayx3[nvisht]/F");
    // pEventTree->Branch("fitlayx4",fitlayx4,"fitlayx4[nvisht]/F");
    pEventTree->Branch("fitlaymom",fitlaymom,"fitlaymom[nvisht]/F");
    pEventTree->Branch("extrapolxx",extrapolxx,"extrapolxx[nvisht]/F");
    pEventTree->Branch("extrapolyy",extrapolyy,"extrapolyy[nvisht]/F");
    pEventTree->Branch("extrapolmom",extrapolmom,"extrapolmom[nvisht]/F");
    pEventTree->Branch("momdiff1",&momdiff1,"momdiff1/F");
    pEventTree->Branch("radialdiff1",&radialdiff1,"radialdiff1/F");
    // pEventTree->Branch("fitlaythe",fitlaythe,"fitlaythe[nvisht]/F");
    // pEventTree->Branch("fitlayphi",fitlayphi,"fitlayphi[nvisht]/F");
    pEventTree->Branch("nvisclst", &nvisclst, "nvisclst/i");
    pEventTree->Branch("clstposxx",clstposxx,"clstposxx[nvisclst]/F");
    pEventTree->Branch("clstposyy",clstposyy,"clstposyy[nvisclst]/F");
    pEventTree->Branch("clstposzz",clstposzz,"clstposzz[nvisclst]/F");
    pEventTree->Branch("clstposzpln",clstposzpln,"clstposzpln[nvisclst]/I");
    
    pEventTree->Branch("strtnhitsx",&strtnhitsx,"strtnhitsx/I");
    pEventTree->Branch("strtchisqx",&strtchisqx,"strtchisqx/F");
    pEventTree->Branch("strtintercptx",&strtintercptx,"strtintercptx/F");
    pEventTree->Branch("strtslopex",&strtslopex,"strtslopex/F");
    pEventTree->Branch("strtnhitsy",&strtnhitsy,"strtnhitsy/I");
    pEventTree->Branch("strtchisqy",&strtchisqy,"strtchisqy/F");
    pEventTree->Branch("strtintercpty",&strtintercpty,"strtintercpty/F");
    pEventTree->Branch("strtslopey",&strtslopey,"strtslopey/F");
    pEventTree->Branch("strtxexpec",&strtxexpec,"strtxexpec/I");
    pEventTree->Branch("strtyexpec",&strtyexpec,"strtyexpec/I");

    pEventTree->Branch("simpleradii",&simpleradii,"simpleradii/F");
    pEventTree->Branch("simplecurv",&simplecurv,"simplecurv/F");
    pEventTree->Branch("simplex0",&simplex0,"simplex0/F");
    pEventTree->Branch("simplez0",&simplez0,"simplez0/F");
    pEventTree->Branch("simplechisqpos",&simplechisqpos,"simplechisqpos/F");
    pEventTree->Branch("simplechisqneg",&simplechisqneg,"simplechisqneg/F");
    pEventTree->Branch("simplechisqcndn",&simplechisqcndn,"simplechisqcndn/F");
    pEventTree->Branch("simpleavgxpos",&simpleavgxpos,"simpleavgxpos/F");
    pEventTree->Branch("simpleavgxneg",&simpleavgxneg,"simpleavgxneg/F");
    pEventTree->Branch("simpleavgxcndn",&simpleavgxcndn,"simpleavgxcndn/F");
    pEventTree->Branch("simpleavgxmeas",&simpleavgxmeas,"simpleavgxmeas/F");
    pEventTree->Branch("simplenhits",&simplenhits,"simplenhits/I");
    pEventTree->Branch("simplexexpec",&simplexexpec,"simplexexpec/I");
    pEventTree->Branch("nhits_last",&nhits_last,"nhits_last/I");
    pEventTree->Branch("nhits_last_m1",&nhits_last_m1,"nhits_last_m1/I");
    pEventTree->Branch("nhits_below",&nhits_below,"nhits_below/I");
    pEventTree->Branch("ftime_last",&ftime_last,"ftime_last/F");
    
    if(isInOut==2) {
      pEventTree->Branch("ntdc1x",&ntdc1x,"ntdc1x/i");
      pEventTree->Branch("tdcID1x",tdcID1x,"tdcID1x[ntdc1x]/I");
      pEventTree->Branch("TDCval1x",TDCval1x,"TDCval1x[ntdc1x]/F");
      pEventTree->Branch("ntstrp1x",&ntstrp1x,"ntstrp1x/i");
      pEventTree->Branch("StrpID1x",StrpID1x,"StrpID1x[ntstrp1x]/I");
      pEventTree->Branch("ntdc2x",&ntdc2x,"ntdc2x/i");
      pEventTree->Branch("tdcID2x",tdcID2x,"tdcID2x[ntdc2x]/I");
      pEventTree->Branch("TDCval2x",TDCval2x,"TDCval2x[ntdc2x]/F");
      pEventTree->Branch("ntstrp2x",&ntstrp2x,"ntstrp2x/i");
      pEventTree->Branch("StrpID2x",StrpID2x,"StrpID2x[ntstrp2x]/I");
      
      pEventTree->Branch("ntdc1y",&ntdc1y,"ntdc1y/i");
      pEventTree->Branch("tdcID1y",tdcID1y,"tdcID1y[ntdc1y]/I");
      pEventTree->Branch("TDCval1y",TDCval1y,"TDCval1y[ntdc1y]/F");
      pEventTree->Branch("ntstrp1y",&ntstrp1y,"ntstrp1y/i");
      pEventTree->Branch("StrpID1y",StrpID1y,"StrpID1y[ntstrp1y]/I");
      pEventTree->Branch("ntdc2y",&ntdc2y,"ntdc2y/i");
      pEventTree->Branch("tdcID2y",tdcID2y,"tdcID2y[ntdc2y]/I");
      pEventTree->Branch("TDCval2y",TDCval2y,"TDCval2y[ntdc2y]/F");
      pEventTree->Branch("ntstrp2y",&ntstrp2y,"ntstrp2y/i");
      pEventTree->Branch("StrpID2y",StrpID2y,"StrpID2y[ntstrp2y]/I");
      pEventTree->Branch("ntrecord1x",&ntrecord1x,"ntrecord1x/i");
      pEventTree->Branch("striprec1x",striprec1x,"striprec1x[ntrecord1x]/I");
      pEventTree->Branch("tdcrec1x",tdcrec1x,"tdcrec1x[ntrecord1x]/F");
      pEventTree->Branch("ntrecord1y",&ntrecord1y,"ntrecord1y/i");
      pEventTree->Branch("striprec1y",striprec1y,"striprec1y[ntrecord1y]/I");
      pEventTree->Branch("tdcrec1y",tdcrec1y,"tdcrec1y[ntrecord1y]/F");
      pEventTree->Branch("ntrecord2x",&ntrecord2x,"ntrecord2x/i");
      pEventTree->Branch("striprec2x",striprec2x,"striprec2x[ntrecord2x]/I");
      pEventTree->Branch("tdcrec2x",tdcrec2x,"tdcrec2x[ntrecord2x]/F");
      pEventTree->Branch("ntrecord2y",&ntrecord2y,"ntrecord2y/i");
      pEventTree->Branch("striprec2y",striprec2y,"striprec2y[ntrecord2y]/I");
      pEventTree->Branch("tdcrec2y",tdcrec2y,"tdcrec2y[ntrecord2y]/F");
    }
    
    if (!pEventTree) {
      cout << "Error allocating Tree !" << endl;
      exit(-1);
    } else {
      if(isInOut==2) {cout << "Output stored in root tree: T4; INODATA" << endl;}
      else {cout << "Output stored in root tree: T1; INORECO" << endl;}
    }
    
    if(isVisOut==1) {
      char outVisFile[300];
      sprintf(outVisFile,"%s.inh",outfile);
      pVisFile = new TFile(outVisFile, "RECREATE"); //VALGRIND
      if (!pVisFile) {
	cout << "Error opening .inh file !" << endl;
	exit(-1);
      } else {
	cout << "Vis Hits stored in: " << outVisFile << endl;
      }
      pVisFile->cd();
      EveCnt=0;//Event Counter
      if(!H) H = new Hits(); //VALGRIND
      if(!Hp) Hp= new HitPos();
      nloops=0;// number of tree fills
      if(!visTree){visTree = new TTree("Hitstree","Geant Hits File");}
      visTree->Branch("Hits_Branch","Hits",&H,1600000,2);
      //H->Clear();
      //Hp->Clear();
      //H->ClearTracks();
      H->ENum=0;
    }

    char namex[200];
    for(int iaj=0; iaj<numberInLA; iaj++) {
      sprintf(namex,"hdifftime1_xy_%i",iaj);
      hdifftime1[iaj] = new TH1D(namex,namex,120,-25.,25.);
      sprintf(namex,"hdifftime2_xy_%i",iaj);
      hdifftime2[iaj] = new TH1D(namex,namex,120,-5.,5.);
      
      sprintf(namex,"hxtime_ext_%i",iaj);
      hxtime_ext[iaj] = new TH1D(namex,namex,120,-25.,25.);
      sprintf(namex,"hytime_ext_%i",iaj);
      hytime_ext[iaj] = new TH1D(namex,namex,120,-25.,25.);
      
      sprintf(namex,"hxpos_ext_%i",iaj);
      hxpos_ext[iaj] = new TH1D(namex,namex,120,-6.,6.);
      sprintf(namex,"hypos_ext_%i",iaj);
      hypos_ext[iaj] = new TH1D(namex,namex,120,-6.,6.);
      sprintf(namex,"hxpos_ext_kalman_%i",iaj);
      hxpos_ext_kalman[iaj] = new TH1D(namex,namex,120,-6.,6.);
      sprintf(namex,"hypos_ext_kalman_%i",iaj);
      hypos_ext_kalman[iaj] = new TH1D(namex,namex,120,-6.,6.);
      
      sprintf(namex,"h_hit_time_ext_%i",iaj);
      h_hit_time_ext[iaj] = new TH1D(namex,namex,120,-25.,25.);
      // for(int jak=0; jak<8; jak++) {
      // 	sprintf(namex,"xtdc_minus_ref_l%i_%i",iaj,jak);
      // 	xtdc_minus_ref[iaj][jak] = new TH1D(namex,namex,120,-20000.,200000);
      // 	sprintf(namex,"ytdc_minus_ref_l%i_%i",iaj,jak);
      // 	ytdc_minus_ref[iaj][jak] = new TH1D(namex,namex,120,-20000.,200000);
      // 	sprintf(namex,"tshift_xtdc_minus_ref%i_%i",iaj,jak);
      // 	tshift_xtdc_minus_ref[iaj][jak] = new TH1D(namex,namex,120,-20000.,200000);
      // 	sprintf(namex,"tshift_ytdc_minus_ref%i_%i",iaj,jak);
      // 	tshift_ytdc_minus_ref[iaj][jak] = new TH1D(namex,namex,120,-20000.,200000);
      // }
    }
    DGap = new TH1D("DetGap"," Plot to see the performance of code at detgap ",10,0,10);  //asm
    ShwXw	= new TH1D("ShwXw ","This is a distribution for X position resolution in shower region " ,80 , -3.14 , 3.14);  //asm
    ShwYw	= new TH1D("ShwYw ","This is a distribution for Y position resolution in shower region " ,80 , -3.14 , 3.14);  //asm
    trk_edge = new TH2D("theta_dmn at trk_edge","theta_dmn at trk_edge" ,20, 0 ,20,200,0,3);  //asm
    trk_gap = new TH1D("difference in energy lost in the  dead space and the tracklength - )"," dpp-m*dxx + c" ,100, -5 ,5);  //asm
    if(!pPosX)pPosX = new TH1D("deltaX", "#Delta X (cm)", 100, -52.5, 52.5);
    if(!pPosY)pPosY = new TH1D("deltaY", "#Delta Y (cm)", 100, -52.5, 52.5);
    if(!pPosZ)pPosZ = new TH1D("deltaZ", "#Delta Z (cm)", 100, -2.5, 2.5);
    if(!pPosXX)pPosXX = new TH2D("deltaXX", "#Delta XX (cm)", 100, -2500, 2500, 100, -52.5, 52.5);
    if(!pPosYY)pPosYY = new TH2D("deltaYY", "#Delta YY (cm)", 100, -500, 500, 100, -52.5, 52.5);
    if(!pPosZZ)pPosZZ = new TH2D("deltaZZ", "#Delta ZZ (cm)", 100, -1000, 1000, 100, -2.5, 2.5);
  } 
}

void MultiSimAnalysisDigi::OpenCollatedRootFile() {
  cout<<"Reading Collated Root Input file for digitization."<<endl;
  cout<<"collatedin "<<collatedIn<<endl;
  if(collatedIn) {
    // collatedRootFile = new TFile("Collated_trig6789_20180614_mical.root","read");
    collatedRootFile = new TFile("Collated_mIcal_trig6789_20201009.root","read");
    if(!collatedRootFile) {
      cout << "Error opening collated file !" << endl;
      exit(-1);
    } else{
      cout<<"Collated file open successful..."<<endl;
      char namex[200];
      for(int iki=0; iki<numberInLA; iki++) {
	sprintf(namex,"inefficiency_corx_l%i",iki);
	inefficiency_corx[iki] = (TH2D*)collatedRootFile->Get(namex);
	sprintf(namex,"inefficiency_uncx_l%i",iki);
	inefficiency_uncx[iki] = (TH2D*)collatedRootFile->Get(namex);
	// cout<<"inefficiency_uncx["<<iki<<"] "<<inefficiency_uncx[iki]<<endl;
	sprintf(namex,"inefficiency_uncy_l%i",iki);
	inefficiency_uncy[iki] = (TH2D*)collatedRootFile->Get(namex);
	// cout<<"inefficiency_uncy["<<iki<<"] "<<inefficiency_uncy[iki]<<endl;
	sprintf(namex,"triggereffi_xevt_l%i",iki);
	triggereffi_xevt[iki] = (TH2D*)collatedRootFile->Get(namex);
	sprintf(namex,"triggereffi_yevt_l%i",iki);
	triggereffi_yevt[iki] = (TH2D*)collatedRootFile->Get(namex);
	sprintf(namex,"strp_xmulsim_cor_l%i",iki);
	strp_xmulsim_cor[iki] = (TH2D*)collatedRootFile->Get(namex);
	sprintf(namex,"strp_ymulsim_cor_l%i",iki);
	strp_ymulsim_cor[iki] = (TH2D*)collatedRootFile->Get(namex);
      }
      for(int iki=0; iki<numberInLA; iki++) {
	for(int ikj=0; ikj<16; ikj++) {
	  for(int ikk=0; ikk<16; ikk++) {
	    sprintf(namex,"strp_xmullaysim_l%i_%i_%i",iki,ikj,ikk);
	    block_xmulsim[iki][ikj][ikk] = (TH2D*)collatedRootFile->Get(namex);
	    sprintf(namex,"strp_ymullaysim_l%i_%i_%i",iki,ikj,ikk);
	    block_ymulsim[iki][ikj][ikk] = (TH2D*)collatedRootFile->Get(namex);
	  }
	}
      } // for(int iki=0; iki<numberInLA; iki++) {
    }
  }
}

void MultiSimAnalysisDigi::CloseInputRootFiles() {
  if (inputRootFile) {
    inputRootFile->cd();
    inputRootFile->Close();
    // delete inputEventTree; inputEventTree=0;
    delete inputRootFile; inputRootFile=0;
    cout << "Input root file close !" << endl;
  } else {
    cout << "No input file to close.... !" << endl;
  }
  
}

void MultiSimAnalysisDigi::CloseOutputRootFiles() {
  if (pRootFile) {
    pRootFile->cd();
    if (pEventTree) { pEventTree->Write(); delete pEventTree; pEventTree=0;}
    pRootFile->Write();
    // for(int iai=0; iai<numberInLA; iai++) {
    //   if(hdifftime1[iai]) {hdifftime1[iai]->Write(); delete hdifftime1[iai]; hdifftime1[iai]=0;}
    //   if(hdifftime2[iai]) {hdifftime2[iai]->Write(); delete hdifftime2[iai]; hdifftime2[iai]=0;}
    //   if(hdifftime1[iai]) {hdifftime1[iai]->Write(); delete hdifftime1[iai]; hdifftime1[iai]=0;}
    pRootFile->Close();
    delete pRootFile; pRootFile=0;
    cout << "Output root file  close !" << endl;
  }  else {
    cout << "Output root file cannot close !" << endl;
  } // if (pRootFile) {

  if (pVisFile) {
    pVisFile->cd();
    visTree->Fill(); // fill tree
    H->Clear();
    Hp->Clear();
    H->ClearTracks();
    pVisFile->Write(); //VALGRIND
    cout<< "Write to .inh file  done"<<endl;
    if (visTree) {delete visTree; visTree=0;}
    pVisFile->Close();
    if (Hp) {delete Hp; Hp=0;}
    if (H) {delete H; H=0;}
    delete pVisFile; pVisFile=0;
    cout << "Visualization file  close !" << endl;
  } else {
    cout << "No output Hit Display Tree !" << endl;
  }

  if (collatedRootFile) {
    collatedRootFile->cd();
    for(int iki=0; iki<numberInLA; iki++) {
      delete inefficiency_corx[iki];
      delete inefficiency_uncx[iki];
      delete inefficiency_uncy[iki];
      delete triggereffi_xevt[iki];
      delete triggereffi_yevt[iki];
      delete strp_xmulsim_cor[iki];
      delete strp_ymulsim_cor[iki];
    }
    collatedRootFile->Close();
    delete collatedRootFile; collatedRootFile=0;
    cout<<"Collated root file closed."<<endl;
  } else {
    cout << "No collated histograms !" << endl;
  }
  
}

void MultiSimAnalysisDigi::SetCorrTimeError(G4double val) {
  cout<<"void MultiSimAnalysisDigi::SetCorrTimeError(G4double "<<val<<")"<<endl;
  CorrTimeError = val;
}

void MultiSimAnalysisDigi::SetUnCorrTimeError(G4double val) {
  cout<<"void MultiSimAnalysisDigi::SetUnCorrTimeError(G4double "<<val<<")"<<endl;
  UnCorrTimeError = val;
}

void MultiSimAnalysisDigi::SetTimeToDigiConvVal(G4double val) {
  cout<<"void MultiSimAnalysisDigi::SetTimeToDigiConvVal(G4double "<<val<<")"<<endl;
  TimeToDigiConv = val;
}

void MultiSimAnalysisDigi::SetSignalSpeedVal(G4double val) {
  cout<<"void MultiSimAnalysisDigi::SetSignalSpeedVal(G4double "<<val<<")"<<endl;
  SignalSpeed = val;
}

void MultiSimAnalysisDigi::SaveGenVisFile() {
  if(isVisOut==1) {
    for(unsigned ij=0;ij<ngent;ij++) {
      H->NParticles++;
      Hp=  H->AddHits(0,0);
      Hp->TrackType=-14;
      Hp->ParCode= pidin[ij];
      Hp->ZZ= (((numberInLA*(parirlay[2]+parlay[2])*cm/m-parlay[2])- poszin[ij]*cm/m))/((parirlay[2]+parlay[2])*2*(1/m));
      Hp->XX= posxin[ij]*cm/m;
      Hp->YY= posyin[ij]*cm/m;
      Hp->pmag= momin[ij];
      Hp->pt= thein[ij];
      Hp->pp= phiin[ij];
    } // for(unsigned ij=0;ij< ngent;ij++) {
  } // if(isVisOut==1) {
}
