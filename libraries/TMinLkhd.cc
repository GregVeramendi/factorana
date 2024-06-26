#include "TMinLkhd.hh"

#include "Riostream.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMinuit.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "Math/DistFunc.h"

// #ifndef KNITRO_H__
// #include  "knitro.h"
// #endif

//#include <io.h>   // For access().
#include <sys/types.h>  // For stat().
#include <sys/stat.h> 
#include <cmath>
#include <iostream>
#include <sstream>
#include <assert.h>

#include "IpIpoptApplication.hpp"
#include "minlkhd_nlp.hh"

// ClassImp(TMinLkhd)

using namespace std;

const long double PI = 3.141592653589793238;

int fexist(const char *filename );
int dexist(const char *dirname );

void calcgausshermitequadrature(int n, double * x, double * w);

void MinuitLkhdFcn(Int_t &npar, Double_t *gin, Double_t &f,
	       Double_t *par, Int_t iflag);

void LkhdFcn(Int_t &npar, Double_t *gin, Double_t &f,
		   Double_t *par, Int_t iflag, Double_t *hess = NULL);


TMinLkhd::TMinLkhd()
  : TNamed()
{
  ClearDataMembers();
  fac_corr =0; 
  nfac=0;
  ntyp=0;
  type_model = -1;
  nquad_points=0;
  fac_npoints.clear();
  adapt_int = 0;
  adapt_int_thresh = 0.5; 
  nparam = 0;
  nfreeparam = 0;
  stage = 0;
  cpurank = 0;
  max_cpu_limit = -1.0;
  stoch_deriv_frac = -1.0;   
  CalcHess = 1;
  HessFactorScores = 0;
  printlvl = 1;
  printmargeffect = 0;
  current_printgroup = 0;
  current_estimationgroup = 0;
  nprintgroup = 0;
  initializing = 0;
  predicting = 0;
  predictobs = 0;
  includePriorFactorScore = 0;
  indexvar = -1;
  weightvar = -1;
  nsubsamples = 0;
  subsample_percent = 1.0;
  nchsamples = 0;
  nbootsamples = 0;
  bootstrapstart = 0;
  HessStdErr = 1;
  sampleposterior = 0;
  simIncData = 0;
  sim_nobs = 100000;
  newflag = 0;
  EstSequentMix = 0;
  current_sample = -1;
  cpurank = 0;
  initBetaOLS = 0;
  initEstOutcomeLoadings = 0;
  initFixedLoadings = 0;
  initVariance = -1.0;
  loadingMultiplier = 0.001;
  bootstrapstart = 0;
}



TMinLkhd::TMinLkhd(const char *name, const char *title,  Int_t nfactors, Int_t ntypes, Int_t fcorr, Int_t fnmix, Int_t nquad, Int_t thisstage)
  : TNamed(name,title)
{

  //Initialize variables  
  nobs=0;
  fac_corr =fcorr;
  fac_nmix =fnmix;
  est_nmix = fac_nmix;
  nfac=nfactors;
  ntyp = ntypes-1;
  type_model = 9999;
  //  nquad_points=nquad;
  nquad_points=abs(nquad);
  fac_npoints.clear();
  adapt_int = (nquad<0);
  adapt_int_thresh = 0.5; 
  stage = thisstage;
  max_cpu_limit = -1.0;
  stoch_deriv_frac = -1.0;   
  ClearDataMembers();
  initializing = 0;
  predicting = 0;
  predictobs = 0;
  includePriorFactorScore = 0;
  indexvar = -1;
  weightvar = -1;
  nsubsamples = 0;
  subsample_percent = 1.0;
  nchsamples = 0;
  nbootsamples = 0;
  bootstrapstart = 0;
  CalcHess = 1;
  HessFactorScores = 0;
  HessStdErr = 1;
  sampleposterior = 0;
  simIncData = 0;
  simWithData = 0;
  sim_nobs = 100000;
  newflag = 0;
  EstSequentMix = 0;
  current_sample = -1;
  cpurank = 0;
  printlvl = 1;
  printmargeffect = 0;
  current_printgroup = 0;
  current_estimationgroup = 0;
  initBetaOLS = 0;
  initFixedLoadings = 0;
  initEstOutcomeLoadings = 0;
  initVariance = -1.0;
  loadingMultiplier = 0.001;

//   newskipobs = 0;
//   newbootstrap = 0;
//   newfixpar = 0;
//   newmodignore = 0;

  // Calculate Gaussian quadrature points and weights
  if (nquad_points>0) {
    
    //Here we are reserving nquad_points elements for each
    //set of constants even when fewer are needed, bc it's easier.
    x.resize(nquad_points*nquad_points);
    w.resize(nquad_points*nquad_points);
    
    //Set constants for one quadrature points (used for adaptive integration)
    x.at(0) = 0.0;
    w.at(0) = sqrt(PI);
    
    //Create arrays for HGQ function
    Double_t * getx = new Double_t[nquad_points];
    Double_t * getw = new Double_t[nquad_points];
    
    Double_t sqrt2 = sqrt(2.0);
    for (UInt_t iquad = 2 ; iquad <= nquad_points; iquad++ ) {
      //      printf("*****Getting constants for %2d GH quadrature points.\n",iquad);
      calcgausshermitequadrature(iquad,getx,getw);
      for (UInt_t ipoint = 0 ; ipoint < iquad; ipoint++) {
	x.at((iquad-1)*nquad_points + ipoint) = sqrt2*getx[ipoint];
	w.at((iquad-1)*nquad_points + ipoint) = getw[ipoint];
	//	printf("quadrature points [%2d]= %22.16g %22.16g\n",iquad,GetHGQx(iquad, ipoint),GetHGQw(iquad, ipoint));
      }
      //      printf("\n\n");
    }

    delete [] getx;
    delete [] getw;
  }

  // initialize the parameters that describe the factors
  if (nfac>0) {

    //    fac_npoints.resize(nfac,nquad_points);
    
    if ((fac_corr!=0)&&(nfac!=2)) {
      cout << "***Error in TMinLkhd::TMinLkhd*** Can only estimate correlated factors for nfac=2!!" << endl;
      assert(0);
    }
    norm_models.resize(nfac,-2);

    // nmix = 1,2,3 implemented (nmix=1 --> no mixture)
    // fac_nmix var-covar matrices, (nmix-1)*nfac means, nmix-1 weights 
    if ( (fac_nmix==1) || (fac_nmix==2) || (fac_nmix==3)) {
      // n variance, n-1 mean, and n-1 weight parameters
      f_nvariance = nfac;
      if (fac_corr!=0) f_nvariance = (nfac*(nfac+1)/2);
      nfac_param = f_nvariance*fac_nmix + (fac_nmix-1)*nfac + (fac_nmix-1);
      nparam = nfac_param;
      parconstrained.assign(nfac_param,-1);
    }
    else {
      cout << "***Error in TMinLkhd::TMinLkhd*** Can only do 2 or 3 mixtures at this point!!" << endl;
      assert(0);
    }
  }
}

void TMinLkhd::ClearDataMembers() {
  if (nobs!=0) {
     nobs = 0;
     //     cout << "deleting data" << endl;
     data.clear();
     var_table.clear();
     //     delete data;
     //     delete var_table;
  }

  //  cout << "clearing vectors" << endl;
  models.clear();
  param.clear();
  param_err.clear();
  param_fixval.clear();
  parfixed.clear();
  parconstrained.clear();
  skipobs.clear();
  bootstrapobs.clear();
}

TMinLkhd::~TMinLkhd() {
  if (nquad_points>0) {
    x.clear();
    w.clear();
//    cout << "deleting quadrature points" << endl;
//     delete x;
//     delete w;
  }
  ClearDataMembers();
}


void TMinLkhd::SetWorkingDir(TString thisdir) {
  if (dexist(thisdir.Data())==0) {
    if (!thisdir.EndsWith("/")) thisdir.Append("/");
    workingdir = thisdir;
    if (cpurank==0) printf("Working Directory set to:%s\n\n",workingdir.Data());
  }
  else if (cpurank==0) {
    printf("Not a valid working directory!!\n\n");
    assert(0);
  }
}


void TMinLkhd::SubSample(UInt_t nsamples, Double_t keep_perc) {
  nsubsamples = nsamples; 
  subsample_percent = keep_perc;

  skipobs.resize(nobs,0);
  
  TString filename;
  Int_t check = 1;
  for (UInt_t isample = 0 ; isample < nsubsamples; isample++) {
    
    std::stringstream out;
    out << isample;
    filename = TString("subsample_").Append(out.str()).Append(".txt");
    filename.Prepend(workingdir);
    
    if (!fexist((char *)filename.Data())) check = 0;
    
  }

  if (check==0) {
    TRandom *r3 = new TRandom3();
    for (UInt_t isample = 0 ; isample < nsubsamples ; isample++) {
      vector<Int_t> obsarray;
      obsarray.resize(nobs,0);
      
      Int_t count = 0;
      UInt_t firstobs = 9999;
      while (count < (1.0-subsample_percent)*nobs) {
	UInt_t ranobs = UInt_t(r3->Rndm()*(nobs));
	if (ranobs==nobs) ranobs = 0;
	if (obsarray[ranobs]==0) {
	  obsarray[ranobs] = 1;
	  count++;
	  if (ranobs<firstobs) firstobs = ranobs;
	}
      }

      Printf("Generating subsample file #%d, first skipped obs=%d\n",isample,firstobs);
      std::stringstream out;
      out << isample;
      filename = TString("subsample_").Append(out.str()).Append(".txt");
      filename.Prepend(workingdir);

      FILE * pFile;
      pFile = fopen (filename.Data(),"w");
      for (UInt_t iobs=0 ; iobs< nobs ; iobs++) {
	fprintf (pFile, "%5d ",obsarray[iobs]);	  
      }
      fprintf (pFile, "\n");
      fclose (pFile);
    }
  }

}

void TMinLkhd::BootStrap(UInt_t nsamples, UInt_t gensample) {
  nbootsamples = nsamples; 
  printlvl=0;
  
  if ((gensample!=0)&&(cpurank==0)) {

    if (nobs==0) {
      cout << "***Error in TMinLkhd::BootStrap*** You must add the data file before you can generate bootstrap samples" << endl;
      assert(0);
    }

    TString filename;
    TRandom *r3 = new TRandom3();
    for (UInt_t isample = 0 ; isample < nbootsamples ; isample++) {
      vector<Int_t> obsarray;
      obsarray.resize(nobs,0);
      
      for (UInt_t iobs = 0 ; iobs < nobs; iobs++) {
	UInt_t ranobs = UInt_t(r3->Rndm()*(nobs));
	if (ranobs==nobs) ranobs=0;
	//	obsarray[iobs] = ranobs;
	obsarray.at(ranobs)++;
      }

      Printf("Generating bootstrap file #%d, first obs=%d\n",isample,obsarray[0]);
      std::stringstream out;
      out << isample;
      filename = TString("newbootstrapsample_").Append(out.str()).Append(".txt");
      filename.Prepend(workingdir);

      if (!fexist((char *)filename.Data())) {
	FILE * pFile;
	pFile = fopen (filename.Data(),"w");
	for (UInt_t iobs=0 ; iobs< nobs ; iobs++) {
	fprintf (pFile, "%9d ",obsarray[iobs]);	  
	}
	fprintf (pFile, "\n");
	fclose (pFile);
      }
      else {
	printf("***Error in TMinLkhd::BootStrap*** Bootstrap sample file %s already exists!! Please delete them before generating new ones.\n\n",filename.Data());
	assert(0);
      }
    }
  }
}


void TMinLkhd::SetData(char * filename, char * vart_file) {

  if (nobs!=0) {
    cout << "Error: Data has already been loaded!!" << endl;
    cout << "Clearing Data and previous models..." << endl;
    ClearDataMembers();
  }

  if (GetMPRank()==0) cout << "Reading in Variable Table: " << vart_file << endl;
  ifstream in1;
  in1.open(vart_file);
  in1 >> nvar;
  if (!in1.good()||(nvar<=0)) {
    cout << "Problem reading nvar in Variable Table" << endl;
    assert(0);
  }

  if (GetMPRank()==0) printf("There are %d variables.\n",nvar);

  var_table.resize(nvar);
  for (UInt_t i = 0 ; i < nvar ; i++) {
    in1 >> var_table.at(i);
    if (!in1.good()) {
      cout << "Problem reading in Variable Table" << endl;
      assert(0);
    }
  }
  in1.close();

//  cout << "Variable Table:" << endl;
//  for (int i = 0 ; i < nvar ; i++) cout << i << ". " << var_table.at(i) << endl;

  ifstream in;
  in.open(filename);

  //First check if data is complete and the number of observations
  Double_t foo;

  in >> foo;
  while (in.good()) {
    for (UInt_t i = 1 ; i < nvar ; i++) {
      in >> foo;
      if (!in.good()) {
	cout << "Problem reading in Data file: data is not a multiple of nvar!!" << endl;
        nobs = 0;
	assert(0);
      }
    }
    nobs++;
    in >> foo;
  }
  in.close();


  if ((nobs>0)&(nvar>0)) {
    in.open(filename);
    if (GetMPRank()==0) cout << "Found " << nobs << " observations in file " << filename << ".\n\n";
    
    data.resize(nobs*nvar);// = new Double_t[nobs*nvar];
    //    data = new Double_t[nobs*nvar];
    // read in data
    for (UInt_t i = 0 ; i < nobs*nvar ; i++)
      in >> data.at(i);
    //	in >> data[i];
    in.close();
    
//     printf("First three are:\nVariable            Obs1     Obs2     Obs3\n");
//     for (UInt_t j = 0 ; j < nvar ; j++) {
//       printf("%4d %20s %8.3f %8.3f %8.3f\n", j+1, var_table.at(j).Data(),data.at(j),data.at(nvar+j),data.at(2*nvar+j));
//       //      printf("%15s %8.3f %8.3f %8.3f\n", var_table.at(j).Data(),data[j],data[nvar+j],data[2*nvar+j]);
//      }

    skipobs.clear();
    skipobs.resize(nobs,0);

    bootstrapobs.clear();
    bootstrapobs.resize(nobs,0);
    //    for (UInt_t iobs = 0 ; iobs < nobs ; iobs++) bootstrapobs[iobs] = iobs;
    for (UInt_t iobs = 0 ; iobs < nobs ; iobs++) bootstrapobs[iobs] = 1;
  }  
}



void TMinLkhd::SetFactorScores() {

  if (nobs==0) {
    cout << "***Error in TMinLkhd::SetFactorScores*** You must add the data file before you can specify the index variable" << endl;
    assert(0);
  }

  // Read in factor-scores if using adaptive integration
  std::vector<Double_t> initv(nfac,0.0);
  fscore.resize(nobs,initv);
  fstderr.resize(nobs,initv);

  TString filename;
    if ((nsubsamples==0)&&(nbootsamples==0)) {
      filename = TString("factor_predictions.txt");
    }
    else {
      std::stringstream out;
      out << GetCurrentSample();
      if (GetCurrentSample()>-1) filename = TString("factor_predictions_").Append(out.str()).Append(".txt");
      else filename = TString("factor_predictions.txt");
    }

    //  TString filename = TString("factor_predictions.txt");
  filename.Prepend(workingdir);
  
  if (fexist((char *)filename.Data())) {
    Int_t obs_id;
    Double_t score, score_err;
    
    ifstream in1;
    in1.open(filename.Data());
    if (GetMPRank()==0) cout << "Reading in factor scores: " << filename.Data() << "\n";
    
    for (UInt_t iobs = 0 ; iobs < nobs; iobs++) {
      in1 >> obs_id;

      if ((indexvar>-1)&&(obs_id!=Int_t(data[iobs*nvar+indexvar]))) {
	  cout << "***ABORTING: Problem reading factor score values\n";
	  printf("IDs do not match! dataid=%d and factorscoreID=%d\n",Int_t(data[iobs*nvar+indexvar]),obs_id);
	  assert(0);
	}
      for (UInt_t ifac = 0; ifac < nfac; ifac++) {
	in1 >> score >> score_err;
	
	if (!in1.good()) {
	  cout << "***ABORTING: Problem reading factor score values\n";
	  assert(0);
	}
	fscore[iobs][ifac] = score;
	fstderr[iobs][ifac] = score_err;
	//      printf("%4d %8.4f %8.4f\n",parnum,lastpar,lastpar_err);
      }
    }
    in1.close();
  }
  else {
    cout << "***ABORTING: " << filename.Data() << " not found, need values of parameters in measurement system!!\n";
    assert(0);
  }
}  



void TMinLkhd::SetObsIndex(const TString index) {

  if (nobs==0) {
    cout << "***Error in TMinLkhd::SetObsIndex*** You must add the data file before you can specify the index variable" << endl;
    assert(0);
  }
  
  indexvar = -1;
  for (UInt_t j = 0 ; j < nvar ; j++) {
    if (!index.CompareTo(var_table.at(j))) {
      indexvar = j;
      break;
    }
  }
  if (indexvar==-1) {
    cout << "***Error in TMinLkhd::SetObsIndex*** Could not find the variable " 
	 << index.Data()  << " in the list of variables. Maybe a typo??"
	 << endl;
    assert(0);
  }
}


void TMinLkhd::UseWeights(const TString weight) {

  if (nobs==0) {
    cout << "***Error in TMinLkhd::UseWeights*** You must add the data file before you can specify the weight variable" << endl;
    assert(0);
  }
  
  weightvar = -1;
  for (UInt_t j = 0 ; j < nvar ; j++) {
    if (!weight.CompareTo(var_table.at(j))) {
      weightvar = j;
      break;
    }
  }
  if (weightvar==-1) {
    cout << "***Error in TMinLkhd::UseWeights*** Could not find the variable " 
	 << weight.Data()  << " in the list of variables. Maybe a typo??"
	 << endl;
    assert(0);
  }
  else printf("**OK, Using weight variable %s for likelihood.\n",weight.Data());
}

void TMinLkhd::SetFactorSpecificQuadPoints(vector<UInt_t> & quadlist)
{
  if (quadlist.size()==nfac) {
    fac_npoints.clear();
    fac_npoints.resize(nfac,nquad_points);
    
    for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {
      if ( (quadlist.at(ifac)>0) & (quadlist.at(ifac)<=nquad_points) & (quadlist.at(ifac)%2==0) ){
	fac_npoints.at(ifac) = quadlist.at(ifac);
      }
      else {
	cout << "***Error in TMinLkhd::AddModelSetFactorSpecificQuadPoints()*** Nquad array contains either non-positive elements, odd elements, or elements that are larger than the maximum specified! Element " << ifac << " = " << quadlist.at(ifac) << endl;
	assert(0);
      }
    }
  }
  else {
    cout << "***Error in TMinLkhd::AddModelSetFactorSpecificQuadPoints()*** Nquad array not the right size!" << endl;
    assert(0);
  }
}


void TMinLkhd::SetAsymmetricQuadPoints(const int asym) {

  
  if ((nquad_points%2==0) && (abs(asym) < (nquad_points/2)) ) {

    //Create arrays for HGQ function
    Double_t * getx = new Double_t[nquad_points + 2*abs(asym)];
    Double_t * getw = new Double_t[nquad_points + 2*abs(asym)];
    
    Double_t sqrt2 = sqrt(2.0);
    printf("*****Getting constants for assymetric %2d GH quadrature points.\n",nquad_points);


    // get positive points 
    calcgausshermitequadrature(nquad_points+2*asym,getx,getw);
    for (UInt_t ipoint = 0 ; ipoint < (nquad_points/2) + asym ; ipoint++) {
      x.at((nquad_points-1)*nquad_points + ipoint) = sqrt2*getx[ipoint];
      w.at((nquad_points-1)*nquad_points + ipoint) = getw[ipoint];
      printf("Assymetric Quadrature points [%2d]= %22.16g %22.16g\n",nquad_points,GetHGQx(nquad_points, ipoint)/sqrt2,GetHGQw(nquad_points, ipoint));
    }

    // get positive points 
    calcgausshermitequadrature(nquad_points-2*asym,getx,getw);
    int getxoffset = nquad_points/2 - asym;
    int xoffset = nquad_points/2 + asym;
    for (UInt_t ipoint = 0 ; ipoint < (nquad_points/2) - asym; ipoint++) {
      x.at((nquad_points-1)*nquad_points + ipoint + xoffset) = sqrt2*getx[ipoint+getxoffset];
      w.at((nquad_points-1)*nquad_points + ipoint + xoffset) = getw[ipoint+getxoffset];
      printf("Assymetric Quadrature points [%2d]= %22.16g %22.16g\n",nquad_points,GetHGQx(nquad_points, ipoint+xoffset)/sqrt2,GetHGQw(nquad_points, ipoint+xoffset));
    }
    
    printf("\n\n");

    delete [] getx;
    delete [] getw;

  }
  else {
      cout << "***Error in TMinLkhd::SetAsymmetricQuadPoints *** Assymetric quadrature point can only be implemented for even number of quadrature points" << endl;
      assert(0);
  }

}

void TMinLkhd::AddModel(const char *name, const char *title, TString modeltype, vector<TString> & moddata,Double_t * normfac, UInt_t nchoice, UInt_t nrank)
{

  if (nprintgroup>9) {
    nprintgroup = 0;
    current_printgroup++;
  }

  if (nobs==0) {
    cout << "***Error in TMinLkhd::AddModel(" << name << ")*** You must add the data file before you can add a model" << endl;
    assert(0);
  }

  // Don't bother adding outcome models if we are just doing the measurement system
  if ((current_estimationgroup>0)&&((stage==1)||(stage==5))) return;

  if (models.size()>0) {
    if ((models.back().GetGroup()>current_estimationgroup) || (models.back().GetGroup()+1<current_estimationgroup)) {
      cout << "***Error in TMinLkhd::AddModel(" << name << ")*** Please add models so groups are in order" << endl;
      cout << name << " Group=" << current_estimationgroup
	   << ", lastgroup=" << models.back().GetGroup() << "\n";
      assert(0);
    }
  }
  else if (current_estimationgroup!=0) {
    cout << "***Error in TMinLkhd::AddModel(" << name << ")*** Please add measurement system first!" << endl;
    assert(0);
  }


  Int_t nparam_thismodel = 0;
  Int_t type = -1;
  if (!modeltype.CompareTo("linear") || !modeltype.CompareTo("Linear")) type = 1;
  if (!modeltype.CompareTo("probit") || !modeltype.CompareTo("Probit")) type = 2;
  if (!modeltype.CompareTo("logit") || !modeltype.CompareTo("Logit")) type = 3;
  if (!modeltype.CompareTo("orderedprobit") || !modeltype.CompareTo("OrderedProbit")) type = 4;
  if (type==-1) {
    cout << "***Error in TMinLkhd::AddModel(" << name << ")*** Could not find model type " 
	 << modeltype.Data() << " in the list of models. Maybe a typo??" 
	 << endl;
    assert(0);
  }

  if (type<3) {
    nchoice = 2;
    nrank = 1;
  }
  if ( ((type==3)||(type==4)) && ((nchoice<2)||(nrank<1)) ) {
    cout << "***Error in TMinLkhd::AddModel(" << name << ")*** Trying to add discrete choice model with less than two categories or less thank two rankings!!\n";
    assert(0);
  }

  if ((type==3)&&(normfac)) {
    cout << "***Error in TMinLkhd::AddModel(" << name << ")*** Normalizing loadings in logit model is not supported.\n";
    assert(0);
  }


  UInt_t arraysize = moddata.size();
  vector<Int_t> intdata(arraysize,-9999);
  if (!moddata[0].CompareTo("none")) intdata[0] = -1;
  if (!moddata[0].CompareTo("type")) intdata[0] = -2;
  if (!moddata[1].CompareTo("none")) intdata[1] = -1;

  for (UInt_t i = 0 ; i < arraysize ; i++) {
    for (UInt_t j = 0 ; j < nvar ; j++) {
      if (!moddata[i].CompareTo(var_table.at(j))) {
	intdata[i] = j;
	break;
      }
    }
    if (intdata[i]==-9999) {
      cout << "***Error in TMinLkhd::AddModel(" << name << ")*** Could not find the variable " 
	   << moddata[i].Data()  << " in the list of variables. Maybe a typo??"
	   << endl;
      assert(0);
    }
  }

  //Count up number of parameters for this model:
  // number of regressors:
  nparam_thismodel += arraysize - 2;

  // number of unnormalized unobservables:
  UInt_t totunobs = nfac;

  // Add type "loadings" to params if this is not the type probability model
  if (intdata[0]!= -2) totunobs += ntyp;
  
  if ((normfac) && ((nfac>0) || (ntyp>0)) ) {
    for (UInt_t i = 0; i < totunobs ; i++) {
      if (normfac[i]<-9998) nparam_thismodel++;
      //      else if ((norm_models[i]<0)&&(fabs(normfac[i])>0.001)) norm_models[i] = models.size();
    }
  }
  else nparam_thismodel += totunobs;

  //Look for first normalization model in measurement system
  //Also set norm_models = -1 if at least one model loads on factor
  if ((current_estimationgroup==0) && (nfac>0)) {
    if (normfac) {
      for (UInt_t i = 0; i < nfac ; i++) {
	//Record first non-zero loading as normalization
	if ( (norm_models[i]<0) && (fabs(normfac[i])>0.001) && (normfac[i]>-9998)) norm_models[i] = models.size();
      }
    }
    else {
      //If a model has no normalization, then set norm_models = -1 to indicate there is at least one model that loads on factor
      for (UInt_t i = 0; i < nfac ; i++) {
	if (norm_models[i]==-2) norm_models[i] = -1;
      }
    }
  }

  
  // variance of error term
  if (type==1) nparam_thismodel++;

  // need to multiply parameters by the number of choices-1 for logit model
  if (type==3) nparam_thismodel *= (nchoice - 1);

  // Choice thresholds for ordered probit model
  if (type==4) nparam_thismodel += (nchoice - 1);

  TModel newmodel(name, title, type,current_estimationgroup,
		  current_printgroup,intdata,nfac,ntyp,normfac,nchoice,nrank);
  models.push_back(newmodel);
  nparam_models.push_back(nparam_thismodel);
  fparam_models.push_back(nparam); //first parameter of this model
  
  nprintgroup++;

  nparam += nparam_thismodel;
  parconstrained.resize(nparam,-1);
  param_fixval.resize(nparam,-9999.0);  
}


void TMinLkhd::AddTypesModel(vector<TString> & typedata, Double_t * typenorm)
{
  // Keep track of which model defines the type probablities
  type_model = models.size();

  //Add type probability model to list of models
  typedata.insert(typedata.begin(), TString("type"));
  AddModel("TypeProb", "Type Probability", "logit", typedata, typenorm, ntyp+1);

}



void TMinLkhd::ConstrainFactorCorrelations() { 
  for (UInt_t imix = fac_nmix-1 ; imix > 0 ; imix--) {
    parconstrained.at(Getifvar(imix,2)) = Getifvar(0,2);
  }
}

// void TMinLkhd::ConstrainLastBetaToModel(const Int_t imod, const TString cov) {
//   // Check that target is not already constrained
//   printf("
// }


// This function equalizes factor loadings between two models
void TMinLkhd::ConstrainLastFactorLoadingToModel(const UInt_t targetmod, UInt_t ifac) {
  // Check that loading is not fixed
  UInt_t imodels[2], ipar[2];

  imodels[0] = targetmod;
  imodels[1] = models.size()-1;
  
  for (int imod = 0 ; imod < 2 ; imod++) {
    if (models.at(imodels[imod]).GetNorm().size()!=0) {
      if (models.at(imodels[imod]).GetNorm().at(ifac-1)>-9998) {
	//ABORT loading is normalized!
	cout << "***Error in TMinLkhd::ConstrainLastFactorLoadingToModel*** Trying to constrain loading to a fixed parameter!" << endl;
	assert(0);
      }
      else if ((ifac>1)&&(models.at(imodels[imod]).GetNorm().at(ifac-2)>-9998)) {
	ipar[imod] = fparam_models.at(imodels[imod]) + models.at(imodels[imod]).GetNreg() + ifac-2;	
      }
      else ipar[imod] = fparam_models.at(imodels[imod]) + models.at(imodels[imod]).GetNreg() + ifac-1;
    }
    else ipar[imod] = fparam_models.at(imodels[imod]) + models.at(imodels[imod]).GetNreg() + ifac-1;
  }
  parconstrained.at(ipar[1]) = ipar[0];
}


void TMinLkhd::LastModel_SetEndogenousReg(TString endogvar) {

  // Get var number:
  int endogvarnum = -9999;
  for (UInt_t ivar = 0 ; ivar < nvar ; ivar++) {
    if (!endogvar.CompareTo(var_table.at(ivar))) {
      endogvarnum = ivar;
      break;
    }
  }
  if (endogvarnum==-9999) {
    cout << "***Error in TMinLkhd::LastModel_SetEndogenousReg(" << models.back().GetName() << ")*** Could not find the variable " 
	 << endogvar.Data()  << " in the list of variables. Maybe a typo??"
	 << endl;
    assert(0);
  }

  models.back().AddEndogenousReg(endogvarnum);
  
}


void TMinLkhd::LastModel_SetRankShareVar(TString sharevar) {

  // Get var number:
  int sharevarnum = -9999;
  for (UInt_t ivar = 0 ; ivar < nvar ; ivar++) {
    if (!sharevar.CompareTo(var_table.at(ivar))) {
      sharevarnum = ivar;
      break;
    }
  }
  if (sharevarnum==-9999) {
    cout << "***Error in TMinLkhd::LastModel_SetRankShareVar(" << models.back().GetName() << ")*** Could not find the variable " 
	 << sharevar.Data()  << " in the list of variables. Maybe a typo??"
	 << endl;
    assert(0);
  }

  models.back().SetRankShareVar(sharevarnum);
}


void TMinLkhd::SetSimModelOrder(std::vector<TString> & modlist) {
  
  sim_modelorder.clear();
  if (type_model<9999) sim_modelorder.push_back(type_model);

  for (UInt_t ilist = 0 ; ilist < modlist.size(); ilist++) {
  
    // Get var number:
    int modvarnum = -9999;
    for (UInt_t ivar = 0 ; ivar < nvar ; ivar++) {
      if (!modlist.at(ilist).CompareTo(var_table.at(ivar))) {
	modvarnum = ivar;
	break;
      }
    }
    if (modvarnum==-9999) {
      cout << "***Error in TMinLkhd::SetSimModelOrder*** Could not find the variable " 
	   << modlist.at(ilist).Data()  << " in the list of variables. Maybe a typo??"
	   << endl;
      assert(0);
    }
    
    int modnum = -9999.0;
    for (UInt_t imod = 0 ; imod < models.size() ; imod++) {
      if (models.at(imod).GetOutcome()==modvarnum) {
	modnum = imod;
	break;
      }
    }
    if (modnum==-9999) {
      cout << "***Error in TMinLkhd::SetSimModelOrder*** Could not find model with outcome " 
	   << modlist.at(ilist).Data()  << " in the list of models. Maybe a typo??"
	   << endl;
      assert(0);
    }
    sim_modelorder.push_back(modnum);
  }

  //Check that all models have been included:
  for (UInt_t imod = 0 ; imod < models.size() ; imod++) {
    int foundmodel = 0;
    for (UInt_t jmod = 0 ; jmod < sim_modelorder.size() ; jmod++) {
      if (sim_modelorder.at(jmod)==imod) {
	foundmodel = 1;
	break;
      }
    }
    if (foundmodel==0) sim_modelorder.push_back(imod);
  }

  if (models.size()!=sim_modelorder.size()) {
      cout << "***Error in TMinLkhd::SetSimModelOrder*** Model lists don't match! WTF?" 
	   << endl;
      assert(0);
  }

  //Print order of simulating models
  if (GetMPRank()==0) {
    printf("Simulation order of models set to:\n");

    for (UInt_t jmod = 0 ; jmod < sim_modelorder.size() ; jmod++) {
      printf("%d: %s\n",jmod,models.at(sim_modelorder.at(jmod)).GetName());
    }
  }
}



void TMinLkhd::LastModel_FixParamValue(TString fixvar, double value, int choice) {

  if (models.size()==0) {
    cout << "***Error in TMinLkhd::LastModel_FixParamValue*** Add Models first!" << endl;
    assert(0);
  }
  int thismodel = models.size() - 1;


  // Check that choice is not greater than nchoice
  int choice_offset = 0; // Only need to change if multinomial logit

  if ( (models[thismodel].GetType()==3) && (models.back().GetNchoice()>2) ) {
    if ( (choice>1) && (choice <= models[thismodel].GetNchoice()) ) choice_offset = choice-2;
    else {
      cout << "***Error in TMinLkhd::LastModel_FixParamValue(" << models[thismodel].GetName() << ")*** choice "
	   << choice  << " is larger than nchoice=" << models[thismodel].GetNchoice()
	   << ". Maybe a typo??"
	   << endl;
      assert(0);
    }
  }


  // Get var number:
  int fixvarnum = -9999;
  for (UInt_t ivar = 0 ; ivar < nvar ; ivar++) {
    if (!fixvar.CompareTo(var_table.at(ivar))) {
      fixvarnum = ivar;
      break;
    }
  }
  if (fixvarnum==-9999) {
    cout << "***Error in TMinLkhd::LastModel_FixParamValue(" << models[thismodel].GetName() << ")*** Could not find the variable " 
	 << fixvar.Data()  << " in the list of variables. Maybe a typo??"
	 << endl;
    assert(0);
  }

  //Find modparoffset
  int modparoffset = -9999;
  vector <Int_t> regs = models[thismodel].GetReg();
  for (Int_t ireg =  0; ireg <  models[thismodel].GetNreg() ; ireg++) {
    if (regs.at(ireg) == fixvarnum) {
      modparoffset = ireg;
      break;
    }
  }
  if (modparoffset==-9999) {
    cout << "***Error in TMinLkhd::LastModel_FixParamValue(" << models[thismodel].GetName() << ")*** Could not find the variable " 
	 << fixvar.Data()  << " in the list of regressors. Maybe a typo??"
	 << endl;
    assert(0);
  }

  int fixparam = fparam_models[thismodel] + choice_offset*(models[thismodel].GetNreg()+nfac+ntyp) + modparoffset;
  FixParamValue(fixparam,value);

}

//void TMinLkhd::LastModel_SetEndogenousMajor(UInt_t modelN) {
//
//  if (modelN<models.size()) {
//    models.back().SetEndogenousMajor(modelN,models,var_table);  
//  }
//  else {
//    cout << "***Error in TMinLkhd::LastModel_SetEndogenousReg*** Model number is bigger than the size of the model array!" << endl;
//    assert(0);
//  }
//}


void TMinLkhd::LastModel_Splitsim(TString splitvar) {

  Int_t intdata = -9999;
  for (UInt_t j = 0 ; j < nvar ; j++) {
    if (!splitvar.CompareTo(var_table.at(j))) {
      intdata = j;
      break;
    }
  }
  if (intdata==-9999) {
    cout << "***Error in TMinLkhd::Splitsim*** Could not find the variable " 
	 << splitvar.Data()  << " in the list of variables. Maybe a typo??"
	 << endl;
    assert(0);
  }
 

  models.back().SplitSim(intdata);
}

void TMinLkhd::ResetFitInfo() {

  //  printf("Starting TMinLkhd::ResetFitInfo()\n");
  //Check that type model exists if the number of types is greater than 0
  if (ntyp > 0) {
    if (type_model==9999) {
      cout << "ERROR (TMinLkhd::ResetFitInfo): ntyp>0, but the model for types was not specified!" << endl;
      assert(0);
    }    
  }

  
  counter =0;
  param.clear();
  param_err.clear();
  parfixed.clear();
  nfreeparam = nparam;

  EstGroup_lkhd.clear();

  if (param.size()<nparam) {
    param.resize(nparam,-9999);
    param_err.resize(nparam,-9999);
    parfixed.resize(nparam,0);
  }

  // Set parameters that are permanently fixed
  for (UInt_t ipar = 0 ; ipar < nparam ; ipar++) {
    if (param_fixval.at(ipar)>-9998.0) {
      param.at(ipar) = param_fixval.at(ipar);
      param_err.at(ipar) = -9999.0;
      FixParPerm(ipar);
    }
  }
  if (GetMPRank()==0) printf("***TMinLkhd::ResetFitInfo() with %d parameters and %d free.\n",int(param.size()),nfreeparam);


  //Remove chains of constraints from parconstrained
  for (UInt_t i = 0 ; i < nparam ; i++) {
    if (parconstrained.at(i)>-1) {
      Int_t target = parconstrained.at(i);
      Int_t counter = 0;
      while ( (parconstrained.at(target)>=0)&&(parconstrained.at(target)<Int_t(nparam))
	      &&(counter<10) ) {
	target = parconstrained.at(target);
	counter++;
      }
      if (parconstrained.at(target)==-1)
	parconstrained.at(i) = target;
      else if (counter>=10) {
	printf("***ERROR:ResetFitInfo() Found chain of constraints longer than 9 links: possible loop?\n");
	assert(0);
      }
      else {
	printf("***ERROR:ResetFitInfo() par #%d constrained to target par #%d, which is not valid!\n",i,target);
	assert(0);
      }
    }
  }


  for (UInt_t imod = 0 ; imod < models.size() ; imod++) {
    UInt_t obscounter = 0;

    vector <Int_t> regs = models[imod].GetReg();
    for (UInt_t iobs = 0 ; iobs < nobs ; iobs++) {

      if (models[imod].GetMissing()==-1) {
	obscounter++;

	if (imod!=type_model) {
	  if ((fabs(data.at(iobs*nvar + models[imod].GetOutcome())+9999)<0.1)&&(models[imod].GetNrank()==1)) {
	    cout << "ERROR (TMinLkhd::ResetFitInfo): Found bad outcome in data!"
		 << " In model " <<  models[imod].GetTitle() << "\n"
		 << "Looking at outcome #" << models[imod].GetOutcome() 
		 << ", obs#:" << iobs
		 << ", Value=" << data.at(iobs*nvar + models[imod].GetOutcome()) << endl;
	    assert(0);
	  }
	  
	  if (models[imod].GetType()==2) {
	    if (( int(data.at(iobs*nvar+models[imod].GetOutcome()))!=1)
		&&( int(data.at(iobs*nvar+models[imod].GetOutcome()))!=0)){
	      cout << "ERROR (TMinLkhd::ResetFitInfo): Found non-binary number for probit outcome!"
		   << " In model " <<  models[imod].GetTitle() << "\n"
		   << "Looking at outcome #" << models[imod].GetOutcome() 
		   << ", obs#:" << iobs
		   << ", Value=" << data.at(iobs*nvar + models[imod].GetOutcome()) << endl;
	      assert(0);
	    }
	  }
	}

	//Now check regressors:
	for (UInt_t ireg = 0 ; ireg < regs.size(); ireg++) {
	  if (fabs(data.at(iobs*nvar + regs.at(ireg))+9999)<0.1) {
	    cout << "ERROR (TMinLkhd::ResetFitInfo): Found bad regressor in data!"
		 << " In model " <<  models[imod].GetTitle() << "\n"
		 << "Looking at variable #" <<  var_table.at(regs.at(ireg)) 
		 << ", obs#:" << iobs
		 << ", Value=" << data.at(iobs*nvar + regs.at(ireg)) << endl;
	    assert(0);
	  }
	}
      }

      else if (data.at(iobs*nvar + models[imod].GetMissing())==1) {
	obscounter++;
	if ((fabs(data.at(iobs*nvar + models[imod].GetOutcome())+9999)<0.1)&&(models[imod].GetNrank()==1)) {
	  cout << "ERROR (TMinLkhd::ResetFitInfo): Found bad outcome in data!"
	       << " In model " <<  models[imod].GetTitle() << "\n"
	       << "Looking at outcome: " <<  var_table.at(models[imod].GetOutcome())
	       << ", Missing (" << models[imod].GetMissing() << ")=" << data.at(iobs*nvar+models[imod].GetMissing())
	       << ", obs#:" << iobs
	       << ", Value=" << data.at(iobs*nvar + models[imod].GetOutcome()) << endl;
	  assert(0);
	}
	
	if (models[imod].GetType()==2) {
	  if (( int(data.at(iobs*nvar + models[imod].GetOutcome()))!=1)
	      &&( int(data.at(iobs*nvar + models[imod].GetOutcome()))!=0)) {
	    cout << "ERROR (TMinLkhd::ResetFitInfo): Found non-binary number for probit outcome!"
		 << " In model " <<  models[imod].GetTitle() << "\n"
		 << "Looking at outcome: " <<  var_table.at(models[imod].GetOutcome())
		 << ", Missing (" << models[imod].GetMissing() << ")=" << data.at(iobs*nvar+models[imod].GetMissing())
		 << ", obs#:" << iobs
		 << ", Value=" << data.at(iobs*nvar + models[imod].GetOutcome()) << endl;
	    assert(0);
	  }
	}

	if (((models[imod].GetType()==3)||(models[imod].GetType()==4))&&(models[imod].GetNrank()==1)) {
	  int validvalue = 0;
	  for (int ichoice = 1 ; ichoice <= models[imod].GetNchoice() ; ichoice++) {
	    if (int(data.at(iobs*nvar + models[imod].GetOutcome())) == ichoice) validvalue = 1;
	  }
	  if (validvalue==0) {
	    cout << "ERROR (TMinLkhd::ResetFitInfo): Found value that is not allowed for multinomial outcome!"
		 << " In model " <<  models[imod].GetTitle() << "\n"
		 << "Looking at outcome: " <<  var_table.at(models[imod].GetOutcome())
		 << ", Missing (" << models[imod].GetMissing() << ")=" << data.at(iobs*nvar+models[imod].GetMissing())
		 << ", obs#:" << iobs
		 << ", Value=" << data.at(iobs*nvar + models[imod].GetOutcome()) 
		 << ", but nchoices=" << models[imod].GetNchoice() << endl;
	    assert(0);
	  }
	}

	//Now check regressors:
	for (UInt_t ireg = 0 ; ireg < regs.size(); ireg++) {
	  if (fabs(data.at(iobs*nvar + regs.at(ireg))+9999)<0.1) {
	    cout << "ERROR (TMinLkhd::ResetFitInfo): Found bad regressor in data!"
		 << " In model " <<  models[imod].GetTitle() << "\n"
		 << "Looking at regressor: " << var_table.at(regs.at(ireg))
		 << ", Missing (" << models[imod].GetMissing() << ")=" << data.at(iobs*nvar+models[imod].GetMissing())
		 << ", obs#:" << iobs
		 << ", Value=" << data.at(iobs*nvar + regs.at(ireg)) << endl;
	    assert(0);
	  }
	}
      }

    }
    nobs_models.push_back(obscounter);
  }

  // Initialize factor distribution
  // nvariance * nmix
  //(nmix - 1)*nfac means
  //nmix - 1 weights
  if (nfac>0) {
    for (UInt_t imix = 0 ; imix < fac_nmix ; imix++) {
      
      Int_t offset = imix*(f_nvariance + nfac + 1);
      //variances
      for (UInt_t ivar = 0; ivar < f_nvariance ; ivar++) {

	// Set param only if it isn't fixed
	  if (ivar < nfac) {
	    Setparam(offset + ivar,1.0);
	    Setparam_err(offset+ivar, 0.5);
	  }
	  else {
	    Setparam(offset + ivar, 0.2);
	    Setparam_err(offset + ivar, 0.1);
	  }
	}

      if (imix < fac_nmix-1) {
	// means
	for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {
	  Setparam(offset + f_nvariance + ifac,  0.5);
	  Setparam_err(offset + f_nvariance + ifac,  0.1);
	}
	
	// weights: =0.0 means all equally weighted
	  Setparam(offset + f_nvariance + nfac, 0.0);
	  Setparam_err(offset + f_nvariance + nfac, 0.1);
      }
    }
  }


  UInt_t ipar = nfac_param;

  vector <Double_t> fnorm(nfac,-9999);
  vector<Double_t> normout_mn(nfac,-9999);

  //  cout << "Calculating mean of normalized models:";
  
  // Check that there is a normalization for each factor
  for (UInt_t ifac =0; ifac < nfac; ifac++) {
    if (norm_models.at(ifac) < 0 ) {
      if ((GetMPRank()==0) && (norm_models.at(ifac)==-1) ) cout << "Factor #" << ifac << " is not normalized: ";
      if ((GetMPRank()==0) && (norm_models.at(ifac)==-2) ) cout << "Factor #" << ifac << " is not in measurement system: ";

      if (fac_nmix>1) {
	printf("Not able to use mixtures without normalizing a loading\n");
	assert(0);
      }
      else {
	// If using normal distribution fix variance to 1.0:
	if (GetMPRank()==0) printf("Fixing variance to 1.0 for factor #%d\n",ifac);
	param_fixval.at(ifac) = 1.0;
	param.at(ifac) = 1.0;
	param_err.at(ifac) = 0.0;

	FixParPerm(ifac);
	//	parfixed.at(ifac) = 1;
	//	nfreeparam--;

	fnorm.at(ifac) = 1.0;
	normout_mn.at(ifac) = 0.0;
      }
      //      assert(0);
    }
    else {
      fnorm.at(ifac) = (models.at(norm_models.at(ifac)).GetNorm()).at(ifac);
      
      Double_t sumdt = 0.0;
      UInt_t ncount = 0;
      
      Int_t imod = norm_models.at(ifac);
      
      for (UInt_t iobs = 0 ; iobs < nobs ; iobs++) {
	if (models.at(imod).GetMissing()==-1) {
	  sumdt += data.at(iobs*nvar + models.at(imod).GetOutcome());
	  ncount++;
	}
	else {
	  if (data.at(iobs*nvar + models.at(imod).GetMissing())==1) {
	    sumdt += data.at(iobs*nvar + models.at(imod).GetOutcome());
	    ncount++;
	  }
	}
      }
      normout_mn.at(ifac) = sumdt/ncount;
      //      printf("ifac#%d(%d)=%f, ",ifac,ncount,normout_mn.at(ifac));
    }
  }
  //  printf("\n\n");


  for (UInt_t imod = 0 ; imod < models.size() ; imod++) {

    //    printf("Calculating stuff for model %d outcome=%d, missing=%d, type=%d, nchoice=%d\n",imod,models[imod].GetOutcome(), models[imod].GetMissing(), models[imod].GetType(), models[imod].GetNchoice());
    //Scale initial values differently depending on model
    Double_t thismult = 1.0;
    if (models[imod].GetType()==3) thismult = 0.1;
    
    // If this is a multinomial logit, we want to do this for each choice if nchoice>2:
    int nchoice = 2; // In this case we don't need to do anything special even if it is logit
    if (models[imod].GetType()==3) nchoice = models[imod].GetNchoice();
    for (int ichoice = 2 ; ichoice <= nchoice ; ichoice++) {
      vector <Int_t> regs = models[imod].GetReg();
      
      Double_t sumdt = 0.0;
      Double_t sumdtsq = 0.0;
      Double_t sumoutreg = 0.0;
      
      vector<Double_t> sumnormout(nfac,0.0);
      vector<Int_t> normoutcount(nfac,0);
      
      // First get mean and sd of outcome
      UInt_t ncount = 0;
      Double_t outcome_sd = 1.0;
      Double_t outcome_mn = 0.0;
      
      if (imod!=type_model) {
	for (UInt_t iobs = 0 ; iobs < nobs ; iobs++) {
	  double outcome = -9999.0;
	  UInt_t lastcount = ncount;
	  if (models[imod].GetMissing()==-1) {
	    outcome = data.at(iobs*nvar + models[imod].GetOutcome());
	    if ((nchoice>2)&&(models[imod].GetType()==3)) outcome = (outcome==ichoice);
	    sumdt += outcome;
	    sumdtsq += outcome*outcome;
	    ncount++;
	  }
	  else {
	    if (data.at(iobs*nvar + models[imod].GetMissing())==1) {
	      outcome = data.at(iobs*nvar + models[imod].GetOutcome());
	      if ((nchoice>2)&&(models[imod].GetType()==3)) outcome = (outcome==ichoice);
	      sumdt += outcome;
	      sumdtsq += outcome*outcome;
	      ncount++;
	    }
	  }
	
	  if (ncount>lastcount) {
	    for (UInt_t ifac = 0 ; ifac < nfac; ifac++) {
	      if (norm_models.at(ifac) > -1) {
		if (models.at(norm_models.at(ifac)).GetMissing()==-1) {
		  sumnormout.at(ifac) += outcome*data.at(iobs*nvar + models[norm_models.at(ifac)].GetOutcome());
		  normoutcount.at(ifac)++;
		}
		else {
		  if (data.at(iobs*nvar + models.at(norm_models.at(ifac)).GetMissing())==1) {
		    sumnormout.at(ifac) += outcome*data.at(iobs*nvar + models[norm_models.at(ifac)].GetOutcome());
		    normoutcount.at(ifac)++;
		  }
		}
	      }
	    }
	  }
	}
    
      
	outcome_sd = sqrt((sumdtsq - sumdt*sumdt/ncount)/ncount);
	outcome_mn = sumdt/ncount;
      }
      if (models[imod].GetType()!=1) outcome_sd = 1.0;
      //      printf("Outcome SD for model %d = %f, %f, %f, %d\n",imod,outcome_sd, sumdtsq, sumdt,ncount);
      
      //Get sign for loading in this model
      vector<Double_t> facloadsign(nfac,0.0);
      for (UInt_t ifac = 0 ; ifac < nfac; ifac++) {
	// If there is no covariance because of conditioning, then just use the same sign normalization
	if (normoutcount.at(ifac)>0) facloadsign.at(ifac) = (fnorm.at(ifac)*(sumnormout.at(ifac)/normoutcount.at(ifac) - normout_mn.at(ifac)*outcome_mn) < 0) ? -1.0 : 1.0;
	else facloadsign.at(ifac) = (fnorm.at(ifac) < 0) ? -1.0 : 1.0;
	//      if ((imod<4)&&(ifac==1)) printf("
      }
      
      // Initialize regressors Get mean, sd and covariance with outcome
      int countconstcovariates = 0;
      for (UInt_t ireg = 0 ; ireg < regs.size(); ireg++) {
	sumdt = 0.0;
	sumdtsq = 0.0;
	sumoutreg = 0.0;
	ncount = 0;
	if (imod!=type_model) {
	  for (UInt_t iobs = 0 ; iobs < nobs ; iobs++) {

	    if (models[imod].GetMissing()==-1) {
	      double outcome = data.at(iobs*nvar + models[imod].GetOutcome());
	      if ((nchoice>2)&&(models[imod].GetType()==3)) outcome = (outcome==ichoice);
	      sumdt += data.at(iobs*nvar + regs[ireg]);
	      sumdtsq += data.at(iobs*nvar + regs[ireg])*data.at(iobs*nvar + regs[ireg]);
	      sumoutreg += data.at(iobs*nvar + regs[ireg])*outcome;
	      ncount++;
	    }
	    else {
	      if (data.at(iobs*nvar + models[imod].GetMissing())==1) {
		double outcome = data.at(iobs*nvar + models[imod].GetOutcome());
		if ((nchoice>2)&&(models[imod].GetType()==3)) outcome = (outcome==ichoice);
		sumdt += data.at(iobs*nvar + regs[ireg]);
		sumdtsq += data.at(iobs*nvar + regs[ireg])*data.at(iobs*nvar + regs[ireg]);
		sumoutreg += data.at(iobs*nvar + regs[ireg])*outcome;
		ncount++;
	      }
	    }
	  }
	}
	
	Double_t reg_sd = 0.0;
	Double_t reg_mn = 0.0;
	Double_t covsign = 0.0;
	
	if (ncount>0) {
	  reg_sd = sqrt((sumdtsq - sumdt*sumdt/ncount)/ncount);
	  reg_mn = sumdt/ncount;
	  covsign = (sumoutreg/ncount - reg_mn*outcome_mn < 0) ? -1.0 : 1.0;
	}
	// Set intercept or beta0
	if ((reg_sd<0.001)&&(imod!=type_model)) {
	  //	printf("Setting par%d=%f\n",ipar,outcome_mn/reg_mn);
	  if (fabs(reg_mn-1.0) < 0.001) {
	    if (countconstcovariates==0) {
	      Setparam(ipar, thismult*loadingMultiplier*outcome_mn/reg_mn);
	      Setparam_err(ipar, 0.1*thismult*loadingMultiplier*outcome_mn/reg_mn);
	      countconstcovariates++;
	    }
	    else {
	      countconstcovariates++;
	      Setparam(ipar, 0.0);
	      Setparam_err(ipar, -9999.0);
	      FixParPerm(ipar);
	    }
	  }
	  else {
	    Setparam(ipar, 0.0);
	    Setparam_err(ipar, -9999.0);
	    FixParPerm(ipar);
	  }
	}
	else if (imod==type_model) {
	  Setparam(ipar, 0.1*loadingMultiplier);
	  Setparam_err(ipar, 0.1*loadingMultiplier); 
	}
	else Setparam(ipar, thismult*loadingMultiplier*covsign*outcome_sd/reg_sd/regs.size());
	Setparam_err(ipar, 0.5*fabs(param[ipar]));
	if (param_err[ipar]<0.1*outcome_sd/regs.size()) Setparam_err(ipar, 0.1*outcome_sd/regs.size());

	ipar++;
      }
      
      // Initialize factor loadings and precision and thresholds
      vector <Double_t> thisnorm = models[imod].GetNorm();
      for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {
	if (thisnorm.size()==0) {
	  Setparam(ipar, thismult*loadingMultiplier*facloadsign.at(ifac)*0.1);
	  Setparam_err(ipar, 0.01);
	  ipar++;
	}
	else {
	  if (thisnorm.at(ifac)<-9998) {
	    Setparam(ipar, thismult*loadingMultiplier*facloadsign.at(ifac)*0.1);
	    Setparam_err(ipar, 0.01);
	    ipar++;
	  }
	}
      }
      // Now set loadings for types
      if (imod!=type_model) {
	for (UInt_t itype = 0 ; itype < ntyp ; itype++) {
	  if (thisnorm.size()==0) {
	    Setparam(ipar, thismult*loadingMultiplier*0.1);
	    Setparam_err(ipar, 0.01);
	    ipar++;
	  }
	  else {
	    if (thisnorm.at(nfac + itype)<-9998) {
	      Setparam(ipar, thismult*loadingMultiplier*0.1);
	      Setparam_err(ipar, 0.01);
	      ipar++;
	    }
	  }
	}
      }

      if (models[imod].GetType()==1) {
	if (initVariance>0) {
	  Setparam(ipar, initVariance);
	  Setparam_err(ipar, 0.5*initVariance);
	}
	else {
	 Setparam(ipar, outcome_sd);
	 Setparam_err(ipar, 0.5*outcome_sd);
	}
	ipar++;
      }
      if (models[imod].GetType()==4) {
	Setparam(ipar, -2.0 + 4.0/double(models[imod].GetNchoice()));
	Setparam_err(ipar, 1.0);
	ipar++;
	
	for (int ithresh = 2 ; ithresh < models[imod].GetNchoice() ; ithresh++) {
	  Setparam(ipar, 4.0/double(models[imod].GetNchoice()));
	  Setparam_err(ipar, 1.0);
	  ipar++;
	}
      }
    } // loop over choices of logit
  } // loop over models
  //   for (int i = 0 ; i < param.size() ; i++) printf("par%d = %f\n",i,param[i]);

  if (fscore.size()==nobs) {
    SetFactorScores();
  }

  
  if (GetMPRank()==0) printf("***TMinLkhd::ResetFitInfo() with %d parameters and %d free.\n",int(param.size()),nfreeparam);
}

Int_t TMinLkhd::Minimize(Int_t printlevel) {
  // Check that the prerequisites are in place
  
  // Reset the fit parameters and results
  ResetFitInfo();
  if (stage==0) return TestCode(printlvl);
  else if (stage==1) return Est_measurementsys(printlvl);
  else if (stage==2) return Est_outcomes(printlvl);
  else if (stage==3) {
    Int_t err = Est_measurementsys(printlvl);
    if (err==0) err = Est_outcomes(printlvl);
    return err;
  }
  else if (stage==4) return Simulate(printlvl);
  else if (stage==5) return PredictFactors(printlvl);
  else if (stage==6) return GenerateCHsamples(printlvl); // generate "bootstrap" samples for Cameron-Heckman standard errors
  else {
    cout << "******ERROR: non-existent stage=" << stage << "\n";
    return -1;
  }
}

Int_t TMinLkhd::TestCode(Int_t printlevel) {

  TString filename = "init_par.txt";
  filename.Prepend(workingdir);
  if (fexist(filename.Data())) {
    Int_t parnum;
    Double_t lastpar;
    Double_t lastpar_err;
    
    ifstream in1;
    in1.open(filename.Data());
    cout << "Reading in initial values for " << nparam << " parameters.\n";
    
    for (UInt_t i = 0 ; i < nparam; i++) {
      in1 >> parnum >> lastpar >> lastpar_err;
      if (!in1.good()) {
	cout << "Problem reading Initial parameter values\n";
	assert(0);
      }
      printf("%4d %8.4f %8.4f\n",parnum,lastpar,lastpar_err);
 
      Setparam(i, lastpar);
      Setparam_err(i, lastpar_err);
    }
    in1.close();
  }
  cout << endl << "********Initial parameter values:" << endl;
  PrintParam();
  //  cout << "point 1" << endl;

//   TStopwatch timer;
//   //Check likelihood at initial parameters
//   Double_t fvalue = 0;
//   Double_t * thisparam = new Double_t[nparam];
//   Double_t * grad = new Double_t[nparam];
//   Int_t num_param = nparam;
  
//   for (UInt_t ipar = 0 ; ipar < nparam ;ipar++) thisparam[ipar] = param[ipar];
//   //  cout << "point 2" << endl;
//   timer.Reset();
//   timer.Start();
//   LkhdFcn(num_param,grad,fvalue,thisparam,2);
//   //  EvalLkhd(fvalue, grad, parinit, 2,0,1);
//   timer.Stop();
//   cout << "Likelihood=" << fvalue << endl;
//   // cout << "Gradiant:" << endl;
//   //  for (UInt_t i = 0 ; i < nparam; i++) cout << i << " "<< grad[i] << endl;
//   Double_t cputime = timer.CpuTime();
//   printf("Calculating the Likelihood took %8.4f seconds.\n",cputime);
  
  Int_t nparam_min = nparam;
  //Doing timing tests:
  Double_t fvalue = 0;
  Double_t * thisparam = new Double_t[nparam_min];
  Double_t * grad = new Double_t[nparam_min];
  Double_t * hess = new Double_t[nparam_min*(nparam_min+1)/2];
  Double_t cputime;
  TStopwatch timer;

  //  cout << "ok evaluating likelihood now!\n";
  Int_t ifreepar = 0;
  for (UInt_t ipar = 0 ; ipar < nparam ;ipar++) {
    if (parfixed[ipar]==0) {
      thisparam[ifreepar] = param[ipar];
      ifreepar++;
    }
  }
  Int_t repeat = 1;

  timer.Reset();
  timer.Start();
  for (int i = 0 ; i < repeat ; i++) LkhdFcn(nparam_min,grad,fvalue,thisparam,1,hess );
  timer.Stop();
   cputime = timer.CpuTime();
  printf("Likelihood at initial par=%8.4f.\n",fvalue);
  printf("Calculating the Likelihood took %8.4f seconds.\n",cputime/Double_t(repeat));

  timer.Reset();
  timer.Start();
  for (int i = 0 ; i < repeat ; i++) LkhdFcn(nparam_min,grad,fvalue,thisparam,2,hess);
  timer.Stop();
   cputime = timer.CpuTime();
  printf("Calculating the Gradient   took %8.4f seconds.\n",cputime/Double_t(repeat));

  if ((CalcHess ==1) || (initializing==1)) {
    timer.Reset();
    timer.Start();
    for (int i = 0 ; i < repeat ; i++) LkhdFcn(nparam_min,grad,fvalue,thisparam,3,hess);
    timer.Stop();
    cputime = timer.CpuTime();
    printf("Calculating the Hessian    took %8.4f seconds.\n",cputime/Double_t(repeat));
  }
  
  delete [] thisparam;
  delete [] grad;
  delete [] hess;


  // ****************************************************************
  // ok try KNITRO:
  // ************************************************************

    /*---- CREATE A NEW KNITRO SOLVER INSTANCE. */

//   int status;
  
//   KTR_context_ptr  kc;
//   kc = KTR_new();
//   if (kc == NULL)
//     {
// 	cout << "*** KTR_new failed, maybe a license issue?\n";
// 	exit( EXIT_FAILURE );
//     }
  
//   status = KTR_load_param_file (kc, "knitro.opt");


//   if (status != 0)
//     {
//       cout << "*** KTR_load_param_file() returned " << status << "\n";
//       return( false );
//     }
//   else cout << "**knitro options set\n";

//     int       _nN;
//     double *  _daXInit = NULL;
//     double *  _daXLo = NULL;
//     double *  _daXUp = NULL;
//     int       _nM;
//     int    *  _naCType = NULL;
//     double *  _daCLo = NULL;
//     double *  _daCUp = NULL;
//     int       _nNnzJ;
//     int    *  _naJacIndexVars = NULL;
//     int    *  _naJacIndexCons = NULL;
//     int       _nNnzH;
//     int    *  _naHessCols = NULL;
//     int    *  _naHessRows = NULL;


//     _nN = nparam; // number of parameters
//     _nM = 0; // number of constraints

//     //---- VARIABLES ARE BOUNDED FROM BELOW.
//     _daXLo  = new double[_nN];
//     _daXUp  = new double[_nN];

//     // Tell MINUIT about the parameter values:
//     for (UInt_t ipar = 0 ; ipar < nparam ; ipar++) {
//       // With limits:
//       // set limits for covariances for first mixture
//       if ((fac_corr!=0) && (ipar>=nfac) && (ipar<f_nvariance)) {
// 	_daXLo[ipar] = -1.0;
// 	_daXUp[ipar] = 1.0;

// // 	_daXLo[ipar] = -0.8; // What we used in Minuit
// // 	_daXUp[ipar] = 0.8;
//       }
//       // set limits for covariances for second mixture
//       else if ((fac_nmix>1) && (fac_corr!=0) && (ipar>=f_nvariance+nfac+1+nfac) && (ipar<f_nvariance+nfac+1+f_nvariance)) {
// 	_daXLo[ipar] = -1.0;
// 	_daXUp[ipar] = 1.0;
// // 	_daXLo[ipar] = -0.8;
// // 	_daXUp[ipar] = 0.8;
//       }
//       // set limits for covariances for third mixture
//       else if ((fac_nmix>2) && (fac_corr!=0) && (ipar>=2*(f_nvariance+nfac+1)+nfac) && (ipar<2*(f_nvariance+nfac+1)+f_nvariance)) {
// 	_daXLo[ipar] = -1.0;
// 	_daXUp[ipar] = 1.0;
// // 	_daXLo[ipar] = -0.8;
// // 	_daXUp[ipar] = 0.8;

//       }
//       else {
//         _daXLo[ipar] = -KTR_INFBOUND;
//         _daXUp[ipar] = KTR_INFBOUND;
//       }
//     }

//     //---- THE CONSTRAINTS IS A LINEAR INEQUALITY.
//     //---- PUT THE CONSTANT TERM IN THE RIGHT-HAND SIDE.
//     // No constraints here
// //     _naCType  = new int[_nM];
// //     _daCLo    = new double[_nM];
// //     _daCUp    = new double[_nM];
// //     _daCLo[0] = -KTR_INFBOUND;
// //     _daCUp[0] = 3.0;
// //     _naCType[0] = KTR_CONTYPE_LINEAR;

//     //---- SPECIFY THE CONSTRAINT JACOBIAN SPARSITY STRUCTURE.
//     _nNnzJ = 0;
// //     _nNnzJ = 3;
// //     _naJacIndexVars = new int[_nNnzJ];
// //     _naJacIndexCons = new int[_nNnzJ];
// //     _naJacIndexCons[ 0] = 0;  _naJacIndexVars[ 0] = 0;
// //     _naJacIndexCons[ 1] = 0;  _naJacIndexVars[ 1] = 1;
// //     _naJacIndexCons[ 2] = 0;  _naJacIndexVars[ 2] = 2;

//     //---- SPECIFY THE HESSIAN OF THE LAGRANGIAN SPARSITY STRUCTURE.
//     //    _nNnzH = 0;
//     _nNnzH = nparam*(nparam+1)/2;
//     _naHessRows = new int[_nNnzH];
//     _naHessCols = new int[_nNnzH];
//     int counter = 0;
//     for (int i = 0 ; i < nparam ; i++) {
//       for (int j = i ; j < nparam ; j++) { 
// 	_naHessRows[counter] = i;   
// 	_naHessCols[counter] = j;
// 	counter++;
//       }
//     }

// //     _naHessRows[0] = 0;  _naHessCols[0] = 0;
// //     _naHessRows[1] = 0;  _naHessCols[1] = 1;
// //     _naHessRows[2] = 0;  _naHessCols[2] = 2;
// //     _naHessRows[3] = 1;  _naHessCols[3] = 1;
// //     _naHessRows[4] = 2;  _naHessCols[4] = 2;

//     //---- INITIAL GUESS FOR x AND lambda.
//     _daXInit = new double[_nN];
//     double *  daLambdaInit = new double[_nM + _nN];
//     for (int i = 0; i < _nN; i++)
//         _daXInit[i] = param[i];
//     for (int i = 0; i < _nM + _nN; i++)
//         daLambdaInit[i] = 0.0;
//     status = KTR_init_problem (kc, _nN,
//                           KTR_OBJGOAL_MINIMIZE, KTR_OBJTYPE_QUADRATIC,
//                           _daXLo, _daXUp,
//                           _nM, _naCType, _daCLo, _daCUp,
//                           _nNnzJ, _naJacIndexVars, _naJacIndexCons,
//                           _nNnzH, _naHessRows, _naHessCols,
//                           _daXInit, daLambdaInit);

//     delete [] _daXLo;
//     delete [] _daXUp;
//     delete [] _naCType;
//     delete [] _daCLo;
//     delete [] _daCUp;
//     delete [] _naJacIndexVars;
//     delete [] _naJacIndexCons;
//     delete [] _naHessRows;
//     delete [] _naHessCols;
//     delete [] daLambdaInit;

//     if (status != 0)
//         {
//         cout << "*** KTR_init_problem() returned " << status << "\n";
//         return( false );
//         }
//     else cout << "**knitro initialized\n";

//     double * jac = new double[_nN];
//     double * daLambda = new double[_nN];
//     double * c = new double[_nN];
//     status =1 ;
//     int counter = 0;
//     while ((status==1)||(status==2)) {
//       LkhdFcn(num_param,grad,fvalue,_daXInit,status);
//       //  EvalLkhd(fvalue, grad, parinit, 2,0,1);
      
//       status = KTR_check_first_ders (kc, _daXInit, 2, 2.0e-5, 2.0e-5,
//   				     0, fvalue, c, grad, jac, NULL);
//       cout << "Pass " << counter << ". status=" << status << "\n";
//       counter++;
//     }


//      cout << "ok now minimizing\n";
//      status =1 ;
//      int counter = 0;
//      while ((status==1)||(status==2)) {
//        cout << "Eval Lkhd. counter= " << counter << "\n";
//        if (counter>0) LkhdFcn(num_param,grad,fvalue,_daXInit,status);
//        //  EvalLkhd(fvalue, grad, parinit, 2,0,1);
    
//        cout << "KTR_solve. counter= " << counter << "\n";
//        status = KTR_solve (kc, _daXInit, daLambda, 0, &fvalue,
//                                    c, grad, jac, NULL, NULL, NULL);

//  //       status = KTR_check_first_ders (kc, thisparam, 2, 2.0e-6, 2.0e-6,
//  // 				     0, fvalue, c, grad, jac, NULL);
//        cout << "Pass " << counter << ". status=" << status << "\n";
//        counter++;
//      }

//      cout << "Minimization done. Status=" << status << "\n";

//    KTR_free (&kc);

//     /*---- SOLVE THE PROBLEM.
//      *----
//      *---- RETURN STATUS CODES ARE DEFINED IN "knitro.h" AND DESCRIBED
//      *---- IN THE KNITRO MANUAL.
//      */
//     nStatus = KTR_solve (kc, x, lambda, 0, &obj,
//                          NULL, NULL, NULL, NULL, NULL, NULL);

//     printf ("\n\n");
//     if (nStatus != 0)
//         printf ("KNITRO failed to solve the problem, final status = %d\n",
//                 nStatus);
//     else
//     {
//         /*---- AN EXAMPLE OF OBTAINING SOLUTION INFORMATION. */
//         printf ("KNITRO successful, feasibility violation    = %e\n",
//                 KTR_get_abs_feas_error (kc));
//         printf ("                   KKT optimality violation = %e\n",
//                 KTR_get_abs_opt_error (kc));
//     }


//     /*---- DELETE THE KNITRO SOLVER INSTANCE. */

//     free (x);
//     free (lambda);



//      //Setup minimizer
//      if (gMinuit) delete gMinuit;
//      gMinuit = new TMinuit(nparam);
//      gMinuit->SetObjectFit(this);
//      gMinuit->SetPrintLevel(3);
//      //    gMinuit->SetPrintLevel(printlevel);
//      gMinuit->SetFCN(MinuitLkhdFcn);
//      //  gMinuit->SetErrorDef(1.0); // ChiSq
//      gMinuit->SetErrorDef(.5);  // ln(likelihood)


//       // Tell MINUIT about the parameter values:
//       for (UInt_t ipar = 0 ; ipar < nparam ; ipar++) {
//         // With limits:
//         // set limits for covariances for first mixture
//         if ((fac_corr!=0) && (ipar>=nfac) && (ipar<f_nvariance)) {
//           gMinuit->DefineParameter(ipar, "", param[ipar], 0.01*param_err[ipar], -0.8 , 0.8);
//         }
//         // set limits for covariances for second mixture
//         else if ((fac_nmix>1) && (fac_corr!=0) && (ipar>=f_nvariance+nfac+1+nfac) && (ipar<f_nvariance+nfac+1+f_nvariance)) {
//           gMinuit->DefineParameter(ipar, "", param[ipar], 0.01*param_err[ipar], -0.8 , 0.8);
//         }
//         // set limits for covariances for third mixture
//         else if ((fac_nmix>2) && (fac_corr!=0) && (ipar>=2*(f_nvariance+nfac+1)+nfac) && (ipar<2*(f_nvariance+nfac+1)+f_nvariance)) {
//           gMinuit->DefineParameter(ipar, "", param[ipar], 0.01*param_err[ipar], -0.8 , 0.8);
//         }
//         //       else if (ipar<nfac) {
//         // 	gMinuit->DefineParameter(ipar, "", param[ipar], 0.01*param_err[ipar], 0.01 , 10.0);
//         //       }
//         else gMinuit->DefineParameter(ipar, "", param[ipar], 0.01*param_err[ipar], 0,0);
//         //    without limits:
//         //    gMinuit->DefineParameter(ipar, "", param[ipar], param_err[ipar], 0,0);
//       }

//       // Set up other stuff
//       Double_t arglist[16];
//       Int_t ierflg = 0;
//       //arglist[0] = 1;
//       //  gMinuit->mnexcm("CALL FCN", arglist, 1, ierflg);
//       //  arglist[0] = 2;
//       arglist[0] = 2;
//       gMinuit->mnexcm("SET STR", arglist, 1, ierflg);
//       //arglist[0] = 2.5e-8;
//       //gMinuit->mnexcm("SET EPS", arglist, 1, ierflg);
//       //  if (printlevel > 0) 
//       gMinuit->mnexcm("SHO EPS", arglist, 0, ierflg);

//       arglist[0] = 1;
//       //gMinuit->mnexcm("SET GRA", arglist, 1, ierflg);
//       //For checking the gradient:
//       gMinuit->mnexcm("SET GRA", arglist, 0, ierflg);
//       gMinuit->mnexcm("SHO GRA",arglist,0,ierflg);

//      //     gMinuit->mnhess();

//     //Calculate the Hessian
//      //     gMinuit->mnexcm("HESSE", arglist, 0, ierflg);
// //     if (printlevel > -1)
// //       std::cout << "Status from HESSE: " << ierflg << std::endl;


       TestGradient();

       TestHessian();

 //   cout << "ok not calculating gradients" << endl;
 //   arglist[0] = 9000000;  //max number of iterations
 //   arglist[1] = .00001;   // tolerance *10^-5
 //     gMinuit->mnexcm("SET NOGRA", arglist, 0, ierflg);
 //     gMinuit->mnexcm("SHO GRA",arglist,0,ierflg);

 //     gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  return 0;
}

Int_t TMinLkhd::TestGradient(){ 
  //  Double_t machprec = 1e-16;
  //  Double_t delta = pow(machprec,(1.0/3.0));
  //  Double_t delta = sqrt(machprec);
  Double_t delta = 1e-7;

  Double_t fvalue = 0;

  UInt_t numberparam = nparam;
  if (predicting==1) numberparam = nfac;
  Double_t * thisparam = new Double_t[numberparam];
  Double_t * grad = new Double_t[numberparam];
  Int_t num_param = numberparam;

  if (predicting==0) {
    Int_t ifreepar = 0;
    for (UInt_t ipar = 0 ; ipar < numberparam ;ipar++) {
      if (parfixed[ipar]==0) {
	thisparam[ifreepar] = param[ifreepar];
	ifreepar++;
      }
    }
  }
  else {
    for (UInt_t ifac = 0 ; ifac < numberparam ;ifac++) thisparam[ifac] = 0.5;
  }
  LkhdFcn(num_param,grad,fvalue,thisparam,2);

  cout << "***** Testing gradient *******\n";

  Int_t ifreepar = 0;
  for (UInt_t ipar = 0; ipar < numberparam ; ipar++) {
    if (parfixed[ipar]==0) {
      Double_t h = delta*(fabs(thisparam[ifreepar])+1.0);
    
      Double_t oldparam = thisparam[ifreepar];
      
      thisparam[ifreepar] = oldparam - h;
      LkhdFcn(num_param,NULL,fvalue,thisparam,1);
      Double_t f1 = fvalue;
      
      thisparam[ifreepar] = oldparam + h;
      LkhdFcn(num_param,NULL,fvalue,thisparam,1);
      thisparam[ifreepar] = oldparam;
    
      Double_t thisgrad = (fvalue - f1)/(2.0*h);
      Double_t diff = (thisgrad - grad[ifreepar])/thisgrad;

      printf("par %5i: calc grad = %11.4e | finite grad = %11.4e | percent diff = %11.4e\n",ipar, grad[ifreepar], thisgrad, diff);
      ifreepar++;
    }
  }

  cout << "***** Finished Testing gradient *******\n";
  return 0;
}

Int_t TMinLkhd::TestHessian(){ 
  //  Double_t machprec = 1e-16;
  //  Double_t delta = pow(machprec,(1.0/3.0));
  //  Double_t delta = sqrt(machprec);
  Double_t deltag = 1e-7;
  //  Double_t deltaf = 1e-3;

  UInt_t numberparam = nparam;
  if (predicting==1) numberparam = nfac;

  Double_t fvalue = 0;
  Double_t * thisparam = new Double_t[numberparam];
  Double_t * testparam = new Double_t[numberparam];
  Double_t * grad = new Double_t[numberparam];
  Double_t * defgrad = new Double_t[numberparam];
  Double_t * hess = new Double_t[numberparam*(numberparam+1)/2];
  Int_t num_param = numberparam;

  Int_t ifreepar = 0;
  for (UInt_t ipar = 0 ; ipar < numberparam ;ipar++) {
    if (predicting==0) {
      if (parfixed[ipar]==0) {
	thisparam[ifreepar] = param[ipar];
	testparam[ifreepar] = param[ipar];
	ifreepar++;
      }
    }
    else {
      thisparam[ipar] = 0.5;
      testparam[ipar] = 0.5;
    }
  }

  LkhdFcn(num_param,defgrad,fvalue,thisparam,3,hess);

  LkhdFcn(num_param,defgrad,fvalue,thisparam,2);
  //  Double_t defF = fvalue;

  cout << "***** Testing Hessian *******\n";
  Int_t hessterm = 0;
  ifreepar = 0;
  Int_t jfreepar = 0;
  for (UInt_t ipar = 0; ipar < numberparam ; ipar++) {
    if (parfixed[ipar]==0) {
      cout << "Checking row " << ipar << "\n";
      jfreepar = 0;
      for (UInt_t jpar = ipar; jpar < numberparam ; jpar++) {
	if (parfixed[jpar]==0) {
	  Double_t origparami = thisparam[ifreepar];
	  Double_t origparamj = thisparam[jfreepar];

	  // First modify ipar
	  Double_t hig = deltag*(fabs(origparami)+1.0);
	  Double_t hjg = deltag*(fabs(origparamj)+1.0);

      
	  //       thisparam[ipar] = origparami - hig;
	  //       LkhdFcn(num_param,grad,fvalue,thisparam,2);
	  //       Double_t g1 = grad[jpar];
      
	  thisparam[ifreepar] = origparami + hig;
	  LkhdFcn(num_param,grad,fvalue,thisparam,2);

	  thisparam[ifreepar] = origparami;
      
	  //      Double_t Hessijgc = (grad[jpar] - g1)/(4.0*hig);

	  Double_t Hessijgf = (grad[jfreepar]-defgrad[jfreepar])/(2.0*hig);

	  // Now modify jpar
	  
	  //       thisparam[jpar] = origparamj - hjg;
	  //       LkhdFcn(num_param,grad,fvalue,thisparam,2);
	  //       g1 = grad[ipar];
	  
	  thisparam[jfreepar] = origparamj + hjg;
	  LkhdFcn(num_param,grad,fvalue,thisparam,2);

	  thisparam[jfreepar] = origparamj;
      
//      Hessijgc += (grad[ipar] - g1)/(4.0*hjg);

	  Hessijgf += (grad[ifreepar]-defgrad[ifreepar])/(2.0*hjg);

//       // now calculate with function calls
//       Double_t hif = deltaf*(fabs(origparami)+1.0);
//       Double_t hjf = deltaf*(fabs(origparamj)+1.0);
//       Double_t Hessijfc = 0.0;
//       if (ipar==jpar) {
// 	thisparam[ipar] = origparami + 2.0*hif;
// 	LkhdFcn(num_param,grad,fvalue,thisparam,1);
// 	Hessijfc += -1.0*fvalue;

// 	thisparam[ipar] = origparami + hif;
// 	LkhdFcn(num_param,grad,fvalue,thisparam,1);
// 	Hessijfc += 16.0*fvalue;
	
// 	thisparam[ipar] = origparami;
// 	LkhdFcn(num_param,grad,fvalue,thisparam,1);
// 	Hessijfc += -30.0*fvalue;

// 	thisparam[ipar] = origparami - hif;
// 	LkhdFcn(num_param,grad,fvalue,thisparam,1);
// 	Hessijfc += 16.0*fvalue;

// 	thisparam[ipar] = origparami - 2.0*hif;
// 	LkhdFcn(num_param,grad,fvalue,thisparam,1);
// 	Hessijfc += -1.0*fvalue;

// 	thisparam[ipar] = origparami;

// 	Hessijfc *= 1.0/(12.0*hif*hif);
//       }
//       else {

// 	thisparam[ipar] = origparami + hif;
// 	thisparam[jpar] = origparamj + hjf;
// 	LkhdFcn(num_param,grad,fvalue,thisparam,1);
// 	Hessijfc += fvalue;

// 	thisparam[ipar] = origparami + hif;
// 	thisparam[jpar] = origparamj - hjf;
// 	LkhdFcn(num_param,grad,fvalue,thisparam,1);
// 	Hessijfc += -1.0*fvalue;

// 	thisparam[ipar] = origparami - hif;
// 	thisparam[jpar] = origparamj + hjf;
// 	LkhdFcn(num_param,grad,fvalue,thisparam,1);
// 	Hessijfc += -1.0*fvalue;

// 	thisparam[ipar] = origparami - hif;
// 	thisparam[jpar] = origparamj - hjf;
// 	LkhdFcn(num_param,grad,fvalue,thisparam,1);
// 	Hessijfc += fvalue;

// 	thisparam[ipar] = origparami;
// 	thisparam[jpar] = origparamj;

// 	Hessijfc *= 1.0/(4.0*hif*hjf);
//       }

//       Double_t Hessijff = 0.0;
//       if (ipar!=jpar) {
// 	thisparam[ipar] = origparami + hif;
// 	thisparam[jpar] = origparamj + hjf;
//       }
//       else thisparam[ipar] = origparami + 2.0*hif;

//       LkhdFcn(num_param,grad,fvalue,thisparam,1);
//       Hessijff += fvalue;
//       thisparam[ipar] = origparami;
//       thisparam[jpar] = origparamj;
      
//       thisparam[ipar] = origparami + hif;
//       LkhdFcn(num_param,grad,fvalue,thisparam,1);
//       Hessijff += -1.0*fvalue;
//       thisparam[ipar] = origparami;      

//       thisparam[jpar] = origparamj + hjf;
//       LkhdFcn(num_param,grad,fvalue,thisparam,1);
//       Hessijff += -1.0*fvalue;
//       thisparam[jpar] = origparamj;      

//       Hessijff += defF;

//       Hessijff *= 1.0/(hif*hjf);

	  Double_t diff;
	  if (fabs(hess[hessterm]) > 1e-7) diff = fabs(Hessijgf - hess[hessterm])/fabs(hess[hessterm]);
	  else if (fabs(Hessijgf)<1e-7)      diff = fabs(Hessijgf - hess[hessterm]);
	  else                             diff = fabs(Hessijgf - hess[hessterm])/fabs(Hessijgf);
	  if (diff < 1e-10) diff = 1e-10;
	  
	  //	  if (diff > 1e-4)
	    printf("Hessian[%3i][%3i]: calc = %11.4e | finite = %11.4e | percent diff = %11.4e\n",ipar,jpar, hess[hessterm], Hessijgf, diff);
	  //      printf("Hessian[%3i][%3i]: calc = %11.4e | finite = %11.4e | percent diff = %11.4e\n",ipar,jpar, hess[hessterm], Hessijgf, diff);
      
//       Double_t diff;
//       if (fabs(Hessijfc) > 1e-7) diff = (Hessijgc - Hessijfc)/fabs(Hessijfc);
//       else if (fabs(Hessijgc)<1e-7)      diff = Hessijgc - Hessijfc;
//       else                             diff = (Hessijgc - Hessijfc)/fabs(Hessijgc);
//       printf("Hessian[%3i][%3i]: central fun = %11.4e | central grad = %11.4e | percent diff = %11.4e\n",ipar,jpar, Hessijfc, Hessijgc, diff);
      
//       Double_t diff;
//       if (fabs(hess[hessterm]) > 1e-7) diff = (Hessijfc - hess[hessterm])/fabs(hess[hessterm]);
//       else if (fabs(Hessijfc)<1e-7)      diff = Hessijfc - hess[hessterm];
//       else                             diff = (Hessijfc - hess[hessterm])/fabs(Hessijfc);
//       printf("Hessian[%3i][%3i]: calc = %11.4e | finite = %11.4e | percent diff = %11.4e\n",ipar,jpar, hess[hessterm], Hessijfc, diff);
      
	  hessterm++;
	  jfreepar++;
	}
      }
      ifreepar++;
    }
  }

  for (UInt_t ipar = 0 ; ipar < numberparam ;ipar++) {
    if (thisparam[ipar] != testparam[ipar]) 
      cout << "Shit param #" << ipar << " changed!\n";
  }


  cout << "***** Finished Testing Hessian *******\n";
  return 0;
}


Int_t TMinLkhd::Est_measurementsys(Int_t printlevel) {
  cout << "***** Estimating Measurement System *******\n";

  // outcome models are not added if stage==0, so we don't have to worry about them here...
  // might change this in the future
  TStopwatch timer;
  //Check likelihood at initial parameters
  Int_t ierflg = 0;
  Double_t fvalue = 0;
  Double_t * thisparam = new Double_t[nparam];
  Double_t * grad = new Double_t[nparam];
  Double_t * finallkhds = new Double_t[fac_nmix];

  Int_t num_param = nparam;

  // Find first outcome model
  Int_t first_outmodel=-1;
  for (UInt_t imod = 0 ; imod < models.size() ; imod++) {
    if ((models[imod].GetGroup()==1)&&(first_outmodel==-1)) first_outmodel = imod;
  }
  
  //Ignore all outcome models in first stage and fix their parameters
  if (first_outmodel!=-1) {
    for (UInt_t imod = first_outmodel ; imod < models.size() ; imod++) SetIgnoreMod(imod);
    for (UInt_t i = fparam_models[first_outmodel] ; i < nparam ; i++) FixPar(i);
    num_param = fparam_models[first_outmodel];
  }
  
  if (first_outmodel==-1) first_outmodel = models.size();
  
  //**** NEED TO FIX THIS ALSO  -- Fix what? 27.03.2023
  
  //Make factor loading list (for all models?)
  // so we can easily scale initial values of loadings. Not used for anything else
  vector<Bool_t> factorloadinglist(num_param,0);
  for (UInt_t imod = 0 ; int(imod) < first_outmodel ; imod++) {
    int nchoice = 1; // In this case we don't need to do anything special even if it is logit
    if (models[imod].GetType()==3) nchoice = models[imod].GetNchoice() - 1;
    for (int ichoice = 0 ; ichoice < nchoice ; ichoice++) {
      Int_t imodFirstLoading = fparam_models[imod] + ichoice*(models[imod].GetNreg()+nfac+ntyp*(imod!=type_model))+ models[imod].GetNreg();
      UInt_t nfreefac = nparam_models[imod] - models[imod].GetNreg() - ntyp*(imod!=type_model);
      if (models[imod].GetType()==1) nfreefac--;
      if (models[imod].GetType()==4) nfreefac-= (models[imod].GetNchoice()-1);
      if (models[imod].GetType()==3) nfreefac = nfac;
      
      for (UInt_t ifac = 0; ifac < nfreefac ; ifac++) {
	factorloadinglist.at(imodFirstLoading+ifac) = 1;
      }
    }
  }
  

  // Read in initial values. If init_par.txt exists, use that.
  // If bootstrapping or subsampling get "fullmodel_par.txt"
  TString filename;
  if ((nsubsamples!=0)||(nbootsamples!=0)) filename = "fullmodel_par.txt";
  else filename = "init_par.txt";
  filename.Prepend(workingdir);
  if (fexist(filename.Data())) {
    Int_t parnum;
    Double_t lastpar;
    Double_t lastpar_err;
    
    ifstream in1;
    in1.open(filename.Data());
    cout << "Reading in " << num_param << " initial values\n";
    
    for (Int_t i = 0 ; i < num_param; i++) {
      in1 >> parnum >> lastpar >> lastpar_err;
      if (!in1.good()) {
	cout << "Problem reading Initial parameter values\n";
	assert(0);
      }
      //      printf("%4d %8.4f %8.4f\n",parnum,lastpar,lastpar_err);

      
      if ((factorloadinglist.at(i)==1)&&(initFixedLoadings==0)) {
	Setparam(i, lastpar*loadingMultiplier);
	Setparam_err(i, lastpar_err);
      }
      else if (factorloadinglist.at(i)==0) {
	Setparam(i, lastpar);
	Setparam_err(i, lastpar_err);
      }
    }
    in1.close();
  }
  else if ((nsubsamples!=0)||(nbootsamples!=0)) {
    cout << "Problem reading fullmodel parameter values\n";
    assert(0);
  }
  
  //  if ( (!fexist(filename.Data())) || (initBetaOLS==1) ) {
  if (initBetaOLS==1) {
    cout << "running individual regressions for initial parameters...\n";
    
    // Now loop through each model and get the initial values for beta
    initializing=1;
    for (Int_t imod = 0 ; imod < first_outmodel ; imod++) SetIgnoreMod(imod);
    
    Int_t firstpar = nfac_param;
    for (UInt_t imod = 0 ; int(imod) < first_outmodel ; imod++) {
      if (imod!=type_model) {
	cout << "**************Minimizing Model #" << imod << "\n";
	RemoveIgnoreMod(imod);
	
	for (Int_t i = 0 ; i < num_param ; i++) {
	  if ((i<firstpar) || (i>=firstpar+models[imod].GetNreg())) FixPar(i);// parfixed[i]=1;
	  else ReleasePar(i); //parfixed[i]=0;
	}
	if (models[imod].GetType()==1) ReleasePar(firstpar+nparam_models[imod]-1); 
	
	if ( (models[imod].GetType()==3) && (models[imod].GetNchoice()>2) ) {
	  for (int ichoice = 1 ; ichoice < models[imod].GetNchoice()-1 ; ichoice++) {
	    for (Int_t i =  0; i <  models[imod].GetNreg() ; i++) {
	      ReleasePar((firstpar + ichoice*(models[imod].GetNreg()+nfac+ntyp*(imod!=type_model))) + i);
	    }
	  }
	}
	if (models[imod].GetType()==4) {
	  for (Int_t i = 0 ; i < models[imod].GetNchoice()-1 ; i++) {
	    ReleasePar(firstpar+nparam_models[imod]-1-i); 
	  }
	}
	//     ierflg = Min_Minuit(-1);
	if (printlvl>0) {
	  PrintParam(imod+1);
	  ierflg = Min_Ipopt(1);
	  PrintParam(imod+1);
	}
	else ierflg = Min_Ipopt(0);
	if ((ierflg!=0)&&(nsubsamples==0)&&(nbootsamples==0)) assert(0);
	
	SetIgnoreMod(imod);
	firstpar += nparam_models[imod];
      }
    }
    initializing=0;
    for (Int_t imod = 0 ; imod < first_outmodel ; imod++) RemoveIgnoreMod(imod);
    for (Int_t i = 0 ; i < num_param ; i++) ReleasePar(i); //parfixed[i] = 0;
  }
  
  
  if (printlvl>0) cout << endl << "********Initial parameter values:" << endl;
  if (printlvl>0) PrintParam(0);

  Int_t ifreepar = 0;
  for (Int_t ipar = 0 ; ipar < num_param ;ipar++) {
    if (parfixed[ipar]==0) {
      thisparam[ifreepar] = param[ipar];
      ifreepar++;
    }
  }

  cout << "Ok, now minimizing all parameters!\n";
  LkhdFcn(ifreepar,grad,fvalue,thisparam,1);
  cout << "Initial Likelihood=" << fvalue << endl;
  
  timer.Reset();
  timer.Start();
  
  if ((EstSequentMix == 1) && (fac_nmix!=1)) {
    cout << "\n**********Estimating mixtures sequentially*********\n";
    cout << "\n ****Estimating with 1 mixture:\n";
    SetEstNMix(1);
    //    ierflg = Min_Minuit(0);
    ierflg += Min_Ipopt(1);
    
    if (ierflg == 0 ) { 
      
      Int_t ifreepar = 0;
      for (UInt_t ipar = 0 ; ipar < nparam ;ipar++) {
	if (parfixed[ipar]==0) {
	  thisparam[ifreepar] = param[ipar];
	  ifreepar++;
	}
      }
      LkhdFcn(ifreepar,grad,fvalue,thisparam,1);
      cout << "Final Likelihood=" << fvalue << endl;
      finallkhds[0] = fvalue;
      
      SetEstNMix(2);
      for (UInt_t ivar = 0 ; ivar < f_nvariance ; ivar++) {
	Setparam(Getifvar(1,ivar), Getfvar(0,ivar));
      }
      for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {
	Setparam(Getifmean(0,ifac), 0.5*Getfvar(0,ifac));
      }
      
      cout << "\n ****Estimating with 2 mixtures:\n";
      if (printlvl>0) cout << endl << "********Initial parameter values:" << endl;
      if (printlvl>0) PrintParam(0);

       if (fac_corr!=0) {
 	cout << "\n ****First we will fix the correlations\n";
	//First we need to send parameters to all the nodes. Yes this is kludgy:
	for (Int_t ipar = 0 ; ipar < num_param ;ipar++) thisparam[ipar] =  param[ipar];
	LkhdFcn(num_param,grad,fvalue,thisparam,1);

 	for (UInt_t imix = 0 ; imix < 2; imix++) FixPar(Getifvar(imix,2));
       	ierflg = Min_Ipopt(1);
 	for (UInt_t imix = 0 ; imix < 2; imix++) ReleasePar(Getifvar(imix,2));
       }
      
      if (fac_nmix>2) {
	//      ierflg = Min_Minuit(0);
      	ierflg = Min_Ipopt(1);

	if (ierflg == 0 ) { 
	  ifreepar = 0;
	  for (UInt_t ipar = 0 ; ipar < nparam ;ipar++) {
	    if (parfixed[ipar]==0) {
	      thisparam[ifreepar] = param[ipar];
	      ifreepar++;
	    }
	  }
	  LkhdFcn(ifreepar,grad,fvalue,thisparam,1);
	  cout << "Final Likelihood=" << fvalue << endl;
	  finallkhds[1] = fvalue;
	  
	  SetEstNMix(3);
	  for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {
	    Setparam(Getifvar(2,ifac),  (Getfvar(0,ifac) > Getfvar(1,ifac)) ? Getfvar(0,ifac) : Getfvar(1,ifac));
	    Setparam(Getifmean(1,ifac),  -1.0*Getfmean(0,ifac));
	  }
	  
	  cout << "\n ****Estimating with 3 mixtures:\n";
	  cout << endl << "********Initial parameter values:" << endl;
	  if ((nsubsamples==0)&&(nbootsamples==0)) PrintParam(0);
	}
      }
    }
  }

  SetEstNMix(fac_nmix);

  //ierflg = Min_Minuit(0);
  //std::cout << "Status from Minuit: " << ierflg << std::endl;
  if (ierflg==0) {
    // ierflg = Min_Minuit(0);
    ierflg = Min_Ipopt(1);

    ifreepar = 0;
    for (Int_t ipar = 0 ; ipar < num_param ;ipar++) {
      if (parfixed[ipar]==0) {
	thisparam[ifreepar] = param[ipar];
	ifreepar++;
      }
    }
    //    for (Int_t ipar = 0 ; ipar < num_param ;ipar++) thisparam[ipar] = param[ipar];
    LkhdFcn(ifreepar,grad,fvalue,thisparam,1);
    cout << "Final Likelihood=" << fvalue << endl;
    finallkhds[fac_nmix-1] = fvalue;

    EstGroup_lkhd.push_back(fvalue);
  }


  std::cout << "Status from Ipopt: " << ierflg << std::endl;
  timer.Stop();
  
  if (ierflg!=0) {
    cout << "*******MINIMIZATION FAILED!!********" << endl;
    if (printlvl>0) PrintParam(0);
    
    //    for (Int_t ipar = 0 ; ipar < num_param ;ipar++) thisparam[ipar] = param[ipar];

    ifreepar = 0;
    for (Int_t ipar = 0 ; ipar < num_param ;ipar++) {
      if (parfixed[ipar]==0) {
	thisparam[ifreepar] = param[ipar];
	ifreepar++;
      }
    }
    LkhdFcn(ifreepar,grad,fvalue,thisparam,2);
    if (printlvl>0) cout << "Last Likelihood=" << fvalue << endl;
    if (printlvl>0) cout << "Gradient= \n";
    for (Int_t ipar = 0 ; ipar < num_param ;ipar++) printf("%4d: %8.4f\n",ipar,grad[ipar]);
    
    FILE * pFile;
    TString filename;
    if ((nsubsamples==0)&&(nbootsamples==0)) {
      filename = TString("last_par.txt");
    }
    else {
      std::stringstream out;
      out << GetCurrentSample();
      filename = TString("last_par_").Append(out.str()).Append(".txt");
    }

    filename.Prepend(workingdir);
    pFile = fopen (filename.Data(),"w");
    for (Int_t ipar=0 ; ipar< num_param ; ipar++)
      {
	fprintf (pFile, "%5d %12.8f %12.8f\n",ipar,Getparam(ipar),Getparam_err(ipar));
      }
    fprintf (pFile, "\n");
    fclose (pFile);
    
    //	assert(0);
  }
  else {
    
    FILE * pFile;
    TString filename;
    if ((nsubsamples==0)&&(nbootsamples==0)) {
      filename = TString("meas_par.txt");
    }
    else {
      std::stringstream out;
      out << GetCurrentSample();
      filename = TString("meas_par_").Append(out.str()).Append(".txt");
    }
    filename.Prepend(workingdir);
    pFile = fopen (filename.Data(),"w");
    for (Int_t ipar=0 ; ipar< num_param ; ipar++)
      {
	fprintf (pFile, "%5d %12.8f %12.8f\n",ipar,Getparam(ipar),Getparam_err(ipar));
      }
    fprintf (pFile, "\n");
    fclose (pFile);
    
    if (printlvl>0) {
      cout << endl << "********Parameter values:" << endl;
      PrintParam(0);
      cout << "Minimization took " << timer.CpuTime() << " s." << endl;
      if ((EstSequentMix == 1) && (fac_nmix!=1)) {
	cout << "\n Estimated nmix sequentially. Likelihoods are:\n";
	for (UInt_t i = 0 ; i < fac_nmix ; i++) printf("%4d: %12.8f \n",i+1,finallkhds[i]);
      }
    }
  }
  

  for (UInt_t imod = 0 ; imod < models.size() ; imod++) RemoveIgnoreMod(imod);
  for (UInt_t i = 0 ; i < nparam ; i++) ReleasePar(i);
  
  delete [] thisparam;
  delete [] grad;
  delete [] finallkhds;
  
  return ierflg;
}


Int_t TMinLkhd::Est_outcomes(Int_t printlevel) {

  TString filename;

  cout << "***** Estimating Outcomes *******\n";

  //Some checks that factor scores are loaded when they are needed:
  if  ((initEstOutcomeLoadings)||(adapt_int)) {
    if (fscore.size()!=nobs){
      cout << "***ABORTING: Factor scores are not loaded when they are needed for estimation. "
	   << "initEstOutcomeLoadings=" << initEstOutcomeLoadings
	   << " and adapt_int=" << adapt_int << "\n";
      assert(0);
    }
  }
  
  // Get first outcome model and parameter
  Int_t first_outmodel=-1;
  Int_t last_group = 0;
  for (UInt_t imod = 0 ; imod < models.size() ; imod++) {
    if ((models[imod].GetGroup()==1)&&(first_outmodel==-1)) {
      first_outmodel = imod;
    }
    if (models[imod].GetGroup()>last_group) last_group=models[imod].GetGroup();
  }

  if (first_outmodel==-1) {
    cout << "ERROR: NO OUTCOME GROUPS!\n";
    return 1;
  }
  // Read in the values of the parameters in the measurement system:

  //  TString filename;
  if ((nsubsamples==0)&&(nbootsamples==0)) {
    filename = TString("meas_par.txt");
  }
  else {
    std::stringstream out;
    out << GetCurrentSample();
    filename = TString("meas_par_").Append(out.str()).Append(".txt");
  }

  filename.Prepend(workingdir);
  if (fexist((char *)filename.Data())) {
    //  if (fexist((char *)"meas_par.txt")) {
    Int_t parnum;
    Double_t lastpar;
    Double_t lastpar_err;
    
    ifstream in1;
    in1.open(filename.Data());
    cout << "Reading in parameters for measurement system:" << filename.Data() << "\n";
    
    for (UInt_t i = 0 ; i < fparam_models[first_outmodel]; i++) {
      in1 >> parnum >> lastpar >> lastpar_err;
      if (!in1.good()) {
	cout << "***ABORTING: Problem reading Initial parameter values\n";
	assert(0);
      }
      //      printf("%4d %8.4f %8.4f\n",parnum,lastpar,lastpar_err);
      
      Setparam(i, lastpar);
      Setparam_err(i, lastpar_err);
    }
    in1.close();
  }
  else {
    cout << "***ABORTING: " << filename.Data() << " not found, need values of parameters in measurement system!!\n";
    assert(0);
  }


  // Read in initial values. If init_par.txt exists, use that.
  filename = "init_par.txt";
  filename.Prepend(workingdir);
  if (fexist(filename.Data())) {
    Int_t parnum;
    Double_t lastpar;
    Double_t lastpar_err;
    
    ifstream in1;
    in1.open(filename.Data());
    cout << "Reading in " << (nparam-fparam_models[first_outmodel]) << " initial values\n";
    
    for (UInt_t ipar = 0 ; ipar < nparam; ipar++) {
      in1 >> parnum >> lastpar >> lastpar_err;
      if (!in1.good()) {
	cout << "Problem reading Initial parameter values\n";
	assert(0);
      }
      //      printf("%4d %8.4f %8.4f\n",parnum,lastpar,lastpar_err);

      if ( ipar >= fparam_models[first_outmodel] ) {
	Setparam(ipar, lastpar*loadingMultiplier);
	Setparam_err(ipar, lastpar_err);
      }	
    }
    in1.close();
  }

  /*
    //temporary code when checkpoint was not printing out all the paramter values
    // Read in initial values. If init_par.txt exists, use that.
  filename = "checkpoint_nohess3.txt";
  filename.Prepend(workingdir);
  if (fexist(filename.Data())) {
    Int_t parnum;
    Double_t lastpar;
    Double_t lastpar_err;
    
    ifstream in1;
    in1.open(filename.Data());
    //cout << "Reading in " << (nparam-fparam_models[first_outmodel]) << " initial values\n";
    cout << "Reading in checkpoint_nohess3.txt\n";

    Int_t nparfree = 0;
    for (UInt_t ipar = fparam_models[first_outmodel] ; ipar < nparam; ipar++) {
      if (GetFixPar(ipar)==0) {
	in1 >> parnum >> lastpar >> lastpar_err;
	if (!in1.good()) {
	  cout << "Problem reading Initial parameter values\n";
	  assert(0);
	}
      //      printf("%4d %8.4f %8.4f\n",parnum,lastpar,lastpar_err);
	nparfree++;
	Setparam(ipar, lastpar);
	Setparam_err(ipar, lastpar_err);
      }
    }
    cout << "Done reading in " << nparfree << " parameters from checkpoint_hess3.txt\n";
    in1.close();
  }
  */
  
  Double_t fvalue = 0;
  Double_t * thisparam = new Double_t[nparam];
  Double_t * grad = new Double_t[nparam];


  if (first_outmodel!=-1) {
    for (UInt_t imod = first_outmodel ; imod < models.size() ; imod++) SetIgnoreMod(imod);
    for (UInt_t i = fparam_models[first_outmodel] ; i < nparam ; i++) FixPar(i);
  }

  Int_t ifreepar = 0;
  for (UInt_t ipar = 0 ; ipar < nparam ;ipar++) {
    if (parfixed[ipar]==0) {
      thisparam[ifreepar] = param[ipar];
      ifreepar++;
    }
  }

  if (printlvl>2) {
    cout << endl << "********Parameter values:" << endl;
    PrintParam(0);
  }
  cout << "Calculating likelihood of measurement system..." << endl;

  TStopwatch timer;
  timer.Reset();
  timer.Start();
  LkhdFcn(ifreepar,grad,fvalue,thisparam,1);
  timer.Stop();
  cout << "Measurement system Likelihood=" << fvalue << endl;
  cout << "Calculating the Likelihood took " << timer.CpuTime() << " seconds.\n";
  
  for (UInt_t imod = 0 ; imod < models.size() ; imod++) RemoveIgnoreMod(imod);
  for (UInt_t i = 0 ; i < nparam ; i++) ReleasePar(i);

  //Run regressions first:
  Int_t ierflg = 0;

  if (initBetaOLS==1) {
    initializing=1;
    for (UInt_t imod = 0 ; imod < models.size() ; imod++) SetIgnoreMod(imod);
    cout << "**************Running individual regressions" << "\n";
    for (UInt_t imod = first_outmodel ; imod < models.size() ; imod++) {
      if (imod!=type_model) {
	if ((printlvl>0)||(imod%10==0)||(imod == UInt_t(first_outmodel))||(imod+1 == models.size())) cout << "****Minimizing Model #" << imod << "\n";
	if (printlvl>1) PrintParam(imod+1);
	RemoveIgnoreMod(imod);
	for (UInt_t i = 0 ; i < nparam ; i++) {
	  if ((i<fparam_models[imod]) || (i>=fparam_models[imod]+models[imod].GetNreg())) FixPar(i);
	  else ReleasePar(i);
	}
	if (models[imod].GetType()==1) ReleasePar(fparam_models[imod]+nparam_models[imod]-1);
	
	if ( (models[imod].GetType()==3) && (models[imod].GetNchoice()>2) ) {
	  for (int ichoice = 1 ; ichoice < models[imod].GetNchoice()-1 ; ichoice++) {
	    for (Int_t i =  0; i <  models[imod].GetNreg() ; i++) {
	    ReleasePar((fparam_models[imod] + ichoice*(models[imod].GetNreg()+nfac+ntyp*(imod!=type_model))) + i);
	    }
	  }
	}
	if (models[imod].GetType()==4) {
	  for (Int_t i = 0 ; i < models[imod].GetNchoice()-1 ; i++) {
	    ReleasePar(fparam_models[imod]+nparam_models[imod]-1-i); 
	  }
	}

	// Only estimate factor loadings for factors that load in the measurment system
	// In other words, do not estimate random effects
	if  (initEstOutcomeLoadings) {
	  cout << "***Estimating loadings using factor scores.\n";
	  for (UInt_t ifac = 0 ; ifac < nfac; ifac++) {
	    //	  if ((norm_models[ifac]>-2)||(fstderr[0][ifac]<0.0)) {
	    if (fstderr[0][ifac]>0.0) {
	      ReleasePar(fparam_models[imod] + models[imod].GetNreg() + ifac); 
	      
	      if ( (models[imod].GetType()==3) && (models[imod].GetNchoice()>2) ) {
		for (int ichoice = 1 ; ichoice < models[imod].GetNchoice()-1 ; ichoice++) {
		  ReleasePar(fparam_models[imod] + ichoice*(models[imod].GetNreg()+nfac+ntyp*(imod!=type_model)) + models[imod].GetNreg() + ifac);
		}
	      }
	    }
	  }
	}

      
	//    ierflg = Min_Minuit(-1);
	ierflg = Min_Ipopt(1);
	if (printlvl>0) PrintParam(imod+1);
	if ((ierflg!=0)&&(nsubsamples==0)&&(nbootsamples==0)) assert(0);

	SetIgnoreMod(imod);
      }
    }
  }
  

  // Now minimize each model groups with the measurement system and factors
  initializing=0;
  for (Int_t imod = 0 ; imod < first_outmodel ; imod++) RemoveIgnoreMod(imod);
  //     for (UInt_t imod = 0 ; imod < models.size() ; imod++) RemoveIgnoreMod(imod);
  //     for (UInt_t i = 0 ; i < nparam ; i++) ReleasePar(i); //parfixed[i] = 0;
  
  
  // Ignore all outcome models and fix all parameters
  for (UInt_t imod = first_outmodel; imod < models.size() ; imod++) SetIgnoreMod(imod);
  for (UInt_t i = 0 ; i < nparam ; i++) FixPar(i);


  vector<Int_t> failed_groups;

  //  TStopwatch timer;
  timer.Reset();
  timer.Start();
  // Estimate Outcomes in each group:
  for (Int_t igroup = 1 ; igroup <= last_group; igroup++) {
    for (UInt_t imod = first_outmodel; imod < models.size() ; imod++) SetIgnoreMod(imod);
    for (UInt_t i = 0 ; i < nparam ; i++) FixPar(i);
    // Now activate all models in this group and release relevant parameters:
    cout << "\n\n*****************Estimating Group #" << igroup << "\n";

    for (UInt_t imod = first_outmodel ; imod < models.size() ; imod++) {
      if (models[imod].GetGroup()==igroup) {
	RemoveIgnoreMod(imod);
	for (UInt_t ipar = fparam_models[imod] ; ipar < fparam_models[imod]+nparam_models[imod] ; ipar++) ReleasePar(ipar);
	if (printlvl>0) PrintParam(imod+1);
      }
    }

    //    ierflg = Min_Minuit(printlevel);
    ierflg = Min_Ipopt(printlevel); 
    //    if (ierflg!=0) assert(0);
    if (ierflg!=0) { 
      failed_groups.push_back(igroup);
      EstGroup_lkhd.push_back(-9999.0);
    }
    else {  
      ifreepar = 0;
      for (UInt_t ipar = 0 ; ipar < nparam ;ipar++) {
	if (parfixed[ipar]==0) {
	  thisparam[ifreepar] = param[ipar];
	  ifreepar++;
	}
      }
      LkhdFcn(ifreepar,grad,fvalue,thisparam,1);
      cout << "Final Likelihood=" << fvalue << endl;
      EstGroup_lkhd.push_back(fvalue);
    }

    for (UInt_t imod = first_outmodel ; imod < models.size() ; imod++) {
      if ((printlvl>0) &&(models[imod].GetGroup()==igroup)) PrintParam(imod+1);
    }

  }
  timer.Stop();

  //Reset everything
  initializing=0;
  for (UInt_t imod = 0 ; imod < models.size() ; imod++) RemoveIgnoreMod(imod);
  for (UInt_t i = 0 ; i < nparam ; i++) ReleasePar(i);

  delete [] thisparam;
  delete [] grad;

  cout << endl << "********Parameter values:" << endl;
  PrintParam();
  cout << "Minimization took " << timer.CpuTime() << " s." << endl;

  if (failed_groups.size()>0) {
    cout << "But some groups failed to converge: ";
    for (UInt_t i = 0 ; i < failed_groups.size() ; i++) cout << ", " << failed_groups[i];
    cout << "\n";
  }

  for (UInt_t i = 0 ; i < EstGroup_lkhd.size() ; i++) printf("Lkhd(group=%d)=%f\n",i,EstGroup_lkhd.at(i));

  FILE * pFile;

  if ((nsubsamples==0)&&(nbootsamples==0)) {
    filename = TString("fullmodel_par.txt");
  }
  else {
    std::stringstream out;
    out << GetCurrentSample();
    filename = TString("fullmodel_par_").Append(out.str()).Append(".txt");
  }

  filename.Prepend(workingdir);
  pFile = fopen ( ((char *)filename.Data()),"w");
  for (UInt_t ipar=0 ; ipar< nparam ; ipar++)
    {
      fprintf (pFile, "%5d %12.8f %12.8f\n",ipar,Getparam(ipar),Getparam_err(ipar));
    }
  fclose (pFile);

  return 0;
}

Int_t TMinLkhd::Simulate(Int_t printlevel) {

  printf("Starting TMinLkhd::Simulate()\n");
  if (simWithData==0) {
  //Set up endogenous regressors in each model:
  for (UInt_t imod = 0 ; imod < models.size() ; imod++) {
    models.at(imod).SetEndogenousRegs(models,var_table);
  }
  printf("Finished setting endogenous variables!\n");
  }
  
  TString filename;
  if ((nsubsamples==0)&&(nbootsamples==0)) {
    filename = TString("fullmodel_par.txt");
  }
  else {
    std::stringstream out;
    out << GetCurrentSample();
    filename = TString("fullmodel_par_").Append(out.str()).Append(".txt");
  }

  filename.Prepend(workingdir);
  if (fexist((char *)filename.Data())) {
  //  if (fexist((char *)"fullmodel_par.txt")) {
    Int_t parnum;
    Double_t lastpar;
    Double_t lastpar_err;
    
    ifstream in1;
    //    in1.open("fullmodel_par.txt");
    in1.open((char *)filename.Data());
    cout << "Reading in " << nparam << " parameters for full model\n";
    
    for (UInt_t i = 0 ; i < nparam; i++) {
      in1 >> parnum >> lastpar >> lastpar_err;
      //      in1 >> parnum >> lastpar;
      if (!in1.good()) {
	cout << "***ABORTING: Problem reading in parameter values\n";
	assert(0);
      }
      //      printf("%4d %8.4f %8.4f\n",parnum,lastpar,lastpar_err);
      
      Setparam(i, lastpar);
      Setparam_err(i, lastpar_err);
    }
    in1.close();
  }
  else {
    cout << "***ABORTING: fullmodel_par.txt not found, need values of parameters to simulate!!\n";
    return 1;
  }

  PrintParam();
  if ((nsubsamples==0)&&(nbootsamples==0)) {
    PrintParamTab();
    PrintParam_Varlist();
  }
  TRandom *r3 = new TRandom3();
  gRandom = r3;
  
  //calculate the weights and means of mixtures
  
  vector <Double_t> w_mix;
  vector <Double_t> fac_mean;
  w_mix.resize(fac_nmix);

  if (nfac>0) {
    
    //Calculate the denominator for the mixture weights
    Double_t w_mix_den = 1.0;
    if (fac_nmix>1) {
      for (UInt_t imix = 0 ; imix < fac_nmix-1 ; imix++) {
	w_mix_den += exp(Getfw(imix));
      }
    }
    
    // Calculate actual weights and means for different mixtures:
    for (UInt_t imix = 0 ; imix < fac_nmix*(nfac>0) ; imix++) {
      //Calculate weight for this mixture (logit-like scheme)
      w_mix[imix] = 1.0;
      if ((fac_nmix>1)&&(nfac>0)) { 
	w_mix[imix] = 1.0/w_mix_den;
	if (imix<fac_nmix-1) w_mix[imix] = exp(Getfw(imix))*w_mix[imix];
      }
      
      //Get means for this mixture
      if (fac_nmix==1) {
	for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) fac_mean.push_back(0.0);
      }
      else if (imix < fac_nmix-1) {
	for (UInt_t ifac = 0 ; ifac < nfac; ifac++) fac_mean.push_back(Getfmean(imix,ifac));
      }
      else {
	// Calculate mean of third mixture such that mean of all mixtures = 0.0:
	for (UInt_t ifac = 0 ; ifac < nfac; ifac++) {
	  Double_t weighted_sum = 0.0;
	  for (UInt_t imix_mean = 0 ; imix_mean < fac_nmix-1 ; imix_mean++) {
	    weighted_sum += -exp(Getfw(imix_mean))*Getfmean(imix_mean,ifac);
	  }
	  fac_mean.push_back(weighted_sum);
	}
      }
    }
    
    cout << "finished calculating means and weights:\n";
    cout << "fac_mean size=" << fac_mean.size() << " and w_mix size=" << w_mix.size() << endl;
    for (UInt_t imix = 0 ; imix < fac_nmix*(nfac>0) ; imix++) {
      cout << "imix=" << imix << ", mu1=";
      
      for (UInt_t ifac = 0 ; ifac < nfac ; ifac++)
	cout << fac_mean.at(imix*nfac+ifac); 
      
      cout << ", weight=" << w_mix[imix] << ".\n";
    }

  } //nfac>0

  // open file here
  FILE * pFile;

  if ((nsubsamples==0)&&(nbootsamples==0)) {
    filename = TString("simulation_data.csv");
  }
  else {
    std::stringstream out;
    out << GetCurrentSample();
    filename = TString("simulation_data_").Append(out.str()).Append(".csv");
  }
  if (simWithData==1) {
    filename = TString("simulation_data_gof.csv");
  }
  
  filename.Prepend(workingdir);
  pFile = fopen ( ((char *)filename.Data()),"w");

  //Write dictionary part here
  //  fprintf(pFile,"dictionary{\n");


  if (indexvar==-1) {
    fprintf(pFile,"xid");
  }
  else {
    fprintf(pFile,"%s",(var_table.at(indexvar)).Data());
  }
  vector<Int_t> varuse;
  varuse.resize(nvar,0);
  Int_t nsimvar = nfac + ntyp;

  for (UInt_t imod = 0 ; imod < models.size() ; imod++) {

    UInt_t thismod = imod;
    if (sim_modelorder.size()==models.size()) thismod = sim_modelorder.at(imod);  
    else if (type_model<9999) {
      // here we rearrange the order so type_model goes first
      if (imod==0) thismod = type_model;
      else thismod = imod-1*(imod<=type_model);
    }

    //    printf("Writing header for model %d\n",thismod);
    
    // simulation gives: missing, I_obs, I_factor, shock/prob, outcome
    fprintf(pFile,", %s",(TString("sim_").Append(TString(models[thismod].GetName()).Append("_miss"))).Data());
    nsimvar++;

    if (models[thismod].GetDetailSim()) {
      int nlogitchoice = 2;
      if (models[thismod].GetType()==3) nlogitchoice = models[thismod].GetNchoice();

      // If it is a multinomial logit model we need to produce nchoice-1 * V and shprob variables
      for (int ichoice = 2 ; ichoice <= nlogitchoice ; ichoice++) {
      
	// Only want this more than once if it is a multinomial logit model
	if ( (models[thismod].GetType()==3) || (ichoice==2) ) {
	  fprintf(pFile,", %s",(TString("sim_").Append(TString(models[thismod].GetName()).Append("_Vobs"))).Data());
	  if ((models[thismod].GetType()==3) && (nlogitchoice>2)) fprintf(pFile,"_C%1d",ichoice);
	  nsimvar +=1;

	  if ( (models[thismod].GetEndogenousReg().size()>0)&&(simWithData==0) ) {
	    fprintf(pFile,", %s",(TString("sim_").Append(TString(models[thismod].GetName()).Append("_Vend"))).Data());
	    if ((models[thismod].GetType()==3) && (nlogitchoice>2)) fprintf(pFile,"_C%1d",ichoice);
	    nsimvar +=1;
	  }
	
	  for (UInt_t ifac = 0; ifac < nfac ; ifac++) {
	    fprintf(pFile,", %s%d",(TString("sim_").Append(TString(models[thismod].GetName()).Append("_Vfac"))).Data(),ifac);
	    if ((models[thismod].GetType()==3) && (nlogitchoice>2)) fprintf(pFile,"_C%1d",ichoice);
	    nsimvar +=1;
	  }

	  if ((ntyp>0)&&(thismod!=type_model)) {
	    for (UInt_t ityp = 0; ityp < ntyp ; ityp++) {
	      fprintf(pFile,", %s%d",(TString("sim_").Append(TString(models[thismod].GetName()).Append("_Vtyp"))).Data(),ityp);
	      if ((models[thismod].GetType()==3) && (nlogitchoice>2)) fprintf(pFile,"_C%1d",ichoice);
	      nsimvar +=1;
	    }
	  }
	}
      }

      // If it is a multinomial logit model we need to produce nchoice - 1 eps variables
      for (int ichoice = 2 ; ichoice <= nlogitchoice ; ichoice++) {
	// Only want this more than once if it is a multinomial logit model
	if ( (models[thismod].GetType()==3) || (ichoice==2) ) {
	  fprintf(pFile,", %s",(TString("sim_").Append(TString(models[thismod].GetName()).Append("_eps"))).Data());
	  if ((models[thismod].GetType()==3) && (nlogitchoice>2)) fprintf(pFile,"_C%1d",ichoice);
	  nsimvar +=1;
	}
      }
    }

    // If it is a discrete choice model we need to produce nchoice-1 prob variables
    if ((models[thismod].GetType()==2) || (models[thismod].GetType()==3) || (models[thismod].GetType()==4)) {

      int nchoice = 2;
      if ((models[thismod].GetType()==3) || (models[thismod].GetType()==4)) nchoice = models[thismod].GetNchoice();

      for (int ichoice = 2 ; ichoice <= nchoice ; ichoice++) {
	// generate probabilities for choices in multinomial models
	fprintf(pFile,", %s",(TString("sim_").Append(TString(models[thismod].GetName()).Append("_prob"))).Data());
	if ( ((models[thismod].GetType()==3) || (models[thismod].GetType()==4)) && (nchoice>2) ) fprintf(pFile,"_C%1d",ichoice);
	nsimvar += 1;
      }
    }

    fprintf(pFile,", %s",(TString("sim_").Append(TString(models[thismod].GetName()))).Data());
    nsimvar++;

    //This is where we split up a model using an indicator variable
    if (models[thismod].GetSplitSim()) {
      fprintf(pFile,", %s",(TString("sim_").Append(TString(models[thismod].GetName()).Append("1_miss"))).Data());
      if (models[thismod].GetDetailSim()) {
	fprintf(pFile,", %s",(TString("sim_").Append(TString(models[thismod].GetName()).Append("1_Vobs"))).Data());
	nsimvar +=1;
	
	if (models[thismod].GetEndogenousReg().size()>0) {
	  fprintf(pFile,", %s",(TString("sim_").Append(TString(models[thismod].GetName()).Append("1_Vend"))).Data());
	  nsimvar +=1;
	}

	for (UInt_t ifac = 0; ifac < nfac ; ifac++) {
	  fprintf(pFile,", %s%d",(TString("sim_").Append(TString(models[thismod].GetName()).Append("1_Vfac"))).Data(),ifac);
	  nsimvar +=1;
	}

	if (thismod!=type_model) {
	  for (UInt_t ityp = 0; ityp < ntyp ; ityp++) {
	    fprintf(pFile,", %s%d",(TString("sim_").Append(TString(models[thismod].GetName()).Append("1_Vtyp"))).Data(),ityp);
	    nsimvar +=1;
	  }
	}
	
	fprintf(pFile,", %s",(TString("sim_").Append(TString(models[thismod].GetName()).Append("1_eps"))).Data());
	nsimvar +=1;

      }
      
      if ((models[thismod].GetType()==2)||(models[thismod].GetDetailSim())) {
	fprintf(pFile,", %s",(TString("sim_").Append(TString(models[thismod].GetName()).Append("1_prob"))).Data());
	nsimvar++;
      }
      
      fprintf(pFile,", %s",(TString("sim_").Append(TString(models[thismod].GetName()).Append("1"))).Data());
      nsimvar +=2;
    } // Get splitsim
    
    // figure out which variables are used:
    vector <Int_t> regs = models[thismod].GetReg();
    for (UInt_t ireg = 0 ; ireg < regs.size(); ireg++) varuse.at(regs[ireg]) = 1;
    if (thismod!=type_model) {
      varuse.at(models[thismod].GetOutcome()) = 1;
    }
    if (models[thismod].GetMissing()>-1) varuse.at(models[thismod].GetMissing()) = 1;
  }

  //Save factor and type values
  for (UInt_t ifac = 0; ifac < nfac ; ifac++) {
    fprintf(pFile,", sim_fac%d",ifac);
  }
  for (UInt_t itype = 1; itype <= ntyp ; itype++) {
    fprintf(pFile,", sim_type%d",itype);
  }

  
  Int_t nvarused = 0;
  if (simIncData) {
    for (UInt_t ivar = 0 ; ivar < nvar ; ivar++) {
      if (varuse.at(ivar) ==1) {
	fprintf(pFile,", %s",(TString("dt_").Append(var_table[ivar])).Data());
	nvarused++;
      }
    }
  }

  //  fprintf(pFile,"sim_flag }\n");
  fprintf(pFile,", sim_flag\n");
  //Finished writing dictionary
  
  cout << "Starting simulation file will have " << (1+nvarused+nsimvar+1) << " variables.\n";

  vector<Double_t> maxlkhd(nobs,0.0);
  Double_t fvalue = 0;
  Double_t * thisparam = new Double_t[nfac];
  Double_t * grad = new Double_t[nfac];
  Double_t * hess = new Double_t[(nfac*(nfac+1))/2];
  Int_t npar_min = nfac;


  // Get maxlkhd of measurement system for each observation
  if (sampleposterior==1) {
    // setup calculation of posterior
    Int_t first_outmodel=-1;
    for (UInt_t imod = 0 ; imod < models.size() ; imod++) {
      if ((models[imod].GetGroup()==1)&&(first_outmodel==-1)) first_outmodel = imod;
    }
    
    //Ignore all outcome models when predicting
    if (first_outmodel!=-1) {
      for (UInt_t imod = first_outmodel ; imod < models.size() ; imod++) SetIgnoreMod(imod);
    }
    if (first_outmodel==-1) first_outmodel = models.size();
    
    // Fix all parameters of models
    for (UInt_t i = 0 ; i < nparam ; i++) FixPar(i);
    
    predicting = 1;
    predictobs = 0;

//     // clear file of predicted factor scores
//     if ((nsubsamples==0)&&(nbootsamples==0)) {
//       FILE * pFile;
//       filename = TString("factor_predictions.txt");
//       filename.Prepend(workingdir);
//       pFile = fopen (filename.Data(),"w");
//       fclose (pFile);
//     }

    Int_t ierflg = 0;
    // get posterior density at mode
    for (UInt_t iobs = 0 ; iobs < nobs ; iobs++) {
      //       if(iobs%100==0) printf("Predicting factors for observation #%d\n",iobs);
      //      predictobs = iobs;
      //      ierflg = Min_Ipopt(0);

      for (UInt_t ifac = 0; ifac <  nfac; ifac++) {
	thisparam[ifac] = fscore[iobs][ifac];
      }
      predictobs = iobs;
      LkhdFcn(npar_min,grad,fvalue,thisparam,1,hess);
      maxlkhd.at(iobs) = fvalue;
      //      printf("*** maxlkhd(obs=%d)=%5.3f\n\n",iobs,fvalue);

      if (ierflg!=0) {
	printf("PREDICT_FACTOR: Failed to minimize factor for observation #%d\n",iobs);
        assert(0);
      }
    }
    printf("Finished Calculating value of likelihood at factor scores\n");
  }

  // Get maxweight
  Double_t maxweight = -1.0;
  if (weightvar!=-1) {
    for (UInt_t iobs = 0 ; iobs < nobs ; iobs++) {
      if (data[iobs*nvar+weightvar]>maxweight) maxweight = data[iobs*nvar+weightvar];
    }
  }

  //This is mostly for goodness of fit where we simulate each model using data and factor scores
  if (simWithData==1) sim_nobs = nobs;
  for (UInt_t igen = 0 ; igen < sim_nobs ; igen++) {
    if (igen%10000==0) printf("Simulating %5d\n",igen);
    //    printf("Simulating %5d: ",igen);
    

    // First draw observation:
    UInt_t obs_drw = -1;
    int acceptdraw = 0;
    UInt_t ndraws = 0;

    if (simWithData==0) {
      while (acceptdraw==0) {
	ndraws++;
	
	obs_drw = UInt_t(r3->Rndm()*(nobs));
	if (obs_drw==nobs) obs_drw--;
	
	if (weightvar==-1) acceptdraw=1;
	else if (r3->Uniform() < data[obs_drw*nvar+weightvar]/maxweight) acceptdraw=1;
      }
    }
    else obs_drw = igen;

    //    if (igen%1000==0)&&(weightvar!=-1) printf("Drawing Xs, after %d draws, weighratio=%5.3f.\n",ndraws,data[obs_drw*nvar+weightvar]/maxweight);

    vector<Double_t> fac_val;

    if (nfac>0) {

      if (simWithData==0) {

	acceptdraw = 0;
	ndraws = 0;
	while (acceptdraw==0) {
	  ndraws++;
	  // First draw factors from unconditional factor distribution
	  vector<Double_t> f_draw;
	  
	  // Draw from Normal distribution
	  for (UInt_t i = 0 ; i < nfac; i++) {
	    f_draw.push_back(r3->Gaus(0.0,1.0));
	    //draw from one of two points if only two quadrature points are used
	    if ((fac_npoints.size()==nfac) && (fac_npoints.at(i)==2)) f_draw[i] = GetHGQx(2, (f_draw[i]<0) );
	  }

	  //Select mixture
	  UInt_t imix = 0;
	  if (fac_nmix>1) {
	    Double_t draw = r3->Uniform();
	    if (draw<w_mix.at(0)) imix = 0;
	    else if (draw < w_mix.at(0) + w_mix.at(1)) imix=1;
	    else imix=2;
	  }
	  
	  fac_val.clear();
	  
	  //Calculate factors 
	  if ((nfac==2)&&(fac_corr!=0)) {
	    fac_val.push_back(fac_mean.at(imix*nfac) + Getfvar(imix,0)*f_draw.at(0));
	    fac_val.push_back(fac_mean.at(imix*nfac+1) + Getfvar(imix,1)*(f_draw.at(0)*Getfvar(imix,2) + f_draw.at(1)*sqrt(1-Getfvar(imix,2)*Getfvar(imix,2)))); 
	  }
	  else {
	    for (UInt_t ifac = 0 ; ifac < nfac; ifac++) fac_val.push_back(fac_mean.at(imix*nfac+ifac) + Getfvar(imix,ifac)*f_draw.at(ifac));
	  }
	  
	  // check if we are sampling from posterior or just factor distribution
	  if (sampleposterior==0) acceptdraw=1;
	  else {
	    predictobs = obs_drw;
	    
	    // get posterior density 
	    
	    for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) thisparam[ifac] = fac_val.at(ifac);
	    
	    LkhdFcn(npar_min,grad,fvalue,thisparam,1,hess);
	    //	  fvalue = 0.0;
	    Double_t postdensity = exp(maxlkhd.at(obs_drw)-fvalue);
	    Double_t draw = r3->Uniform();
	    
	    //	  printf("%8d ",obs_drw);
	    //	  for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) printf("%7.3f ",fac_val.at(ifac));
	    //	  printf("lkhd=%5.3f maxlkhd=%5.3f prob=%5.4f draw=%5.4f\n", fvalue, maxlkhd.at(obs_drw), postdensity, draw);
	    
	    if (draw<postdensity) acceptdraw=1;
	    //	  if (1) acceptdraw=1;
	    else {
	      // reject draw
	      if (ndraws>10000) {
		printf("***Simulating: Couldn't find factors that pass acceptance sampling! \n Using factors that maximize the likelihood of measurement system for observation %d\n",obs_drw);
		for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {
		  //		printf("Setting %d factor for %d observation where arrays are size %d and %d.\n",ifac,obs_drw,fac_val.size(),fscore.size());
		  fac_val.at(ifac) = fscore[obs_drw][ifac];  
		  acceptdraw=1;
		}
	      }
	    
	      // if (0) {
	      //   Int_t ierflg = 0;
	      //   ierflg = Min_Ipopt(0);
	      
	      //   if (ierflg!=0) {
	      // 	printf("PREDICT_FACTOR: Failed to minimize factor for observation #%d\n",iobs);
	      // 	assert(0);
	      //   }
	      
	      //   printf("***Simulating: Couldn't find factors that pass acceptance sampling! Using factors that maximize the likelihood for observation %d, ",obs_drw);
	      //   for (UInt_t ifac = 0; ifac <  nfac; ifac++) {
	      // 	fac_val.at(ifac) = facprediction.at(ifac);
	      // 	printf("fac%d = %8.3f,  ", ifac, fac_val.at(ifac));
	      //   }
	      //   printf("\n");
	      //   acceptdraw=1;
	      //	    }
	    }
	  }
	} // while (!acceptdraw)
      //      if (igen%1000==0) printf("Drew Theta(%d), after %d draws.\n",obs_drw,ndraws);
      }
      else {
	fac_val.clear();
	for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) fac_val.push_back(fscore[obs_drw][ifac]); 
      }
    }

    if (indexvar==-1) {
      fprintf(pFile,"%10d",obs_drw+1);
    }
    else {
      fprintf(pFile,"%10d", Int_t(data[obs_drw*nvar+indexvar]));
    }

    //Now draw the type
    //    char * typemodeldata;
    if (type_model<9999) {
      // typemodeldata = new char[1200];
      //Need to save output of type model to write to simulation file later
      //      std::FILE* tmpf = fopen("typetmpfile.txt","w");
	//std::tmpfile();
      //      printf("Simulating type model\n");
      //      printf("%5d, ",type_model);
      models[type_model].Sim(obs_drw*nvar,data,models,param,fparam_models[type_model], fac_val, pFile, simWithData);

      //      printf("DONE Simulating type model\n");
      //      std::rewind(tmpf);
//      std::fgets(typemodeldata, sizeof(typemodeldata), tmpf);
//      std::fclose(tmpf);
//      printf("saved temporary info! %s\n",typemodeldata);
      
      //Save type to simulation file
      Int_t typedraw = Int_t(models[type_model].GetSimResult());
      //      printf("The type is %d\n",typedraw);
      for (UInt_t itype = 0 ; itype < ntyp; itype++) {
	//Types are drawn starting from 1, which is the omitted category
	if (Int_t(itype+2)==typedraw) {
	  //	  fprintf(pFile,", %10d",1);
	  fac_val.push_back(1.0);
	}
	else {
	  //	  fprintf(pFile,", %10d",0);
	  fac_val.push_back(0.0);
	}
      }
    }

    //    printf("Starting loop simulating models\n");
    // Now simulate the models:
    for (UInt_t imod = 0 ; imod < models.size() ; imod++) {
      UInt_t thismod = imod;
      if (sim_modelorder.size()==models.size()) thismod = sim_modelorder.at(imod);

      if (thismod != type_model) {
	//	printf("%5d, ",thismod);
	      
	// simulation gives: missing, I_obs, I_factor, shock/prob, outcome
	models[thismod].Sim(obs_drw*nvar,data,models,param,fparam_models[thismod],fac_val, pFile, simWithData);	
      }
//      else {
//	// Finish copying simulation output of type model from tmp file to simulation file
//	if (sizeof(typemodeldata)>0) {
//	  int line_length = std::strlen(typemodeldata);
//          TString format = (TString("\"%").Append(line_length)).Append("s\"");
//	  printf("Printing out model data length %d: %s\n",line_length, typemodeldata);
//	  std::fprintf(pFile, (char *)format.Data(), typemodeldata);
//	  delete typemodeldata;
//	}
//	else {
//	  cout << "ERROR (TMinLkhd::Simulate): typmodeldata is zero length!" << endl;
//	  assert(0);
//	}
//    }
    }

    for (UInt_t ifac = 0 ; ifac < nfac + ntyp; ifac++) fprintf(pFile,", %10.5f",fac_val.at(ifac));


    if (simIncData) for (int ivar = 0 ; ivar < nvarused ; ivar++) fprintf(pFile,", %10.5f",-9999.0);
    fprintf(pFile,", %10d \n",1);

    // Clear outcome data so they aren't used for other observations
    for (UInt_t imod = 0 ; imod < models.size() ; imod++) {
      models.at(imod).ClearSimResult();
    }
    //    printf("\n");
  }

  // Now append data to file
  if (simIncData) {
    for (UInt_t iobs = 0 ; iobs < nobs ; iobs++) {
      //fill in -9999 for sim variables
      fprintf(pFile,", %10d",iobs+1);
      for (Int_t isimvar = 0 ; isimvar < nsimvar ; isimvar++) fprintf(pFile,", %10.5f",-9999.0);
      // Now data
      for (UInt_t ivar = 0 ; ivar < nvar ; ivar++) {
	if (varuse.at(ivar) ==1) fprintf(pFile,", %10.5f",data[iobs*nvar+ivar]);
      }
      fprintf(pFile,", %10d\n",0);
    }
  }

  fclose (pFile);


  delete [] thisparam;
  delete [] grad;
  delete [] hess;
  
  return 0;
}

Int_t TMinLkhd::PredictFactors(Int_t printlevel) {

  printf("Predicting factors...\n");

  if (type_model<9999) {
    cout << "ERROR (TMinLkhd::PredictFactors): Not currently possible to predict factors when there are unobserved types." << endl;
    assert(0);
  }
  
  TString filename;
  if ((nsubsamples==0)&&(nbootsamples==0)) {
    filename = TString("meas_par.txt");
  }
  else {
    std::stringstream out;
    out << GetCurrentSample();
    filename = TString("meas_par_").Append(out.str()).Append(".txt");
  }

  filename.Prepend(workingdir);
  if (fexist((char *)filename.Data())) {
  //  if (fexist((char *)"fullmodel_par.txt")) {
    Int_t parnum;
    Double_t lastpar;
    Double_t lastpar_err;
    
    ifstream in1;
    //    in1.open("fullmodel_par.txt");
    in1.open((char *)filename.Data());
    cout << "Reading in " << nparam << " parameters for full model\n";
    
    for (UInt_t i = 0 ; i < nparam; i++) {
      in1 >> parnum >> lastpar >> lastpar_err;
      //      in1 >> parnum >> lastpar;
      if (!in1.good()) {
	cout << "***ABORTING: Problem reading in parameter values\n";
	assert(0);
      }
      //      printf("%4d %8.4f %8.4f\n",parnum,lastpar,lastpar_err);
      
      Setparam(i, lastpar);
      Setparam_err(i, lastpar_err);
    }
    in1.close();
  }
  else {
    cout << "***ABORTING: meas_par.txt not found, need values of parameters to simulate!!\n";
    return 1;
  }

  if (includePriorFactorScore==1) {
    if (fac_nmix==1) {
      Double_t normfac[10] = {0.0};
      
      for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {
	normfac[ifac] = 1.0;
	vector<Int_t> intdata(2,-1);
	TModel newmodel("Prior","Factor Prior", 1, 0,
			current_printgroup,intdata,nfac,ntyp,normfac);
	models.push_back(newmodel);
	nparam_models.push_back(1);
	fparam_models.push_back(nparam);
	nobs_models.push_back(nobs);
	nparam++;
	normfac[ifac] = 0.0;

	param.resize(nparam,param[ifac]);
	param_err.resize(nparam,param_err[ifac]);
	parfixed.resize(nparam,0);
	printf("Finished adding Prior #%d.\n",ifac);
      }
      parconstrained.resize(nparam,-1);
      param_fixval.resize(nparam,-9999.0);
    }
    else {
      cout << "ERROR TMinLkhd::PredictFactors: Cannot include prior if using mixtures of normals!!"
	   << endl;
      assert(0);
    }
  }

  PrintModels();
  PrintParam();


  // Find first outcome model
  Int_t first_outmodel=-1;
  for (UInt_t imod = 0 ; imod < models.size() ; imod++) {
    if ((models[imod].GetGroup()==1)&&(first_outmodel==-1)) first_outmodel = imod;
  }

  //Ignore all outcome models when predicting
  if (first_outmodel!=-1) {
    for (UInt_t imod = first_outmodel ; imod < models.size() ; imod++) SetIgnoreMod(imod);
  }
  if (first_outmodel==-1) first_outmodel = models.size();

  // Fix all parameters of models
  for (UInt_t i = 0 ; i < nparam ; i++) FixPar(i);
  
  predicting = 1;

  predictobs = 0;

//   Int_t ierflg = 0;
//   ierflg = Min_Ipopt(1);

//   TestGradient();

//   TestHessian();


//   // clear file of predictions
  FILE * pFile;

  //  TString filename;
  if ((nsubsamples==0)&&(nbootsamples==0)) {
    filename = TString("factor_predictions.txt");
  }
  else {
    std::stringstream out;
    out << GetCurrentSample();
    filename = TString("factor_predictions_").Append(out.str()).Append(".txt");
  }

  //  filename = TString("factor_predictions.txt");
  filename.Prepend(workingdir);
  pFile = fopen (filename.Data(),"w");
  fclose (pFile);
  
  Int_t ierflg = 0;
  for (UInt_t iobs = 0 ; iobs < nobs ; iobs++) {
    if(iobs%100==0) printf("Predicting factors for observation #%d\n",iobs);
    predictobs = iobs;
    ierflg = Min_Ipopt(0);
    if (ierflg!=0) {
      printf("PREDICT_FACTOR: Failed to minimize factor for observation #%d\n",iobs);
        assert(0);
    }
  }

  predicting = 0;

  return 0;
}


Int_t TMinLkhd::GenerateCHsamples(Int_t printlevel) {
 
  if (nchsamples==0) {
    cout << "TMinLkhd::GenerateCHsamples(): ***ABORTING: Number of samples not defined!\n";
    assert(0);
  }
  
  TString filename;
  filename = TString("fullmodel_par.txt");

  filename.Prepend(workingdir);
  if (fexist((char *)filename.Data())) {
    Int_t parnum;
    Double_t lastpar;
    Double_t lastpar_err;
    
    ifstream in1;
    //    in1.open("fullmodel_par.txt");
    in1.open((char *)filename.Data());
    cout << "Reading in " << nparam << " parameters for full model\n";
    
    for (UInt_t i = 0 ; i < nparam; i++) {
      in1 >> parnum >> lastpar >> lastpar_err;
      if (!in1.good()) {
	cout << "***ABORTING: Problem reading in parameter values\n";
	assert(0);
      }
      //      printf("%4d %8.4f %8.4f\n",parnum,lastpar,lastpar_err);
      
      Setparam(i, lastpar);
      Setparam_err(i, lastpar_err);
    }
    in1.close();
  }
  else {
    cout << "***ABORTING: fullmodel_par.txt not found, need values of parameters to generate CH samples!!\n";
    return 1;
  }

  PrintParam();

  TRandom *r3 = new TRandom3();
  gRandom = r3;
  FILE * pFile;

  for (UInt_t isample = 0 ; isample < nchsamples ; isample++) {
      
    std::stringstream out;
    out << isample;
    filename = TString("fullmodel_par_").Append(out.str()).Append(".txt");
    filename.Prepend(workingdir);
    pFile = fopen ( ((char *)filename.Data()),"w");

    for (UInt_t ipar=0 ; ipar< nparam ; ipar++) {
      Double_t newparam = Getparam(ipar);
      if (parfixed.at(ipar)==0) newparam += r3->Gaus(0.0,1.0)*Getparam_err(ipar);
      fprintf (pFile, "%5d %12.8f %12.8f\n",ipar,newparam,Getparam_err(ipar));
    }
    fclose (pFile);

  }

  return 0;
}




Int_t TMinLkhd::Min_Minuit(Int_t printlevel) {

  Int_t nfixed = 0;
  for (UInt_t i = 0 ; i < nparam; i++) nfixed += (parfixed.at(i)>0);

  cout << "Minimizing " << (nparam-nfixed) << " parameters while fixing " << nfixed << ".\n";
  //  cout << "Super-special counter of free parameters says there are: " << nfreeparam << ".\n";
  Int_t nparam_min = nparam - nfixed;

  //Setup minimizer
  if (gMinuit) delete gMinuit;
  gMinuit = new TMinuit(nparam_min);
  gMinuit->SetObjectFit(this);
  gMinuit->SetPrintLevel(printlevel);
  gMinuit->SetFCN(MinuitLkhdFcn);
  //  gMinuit->SetErrorDef(1.0); // ChiSq
  gMinuit->SetErrorDef(.5);  // ln(likelihood)
  

  // Tell MINUIT about the parameter values:
  Int_t ipar = 0;
  for (UInt_t i = 0 ; i < nparam ; i++) {
    if (parfixed[i]==0) {

     if ((fac_corr!=0) && (i>=nfac) && (i<f_nvariance)) {
       gMinuit->DefineParameter(ipar, "", param[i], 0.01*param_err[i], -0.8 , 0.8);
     }
     // set limits for covariances for second mixture
     else if ((fac_nmix>1) && (fac_corr!=0) && (i>=f_nvariance+nfac+1+nfac) && (i<f_nvariance+nfac+1+f_nvariance)) {
       gMinuit->DefineParameter(ipar, "", param[i], 0.01*param_err[i], -0.8 , 0.8);
     }
     // set limits for covariances for third mixture
     else if ((fac_nmix>2) && (fac_corr!=0) && (i>=2*(f_nvariance+nfac+1)+nfac) && (i<2*(f_nvariance+nfac+1)+f_nvariance)) {
       gMinuit->DefineParameter(ipar, "", param[i], 0.01*param_err[i], -0.8 , 0.8);
     }
     //       else if (i<nfac) {
     // 	gMinuit->DefineParameter(ipar, "", param[i], 0.01*param_err[i], 0.01 , 10.0);
     //       }
     else gMinuit->DefineParameter(ipar, "", param[i], 0.01*param_err[i], 0,0);
     //    without limits:
     //    gMinuit->DefineParameter(ipar, "", param[i], param_err[i], 0,0);
      ipar++;
    }
  }
   
  // Set up other stuff
  Double_t arglist[16];
  Int_t ierflg = 0;
  //arglist[0] = 1;
  //  gMinuit->mnexcm("CALL FCN", arglist, 1, ierflg);
  //  arglist[0] = 2;
  arglist[0] = 2;
  gMinuit->mnexcm("SET STR", arglist, 1, ierflg);
  //arglist[0] = 2.5e-8;
  //gMinuit->mnexcm("SET EPS", arglist, 1, ierflg);
  //  if (printlevel > 0) 
  gMinuit->mnexcm("SHO EPS", arglist, 0, ierflg);

  arglist[0] = 1;
  gMinuit->mnexcm("SET GRA", arglist, 1, ierflg);
  //  gMinuit->mnexcm("SET NOGRA", arglist, 1, ierflg);
  // For checking the gradient:
  //  gMinuit->mnexcm("SET GRA", arglist, 0, ierflg);
  gMinuit->mnexcm("SHO GRA",arglist,0,ierflg);

  //Doing timing tests:
  Double_t fvalue = 0;
  Double_t * thisparam = new Double_t[nparam_min];
  Double_t * grad = new Double_t[nparam_min];
  Double_t * hess = new Double_t[nparam_min*(nparam_min+1)/2];
  Double_t cputime;
  TStopwatch timer;

  //  cout << "ok evaluating likelihood now!\n";
  Int_t ifreepar = 0;
  for (UInt_t ipar = 0 ; ipar < nparam ;ipar++) {
    if (parfixed[ipar]==0) {
      thisparam[ifreepar] = param[ipar];
      ifreepar++;
    }
  }
  timer.Reset();
  timer.Start();
  LkhdFcn(nparam_min,grad,fvalue,thisparam,1,hess );
  timer.Stop();
   cputime = timer.CpuTime();
  printf("Likelihood at initial par=%8.4f.\n",fvalue);
  printf("Calculating the Likelihood took %8.4f seconds.\n",cputime);

  timer.Reset();
  timer.Start();
  LkhdFcn(nparam_min,grad,fvalue,thisparam,2,hess);
  timer.Stop();
   cputime = timer.CpuTime();
  printf("Calculating the Gradient   took %8.4f seconds.\n",cputime);

  if ((CalcHess ==1) || (initializing==1)) {
    timer.Reset();
    timer.Start();
    LkhdFcn(nparam_min,grad,fvalue,thisparam,3,hess);
    timer.Stop();
    cputime = timer.CpuTime();
    printf("Calculating the Hessian    took %8.4f seconds.\n",cputime);
  }
  
  delete [] thisparam;
  delete [] grad;
  delete [] hess;

  //Minimize
  arglist[0] = 9000000;  //max number of iterations
  arglist[1] = .00001;   // tolerance *10^-5
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  if (printlevel > -1)
    std::cout << "Status from MIGRAD: " << ierflg << std::endl;

  if (ierflg!=0) {
    cout << "OK, Minimizing without gradients" << endl;

    gMinuit->mnexcm("SET NOGRA", arglist, 0, ierflg);
    gMinuit->mnexcm("SHO GRA",arglist,0,ierflg);

    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    if (printlevel > -1)
      std::cout << "Final Status from MIGRAD: " << ierflg << std::endl;
  }

  Double_t thisPar, thisParErr;
  ipar = 0;
  for (UInt_t i = 0 ; i < nparam ; i++) {
    if (parfixed[i]==0) {    
      gMinuit->GetParameter(ipar, thisPar, thisParErr);
      if (i<nfac) {
	Setparam(i, fabs(thisPar));
	Setparam_err(i, thisParErr);
      }
      else if ((fac_nmix>1)&&(i>=f_nvariance+nfac+1)&&(i<f_nvariance+nfac+1+nfac)) {
	Setparam(i, fabs(thisPar));
	Setparam_err(i, thisParErr);
      }
      else if ((fac_nmix>2)&&(i>=2*(f_nvariance+nfac+1))&&(i<2*(f_nvariance+nfac+1)+nfac)) {
	Setparam(i, fabs(thisPar));
	Setparam_err(i, thisParErr);
      }
      // Force sign on means:
      else if ((fac_nmix>1)&&(i>=f_nvariance)&&(i<f_nvariance+nfac)) {
	Setparam(i, fabs(thisPar));
	Setparam_err(i, thisParErr);
      }
      else if ((fac_nmix>2)&&(i>=2*f_nvariance+nfac+1)&&(i<2*f_nvariance+nfac+1+nfac)) {
	//	param[i] = -1.0*fabs(thisPar);
	if (i==2*f_nvariance+nfac+1) Setparam(i, -1.0*fabs(thisPar));
	Setparam_err(i, thisParErr);
      }
      else {
	Setparam(i, thisPar);
	Setparam_err(i, thisParErr);
      }


//       Setparam(i, thisPar);
//       Setparam_err(i, thisParErr);
      ipar++;
    }
  }
  return ierflg;

  
  //  if (ierflg==0) {
    
    //Calculate the Hessian
//     gMinuit->mnexcm("HESSE", arglist, 0, ierflg);
//     if (printlevel > -1)
//       std::cout << "Status from HESSE: " << ierflg << std::endl;
    
    // Get MINOS errors 
    //   arglist[0] = 20000;
    //   for (int i = 1 ; i < 16 ; i++) arglist[i] = i;
    //  gMinuit->mnexcm("MINOS", arglist, 16, ierflg);
    //  gMinuit->mnexcm("MINOS", arglist, 0, ierflg);
    //     std::cout << "Status from MINOS: " << ierflg << std::endl;
    
    // Print results
//     Double_t amin,edm,errdef;
//     Int_t nvpar,nparx,icstat;
//     gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
//     if (printlevel > 0) gMinuit->mnprin(3,amin);
    
//     // Print results
//     Double_t thisPar, thisParErr;
//     //  Double_t thisParPlusErr, thisParMinusErr, thisParParabErr, thisParGcc;
//     for (UInt_t i = 0 ; i < nparam ; i++) {
//       gMinuit->GetParameter(i, thisPar, thisParErr);
//       //    gMinuit->mnerrs(i,thisParPlusErr,thisParMinusErr,
//       //		    thisParParabErr,thisParGcc);
//       //    if (printlevel > 0) std::cout << "param " << i << ": " << thisPar
//       //				  << " +/- " << thisParErr
//       //      //				  << " - " << thisParMinusErr 
//       //      //				  << " + " << thisParPlusErr
//       //				  << std::endl;
      
//       Setparam(i] = thisPar;
//       Setparam_err(i] = thisParErr;
//     }
//     cout << endl << "********Minimized parameter values:" << endl;
//     PrintParam();
//   }
//   else cout << "*******MINIMIZATION FAILED!!********" << endl;
//   cout << "Minimization took " << timer.CpuTime() << " s." << endl;

}


// Int_t TMinLkhd::Min_Knitro(Int_t printlevel) {

//   Int_t nparam_min = 0;
//   UInt_t ifreepar = 0;
//   vector<UInt_t> freeconstraints;
//   vector<Int_t> revfreeparlist;

//   if (predicting==0) {
//     //Calculate number of fixed parameters and make a lookup table of freeparameters
//     revfreeparlist.resize(nparam,-1);
//     Int_t nfixed = 0;
    
//     for (UInt_t i = 0 ; i < nparam; i++) {
//       if (parfixed.at(i)==0) {
// 	revfreeparlist.at(i) = ifreepar;
// 	ifreepar++;
//       }
//       nfixed += parfixed.at(i);
//     }
    
    
//     for (UInt_t i = 0 ; i < nparam; i++) {
      
//       // Is this parameter constrained?
//       if (parconstrained.at(i)!= -1) {
	
// 	// Check to see if both of the constrained parameters are free in this estimation
// 	if ((parfixed.at(i)==0)&&(parfixed.at(parconstrained.at(i))==0)) {
// 	  // We have a constraint!
// 	  freeconstraints.push_back(i);
// 	}
//       }
//     }


//     if (printlevel>0) printf("Minimizing %d parameters with %d constraints and while fixing %d.\n",nparam-nfixed,Int_t(freeconstraints.size()),nfixed);
//     //  cout << "Super-special counter of free parameters says there are: " << nfreeparam << ".\n";

//     nparam_min = nparam - nfixed;
//   }
//   else {
//     nparam_min = nfac;
//     revfreeparlist.resize(nfac,-1);
//     for (UInt_t i = 0 ; i < nfac; i++) revfreeparlist.at(i) = i;
//   }


//     /*---- CREATE A NEW KNITRO SOLVER INSTANCE. */

//   int status;
  
//   KTR_context_ptr  kc;
//   kc = KTR_new();
//   if (kc == NULL)
//     {
// 	cout << "*** KTR_new failed, maybe a license issue?\n";
// 	exit( EXIT_FAILURE );
//     }
  
//   status = KTR_load_param_file (kc, "knitro.opt");

//   if (printlevel==0) KTR_set_int_param (kc, KTR_PARAM_OUTLEV, KTR_OUTLEV_NONE);

//   if (status != 0)
//     {
//       cout << "*** KTR_load_param_file() returned " << status << "\n";
//       return( false );
//     }
//   else  if (printlevel>0) cout << "**knitro options set\n";

//     int       _nN;
//     double *  _daXInit = NULL;
//     double *  _daXLo = NULL;
//     double *  _daXUp = NULL;
//     int       _nM;
//     int    *  _naCType = NULL;
//     double *  _daCLo = NULL;
//     double *  _daCUp = NULL;
//     int       _nNnzJ;
//     int    *  _naJacIndexVars = NULL;
//     int    *  _naJacIndexCons = NULL;
//     int       _nNnzH;
//     int    *  _naHessCols = NULL;
//     int    *  _naHessRows = NULL;


//     _nN = nparam_min; // number of parameters
//     _nM = freeconstraints.size(); // number of constraints
//     _nNnzJ = 2*_nM;

//     //---- VARIABLE BOUNDS
//     _daXLo  = new double[_nN];
//     _daXUp  = new double[_nN];

//     //---- PUT THE CONSTANT TERM IN THE RIGHT-HAND SIDE.
//     if (_nM>0) {
//       _naCType  = new int[_nM];
//       _daCLo    = new double[_nM];
//       _daCUp    = new double[_nM];
//       //---- SPECIFY THE CONSTRAINT JACOBIAN SPARSITY STRUCTURE.
//       _naJacIndexVars = new int[_nNnzJ];
//       _naJacIndexCons = new int[_nNnzJ];
      
//       for (int i = 0 ; i < _nM; i++) {
// 	_daCLo[i] = 0.0;
// 	_daCUp[i] = 0.0;
// 	_naCType[i] = KTR_CONTYPE_LINEAR;
	
// 	_naJacIndexCons[2*i] = i;  
// 	_naJacIndexVars[2*i] = revfreeparlist.at(freeconstraints.at(i));
	
// 	_naJacIndexCons[2*i+1] = i;  
// 	_naJacIndexVars[2*i+1] = revfreeparlist.at(parconstrained.at(freeconstraints.at(i)));
//       }
//     }

//     //---- SPECIFY THE HESSIAN OF THE LAGRANGIAN SPARSITY STRUCTURE.
//     //    _nNnzH = 0;
//     _nNnzH = nparam_min*(nparam_min+1)/2;
//     _naHessRows = new int[_nNnzH];
//     _naHessCols = new int[_nNnzH];
//     Int_t counter = 0;
//     for (Int_t irow = 0 ; irow < nparam_min ; irow++) {
//       for (Int_t icol = irow ; icol < nparam_min ; icol++) {
// 	_naHessRows[counter] = irow;  _naHessCols[counter] = icol;
// 	counter++;
//       }
//     }

//     //---- INITIAL GUESS FOR x AND lambda.
//     _daXInit = new double[_nN];
//     double * daLambdaInit = NULL;
// //     double *  daLambdaInit = new double[_nM + _nN];
// //     for (int i = 0; i < _nM + _nN; i++)
// //         daLambdaInit[i] = 0.0;


//   // Tell knitro about the parameter values and bounds:
//     if (predicting==0) {
//       Int_t ipar = 0;
//       for (UInt_t i = 0 ; i < nparam ; i++) {
// 	if (parfixed[i]==0) {
	  
// 	  if ((fac_corr!=0) && (i>=nfac) && (i<f_nvariance)) {
// 	    _daXLo[ipar] = -1.0;
// 	    _daXUp[ipar] = 1.0;
// 	    _daXInit[ipar] = param[i];
// 	  }
// 	  // set limits for covariances for second mixture
// 	  else if ((fac_nmix>1) && (fac_corr!=0) && (i>=f_nvariance+nfac+1+nfac) && (i<f_nvariance+nfac+1+f_nvariance)) {
// 	    _daXLo[ipar] = -1.0;
// 	    _daXUp[ipar] = 1.0;
// 	    _daXInit[ipar] = param[i];
// 	  }
// 	  // set limits for covariances for third mixture
// 	  else if ((fac_nmix>2) && (fac_corr!=0) && (i>=2*(f_nvariance+nfac+1)+nfac) && (i<2*(f_nvariance+nfac+1)+f_nvariance)) {
// 	    _daXLo[ipar] = -1.0;
// 	    _daXUp[ipar] = 1.0;
// 	    _daXInit[ipar] = param[i];
// 	  }
// 	  else {
// 	    _daXLo[ipar] = -KTR_INFBOUND;
// 	    _daXUp[ipar] = KTR_INFBOUND;
// 	    _daXInit[ipar] = param[i];
// 	  }
// 	  ipar++;
// 	}
//       }
//     }
//     else {
//       // set initial value and bounds for factor prediction
//       for (UInt_t ipar = 0; ipar<nfac; ipar++) {
// 	_daXLo[ipar] = -KTR_INFBOUND;
// 	_daXUp[ipar] = KTR_INFBOUND;
// 	_daXInit[ipar] = 0.01;
//       }
//     }


//   double * jac = NULL;
//   double * c = NULL;

//   if (_nM>0) {
//     c = new double[_nM];
//     jac = new double[_nNnzJ];

//     for (int i = 0 ; i < _nM; i++) {
//       jac[2*i] = 1.0;
//       jac[2*i+1] = -1.0;
//       printf("(%d, %d): %f\n",_naJacIndexCons[2*i], _naJacIndexVars[2*i], jac[2*i]);
//       printf("(%d, %d): %f\n",_naJacIndexCons[2*i+1], _naJacIndexVars[2*i+1], jac[2*i+1]);
//     }
//   }
//   double * daLambda = new double[_nN+_nM];

//   Double_t fvalue = 0;
//   Double_t * thisparam = new Double_t[_nN];
//   Double_t * grad = new Double_t[_nN];
//   Double_t * hess = new Double_t[_nNnzH];

//   // Check gradient:
// //     status =1 ;
// //     int counter = 0;
// //     while ((status==1)||(status==2)) {
// //       LkhdFcn(num_param,grad,fvalue,_daXInit,status);
// //       //  EvalLkhd(fvalue, grad, parinit, 2,0,1);

// //       status = KTR_check_first_ders (kc, _daXInit, 2, 2.0e-6, 2.0e-6,
// // 				     0, fvalue, c, grad, jac, NULL);
// //       cout << "Pass " << counter << ". status=" << status << "\n";
// //       counter++;
// //     }

//   if (predicting==0) {
//     ifreepar = 0;
//     for (UInt_t ipar = 0 ; ipar < nparam ;ipar++) {
//       if (parfixed[ipar]==0) {
// 	thisparam[ifreepar] = param[ipar];
// 	ifreepar++;
//       }
//     }
    
//     for (UInt_t i = 0 ; i < nparam ; i++) {
//       if (parfixed[i]==0) {    
// 	Setparam_err(i] = 0.0;
//       }
//     }
//   }
//   else {
//     for (UInt_t ifac = 0; ifac < nfac ; ifac++) thisparam[ifac] = 0.01; 
//   }

//   Double_t cputime;
//   TStopwatch timer;

//   if (printlevel>0)  {
//     timer.Reset();
//     timer.Start();
//     LkhdFcn(nparam_min,grad,fvalue,thisparam,1,hess );
//     timer.Stop();
//     cputime = timer.CpuTime();
//     printf("Likelihood at initial par=%8.4f.\n",fvalue);
//     printf("Calculating the Likelihood took %8.4f seconds.\n",cputime);
    
//     timer.Reset();
//     timer.Start();
//     LkhdFcn(nparam_min,grad,fvalue,thisparam,2,hess);
//     timer.Stop();
//     cputime = timer.CpuTime();
//     printf("Calculating the Gradient   took %8.4f seconds.\n",cputime);
    
//     timer.Reset();
//     timer.Start();
//     LkhdFcn(nparam_min,grad,fvalue,thisparam,3,hess);
//     timer.Stop();
//     cputime = timer.CpuTime();
//     printf("Calculating the Hessian    took %8.4f seconds.\n",cputime);
//   }

//   if (printlevel>0) cout << "******************ok now minimizing with knitro\n";
//   //  cout << "ok evaluating likelihood now!\n";
//   //  for (UInt_t ipar = 0 ; ipar < nparam ;ipar++) thisparam[ipar] = param[ipar];

//   status = KTR_init_problem (kc, _nN,
// 			     KTR_OBJGOAL_MINIMIZE, KTR_OBJTYPE_GENERAL,
// 			     _daXLo, _daXUp,
// 			     _nM, _naCType, _daCLo, _daCUp,
// 			     _nNnzJ, _naJacIndexVars, _naJacIndexCons,
// 			     _nNnzH, _naHessRows, _naHessCols,
// 			     _daXInit, daLambdaInit);
//   if (status != 0)
//     {
//       cout << "*** KTR_init_problem() returned " << status << "\n";
//       return( status);
//     }
//   else if (printlevel>0) cout << "**knitro initialized\n";


//   delete [] _daXLo;
//   delete [] _daXUp;
//   if (_nM>0) {
//     delete [] _naCType;
//     delete [] _daCLo;
//     delete [] _daCUp;
//     delete [] _naJacIndexVars;
//     delete [] _naJacIndexCons;
//   }
//   delete [] _naHessRows;
//   delete [] _naHessCols;
//   //  delete [] daLambdaInit;
//   delete [] _daXInit;


//   status =1 ;
//   counter = 0;
//   Int_t badstatuscounter = 0;
//   while ((status==1)||(status==2)||(status==3)) {
//     if (counter>0) { 

//       LkhdFcn(_nN,grad,fvalue,thisparam,status,hess);

//       for (int i = 0 ; i < _nM ; i++) {
// 	c[i] = thisparam[revfreeparlist.at(freeconstraints.at(i))] 
// 	  - thisparam[revfreeparlist.at(parconstrained.at(freeconstraints.at(i)))];
//       }
//     }
//     status = KTR_solve (kc, thisparam , daLambda, 0, &fvalue,
//   			c, grad, jac, hess, NULL, NULL);    


// //     if (counter>0) LkhdFcn(_nN,grad,fvalue,thisparam,status,NULL);
// //     status = KTR_solve (kc, _daXInit , daLambda, 0, &fvalue,
// //  			c, grad, jac, NULL, NULL, NULL);    
//     counter++;
//     if ((status==-100)||(status==-101)||(status==-102)) {
//       printf("KNITRO FAILED WITH STATUS %d, restart %d...\n",status,badstatuscounter);
//       Double_t nudge = 0.001;
//       if (badstatuscounter==1) nudge = -0.001;
//       else if (badstatuscounter==2) nudge = 0.01;
//       else if (badstatuscounter==3) nudge = -0.01;
//       else if (badstatuscounter==4) nudge =  0.1;
      
//       if (badstatuscounter<5) {
// 	for (Int_t ipar = 0 ; ipar < _nN ;ipar++) thisparam[ipar] += nudge;
// 	status= KTR_restart(kc,thisparam,daLambda);
// 	status = 1;
// 	badstatuscounter++;
//       }
//     }
//   }
  
// // int  KNITRO_API KTR_solve ( KTR_context_ptr      kc,
// //                                  double * const  x,
// //                                  double * const  lambda,
// //                            const int             evalStatus,
// //                                  double * const  obj,
// //                            const double * const  c,
// //                                  double * const  objGrad,
// //                                  double * const  ac,
// //                            const double * const  hess,
// //                                  double * const  hessVector,
// //                                  void   * const  userParams);

//   if (status==0) {
//     if (printlevel>0) cout << "Minimization done. Lkhd= " << fvalue << "\n";

//     if (printlevel>0) cout << "Calculating Standard Errors...\n";

//     for (UInt_t i = 0 ; i < nparam ; i++) {
//       if (parfixed[i]==0) {
// 	Setparam_err(i] = 0.0;
//       }
//     }

//     Int_t diagneg = 0;
//     if (HessStdErr>0) {
//       LkhdFcn(_nN,grad,fvalue,thisparam,3,hess);
      
//       TMatrixD invhess(_nN,_nN);
//       counter = 0;
//       for (Int_t irow = 0 ; irow < _nN ; irow++) {
// 	for (Int_t icol = irow ; icol < _nN ; icol++) {
// 	  invhess(irow,icol) = hess[counter];
// 	  if (irow!=icol) invhess(icol,irow) = hess[counter];
// 	  counter++;
// 	}
//       }
      
//       if (printlevel>0) cout << "Inverting Hessian...\n";
//       invhess.Invert();
//       //      Int_t diagneg = 0;
//       for (Int_t irow = 0 ; irow < _nN ; irow++) {
// 	if (invhess(irow,irow)<0) { 
// 	  diagneg++;
// 	  cout << "Found negative diagonal element in row/parameter number  " << revfreeparlist.at(irow) << "\n";
// 	}
//       }
      
//       if (diagneg>0) cout << "OOPS, " << diagneg << "/" << _nN << " diagonal elements are negative!!\n";
      
//       if (diagneg==0) {

// 	if (predicting ==0) {
// 	  Int_t ipar = 0;
// 	  for (UInt_t i = 0 ; i < nparam ; i++) {
// 	    if (parfixed[i]==0) {    
// 	      if (i<nfac) {
// 		param[i] = fabs(thisparam[ipar]);
// 		if (invhess(ipar,ipar)>=0) Setparam_err(i] = sqrt(invhess(ipar,ipar));
// 	      }
// 	      else if ((fac_nmix>1)&&(i>=f_nvariance+nfac+1)&&(i<f_nvariance+nfac+1+nfac)) {
// 		param[i] = fabs(thisparam[ipar]);
// 		if (invhess(ipar,ipar)>=0) Setparam_err(i] = sqrt(invhess(ipar,ipar));
// 	      }
// 	      else if ((fac_nmix>2)&&(i>=2*(f_nvariance+nfac+1))&&(i<2*(f_nvariance+nfac+1)+nfac)) {
// 		param[i] = fabs(thisparam[ipar]);
// 		if (invhess(ipar,ipar)>=0) Setparam_err(i] = sqrt(invhess(ipar,ipar));
// 	      }
// 	      // Force sign on means:
// 	      else if ((fac_nmix>1)&&(i>=f_nvariance)&&(i<f_nvariance+nfac)) {
// 		param[i] = fabs(thisparam[ipar]);
// 		if (invhess(ipar,ipar)>=0) Setparam_err(i] = sqrt(invhess(ipar,ipar));
// 	      }
// 	      else if ((fac_nmix>2)&&(i>=2*f_nvariance+nfac+1)&&(i<2*f_nvariance+nfac+1+nfac)) {
// 		param[i] = -1.0*fabs(thisparam[ipar]);
// 		//if (i==2*f_nvariance+nfac+1) param[i] = -1.0*fabs(thisparam[ipar]);
// 		if (invhess(ipar,ipar)>=0) Setparam_err(i] = sqrt(invhess(ipar,ipar));
// 	      }
// 	      else {
// 		param[i] = thisparam[ipar];
// 		if (invhess(ipar,ipar)>=0) Setparam_err(i] = sqrt(invhess(ipar,ipar));
// 	      }
// 	      ipar++;
// 	    }
// 	  }
// 	} // Predicting==0
// 	else {

// 	  // append factor predictions for this observation
// 	  FILE * pFile;
// 	  TString filename;
// 	  filename = TString("factor_predictions.txt");
// 	  filename.Prepend(workingdir);
// 	  pFile = fopen(filename.Data(),"a");

// 	  if (simIncData==0) {
// 	    if (indexvar==-1) {
// 	      fprintf(pFile,"%10d ", predictobs);
// 	    }
// 	    else {
// 	      fprintf(pFile,"%10d ", Int_t(data[predictobs*nvar+indexvar]));
// 	    }
// 	    // append the factor predictions
// 	    for (UInt_t ifac = 0; ifac < nfac ; ifac++) fprintf(pFile,"%10.5f ",thisparam[ifac]);
// 	    //	    for (UInt_t ifac = 0; ifac < nfac ; ifac++) fprintf(pFile,"%10.5f %10.5f ",thisparam[ifac],sqrt(invhess(ifac,ifac)));
// 	  }
// 	  else {
// 	    for (UInt_t ivar = 0 ; ivar < nvar ; ivar++) {
// 	      fprintf(pFile,"%10.5f ",data[predictobs*nvar+ivar]);
// 	    }
// 	    for (UInt_t ifac = 0; ifac < nfac ; ifac++) fprintf(pFile,"%10.5f ",thisparam[ifac]);
// 	  }

// 	  fprintf(pFile,"\n");
// 	  fclose (pFile);
// 	}
//       }
//     }

//     if ( (((diagneg>0)&&(HessStdErr==2)) || (HessStdErr==0)) && (predicting==0) ) {

//       //      if ((nsubsamples==0)&&(nbootsamples==0)) assert(0);

//       Int_t ipar = 0;
//       for (UInt_t i = 0 ; i < nparam ; i++) {
// 	if (parfixed[i]==0) {    
// 	  if (i<nfac) {
// 	    param[i] = fabs(thisparam[ipar]);
// 	  }
// 	  else if ((fac_nmix>1)&&(i>=f_nvariance+nfac+1)&&(i<f_nvariance+nfac+1+nfac)) {
// 	    param[i] = fabs(thisparam[ipar]);
// 	  }
// 	  else if ((fac_nmix>2)&&(i>=2*(f_nvariance+nfac+1))&&(i<2*(f_nvariance+nfac+1)+nfac)) {
// 	    param[i] = fabs(thisparam[ipar]);
// 	  }
// 	  // Force sign on means:
// 	  else if ((fac_nmix>1)&&(i>=f_nvariance)&&(i<f_nvariance+nfac)) {
// 	    param[i] = fabs(thisparam[ipar]);
// 	  }
// 	  else if ((fac_nmix>2)&&(i>=2*f_nvariance+nfac+1)&&(i<2*f_nvariance+nfac+1+nfac)) {
// 	    param[i] = -1.0*fabs(thisparam[ipar]);
// 	    //if (i==2*f_nvariance+nfac+1) param[i] = -1.0*fabs(thisparam[ipar]);
// 	  }
// 	  else {
// 	    param[i] = thisparam[ipar];
// 	  }
// 	  Setparam_err(i] = CalcStdErrLkhdRatio(i);

// 	  ipar++;
// 	}
//       }
//     }

//   }
//   else cout << "*******MINIMIZATION FAILED!!********" << endl;

//   KTR_free (&kc);
  
//   delete [] jac;
//   delete [] daLambda;
//   delete [] c;
//   delete [] thisparam;
//   delete [] grad;
//   delete [] hess;
  
//   return status;
// }


Int_t TMinLkhd::Min_Ipopt(Int_t printlevel) {

  Int_t nparam_min = 0;
  UInt_t ifreepar = 0;
  vector<Int_t> revfreeparlist;

  if (predicting==0) {
    //Calculate number of fixed parameters and make a lookup table of freeparameters
    revfreeparlist.resize(nparam,-1);
    Int_t nfixed = 0;
    
    for (UInt_t i = 0 ; i < nparam; i++) {
      if (parfixed.at(i)==0) {
	revfreeparlist.at(i) = ifreepar;
	ifreepar++;
      }
      nfixed += (parfixed.at(i)>0);
    }


    printf("Minimizing %d parameters while fixing %d (nfreeparam=%d).\n",nparam-nfixed,nfixed,nfreeparam);
    //        if (printlevel>0) printf("Minimizing %d parameters while fixing %d (nfreeparam=%d).\n",nparam-nfixed,nfixed,nfreeparam);

    nparam_min = nparam - nfixed;
  }
  else {
    nparam_min = nfac;
    revfreeparlist.resize(nfac,-1);
    for (UInt_t i = 0 ; i < nfac; i++) revfreeparlist.at(i) = i;
  }


    int       _nN;
    double *  _daXInit = NULL;
    double *  _daXLo = NULL;
    double *  _daXUp = NULL;
    int       _nNnzH;
    int    *  _naHessCols = NULL;
    int    *  _naHessRows = NULL;

    _nN = nparam_min; // number of parameters

    //---- VARIABLE BOUNDS
    _daXLo  = new double[_nN];
    _daXUp  = new double[_nN];


    //---- SPECIFY THE HESSIAN OF THE LAGRANGIAN SPARSITY STRUCTURE.
    //    _nNnzH = 0;
    _nNnzH = nparam_min*(nparam_min+1)/2;
    _naHessRows = new int[_nNnzH];
    _naHessCols = new int[_nNnzH];
    Int_t counter = 0;
    for (Int_t irow = 0 ; irow < nparam_min ; irow++) {
      for (Int_t icol = irow ; icol < nparam_min ; icol++) {
	_naHessRows[counter] = irow;  _naHessCols[counter] = icol;
	counter++;
      }
    }

    //---- INITIAL GUESS FOR x AND lambda.
    _daXInit = new double[_nN];

    Double_t fvalue = 0;
    Double_t * thisparam = new Double_t[_nN];
    Double_t * grad = new Double_t[_nN];
    Double_t * hess = new Double_t[_nNnzH];


  // Tell Ipopt about the parameter values and bounds:
    if (predicting==0) {
      Int_t ipar = 0;
      for (UInt_t i = 0 ; i < nparam ; i++) {
	if (parfixed[i]==0) {
	  
	  if ((fac_corr!=0) && (i>=nfac) && (i<f_nvariance)) {
	    _daXLo[ipar] = -1.0;
	    _daXUp[ipar] = 1.0;
	    _daXInit[ipar] = param[i];
	  }
	  // set limits for covariances for second mixture
	  else if ((fac_nmix>1) && (fac_corr!=0) && (i>=f_nvariance+nfac+1+nfac) && (i<f_nvariance+nfac+1+f_nvariance)) {
	    _daXLo[ipar] = -1.0;
	    _daXUp[ipar] = 1.0;
	    _daXInit[ipar] = param[i];
	  }
	  // set limits for covariances for third mixture
	  else if ((fac_nmix>2) && (fac_corr!=0) && (i>=2*(f_nvariance+nfac+1)+nfac) && (i<2*(f_nvariance+nfac+1)+f_nvariance)) {
	    _daXLo[ipar] = -1.0;
	    _daXUp[ipar] = 1.0;
	    _daXInit[ipar] = param[i];
	  }
	  else {
	    _daXLo[ipar] = -2e19;
	    _daXUp[ipar] = 2e19;
	    _daXInit[ipar] = param[i];
	  }
	  ipar++;
	}
      }
    }
    else {
      // set initial value and bounds for factor prediction
      for (UInt_t ipar = 0; ipar<nfac; ipar++) {
	_daXLo[ipar] = -2e19;
	_daXUp[ipar] = 2e19;
	_daXInit[ipar] = 0.01;
      }
    }

  if (predicting==0) {
    ifreepar = 0;
    for (UInt_t ipar = 0 ; ipar < nparam ;ipar++) {
      if (parfixed[ipar]==0) {
	thisparam[ifreepar] = param[ipar];
	ifreepar++;
      }
    }
    
    for (UInt_t i = 0 ; i < nparam ; i++) {
      if (parfixed[i]==0) {    
	Setparam_err(i,  0.0);
      }
    }
  }
  else {
    for (UInt_t ifac = 0; ifac < nfac ; ifac++) thisparam[ifac] = 0.01; 
  }


  Double_t cputime;
  TStopwatch timer;

  if (printlevel>0)  {
    timer.Reset();
    timer.Start();
    printf("Evaluating Likelihood\n");
    LkhdFcn(nparam_min,grad,fvalue,thisparam,1,hess );
    timer.Stop();
    cputime = timer.CpuTime();
    printf("Likelihood at initial par=%8.4f.\n",fvalue);
    printf("Calculating the Likelihood took %8.4f seconds.\n",cputime);
    
    timer.Reset();
    timer.Start();
    LkhdFcn(nparam_min,grad,fvalue,thisparam,2,hess);
    timer.Stop();
    cputime = timer.CpuTime();
    printf("Calculating the Gradient   took %8.4f seconds.\n",cputime);
    // printf("The gradient is:\n");
    // for (int ipar = 0 ; ipar < nparam_min ; ipar++) printf("Par Grad #%4d: %8.5f\n",ipar,grad[ipar]);


    if ((CalcHess ==1) || (initializing==1)) {
      timer.Reset();
      timer.Start();
      LkhdFcn(nparam_min,grad,fvalue,thisparam,3,hess);
      timer.Stop();
      cputime = timer.CpuTime();
      printf("Calculating the Hessian    took %8.4f seconds.\n",cputime);
    }
      // printf("The Hessian is:\n");
    // for (int ipar = 0 ; ipar < nparam_min*(nparam_min+1)/2 ; ipar++) printf("Hess element #%4d: %8.5f\n",ipar,hess[ipar]);

  }

  if (printlevel>0) cout << "******************ok now minimizing with ipopt\n";
  //  cout << "ok evaluating likelihood now!\n";
  //  for (UInt_t ipar = 0 ; ipar < nparam ;ipar++) thisparam[ipar] = param[ipar];


    /*---- CREATE A NEW IPOPT SOLVER INSTANCE. */

  SmartPtr<MINLKHD_NLP> mynlp = new MINLKHD_NLP(_nN,_daXLo, _daXUp, _daXInit,
						_nNnzH, _naHessRows, _naHessCols);

  //  SmartPtr<TNLP> mynlp = new HS071_NLP();

  SmartPtr<IpoptApplication> app = new IpoptApplication();

  //  app->Options()->SetNumericValue("tol", 1e-7);
  if (predicting==1){
    app->Options()->SetIntegerValue("print_level", printlevel+1);
  }
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("nlp_scaling_method","none");
  app->Options()->SetNumericValue("obj_scaling_factor",1.0);
  //Check gradient:
  //  app->Options()->SetStringValue("derivative_test","second-order");
  app->Options()->SetStringValue("linear_solver", "ma57");
  //  app->Options()->SetStringValue("linear_solver", "mumps");

  if ( (CalcHess == 0) && (initializing==0) )  app->Options()->SetStringValue("hessian_approximation", "limited-memory"); //default is exact
  //   app->Options()->SetStringValue("hessian_approximation", "exact"); //default is exact

//  UInt_t estimating_outcomes = 0;
//  for (UInt_t imod = 0 ; imod < models.size() ; imod++) {
//    // Use the fact that all outcome models are ignored when estimating measurement system
//    if ((models[imod].GetGroup()>0)&&(models[imod].GetIgnore()==0)) estimating_outcomes = 1;
//  }

  if (max_cpu_limit>0.0) app->Options()->SetNumericValue("max_cpu_time", max_cpu_limit);

  app->Options()->SetStringValue("output_file", "ipopt.out");
  app->Options()->SetIntegerValue("max_iter", 1000000);

  // The following overwrites the default name (ipopt.opt) of the
  // options file
  // app->Options()->SetStringValue("option_file_name", "costmin.opt");
  
  // Initialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
    return (int) status;
  }

  delete [] _daXLo;
  delete [] _daXUp;
  delete [] _naHessRows;
  delete [] _naHessCols;
  delete [] _daXInit;



  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(mynlp);
  Int_t badstatuscounter = 0;
  while ((status != Solve_Succeeded) && (status != Solved_To_Acceptable_Level) && (badstatuscounter<6)) {

    printf("IPOPT FAILED WITH STATUS %d, restart %d...\n",status,badstatuscounter);
    Double_t nudge = 0.9;
    if (badstatuscounter==1) nudge = 1.1;
    else if (badstatuscounter==2) nudge = 0.8;
    else if (badstatuscounter==3) nudge = 1.2;
    else if (badstatuscounter==4) nudge =  0.7;
    else if (badstatuscounter==5) nudge =  1.3;

    mynlp->restart(nudge);
    status = app->OptimizeTNLP(mynlp);

    badstatuscounter++;
  }

  if ((status != Solve_Succeeded) && (status != Solved_To_Acceptable_Level)) {
    std::cout << std::endl << std::endl << "*** Minimization failed with code " << status << std::endl;
    return (int) status;
  }
  

//   status =1 ;
//   counter = 0;
//   Int_t badstatuscounter = 0;
//   while ((status==1)||(status==2)||(status==3)) {
//     if (counter>0) { 

//       LkhdFcn(_nN,grad,fvalue,thisparam,status,hess);

//       for (int i = 0 ; i < _nM ; i++) {
// 	c[i] = thisparam[revfreeparlist.at(freeconstraints.at(i))] 
// 	  - thisparam[revfreeparlist.at(parconstrained.at(freeconstraints.at(i)))];
//       }
//     }
//     status = KTR_solve (kc, thisparam , daLambda, 0, &fvalue,
//   			c, grad, jac, hess, NULL, NULL);    


// //     if (counter>0) LkhdFcn(_nN,grad,fvalue,thisparam,status,NULL);
// //     status = KTR_solve (kc, _daXInit , daLambda, 0, &fvalue,
// //  			c, grad, jac, NULL, NULL, NULL);    
//     counter++;
//     if ((status==-100)||(status==-101)||(status==-102)) {
//       printf("KNITRO FAILED WITH STATUS %d, restart %d...\n",status,badstatuscounter);
//       Double_t nudge = 0.001;
//       if (badstatuscounter==1) nudge = -0.001;
//       else if (badstatuscounter==2) nudge = 0.01;
//       else if (badstatuscounter==3) nudge = -0.01;
//       else if (badstatuscounter==4) nudge =  0.1;
      
//       if (badstatuscounter<5) {
// 	for (Int_t ipar = 0 ; ipar < _nN ;ipar++) thisparam[ipar] += nudge;
// 	status= KTR_restart(kc,thisparam,daLambda);
// 	status = 1;
// 	badstatuscounter++;
//       }
//     }
//   }
  
// int  KNITRO_API KTR_solve ( KTR_context_ptr      kc,
//                                  double * const  x,
//                                  double * const  lambda,
//                            const int             evalStatus,
//                                  double * const  obj,
//                            const double * const  c,
//                                  double * const  objGrad,
//                                  double * const  jac,
//                            const double * const  hess,
//                                  double * const  hessVector,
//                                  void   * const  userParams);

  if ((status==1)||(status==0)) {

    mynlp->getparam(thisparam);

    if (printlevel>0) cout << "Minimization done. Calculating Standard Errors...\n";

    for (UInt_t i = 0 ; i < nparam ; i++) {
      if (parfixed[i]==0) {
	Setparam_err(i,  0.0);
      }
    }

    Int_t diagneg = 0;
    if (HessStdErr>0) {

      LkhdFcn(_nN,grad,fvalue,thisparam,3,hess);
      if (printlevel>0) cout << "Lkhd= " << fvalue << "\n";

    
      TMatrixD invhess(_nN,_nN);
      counter = 0;
      for (Int_t irow = 0 ; irow < _nN ; irow++) {
	for (Int_t icol = irow ; icol < _nN ; icol++) {
	  invhess(irow,icol) = hess[counter];
	  if (irow!=icol) invhess(icol,irow) = hess[counter];
	  counter++;
	}
      }
      
      if (printlevel>0) cout << "Inverting Hessian...\n";
      invhess.Invert();
      //      Int_t diagneg = 0;
      for (Int_t irow = 0 ; irow < _nN ; irow++) {
	if (invhess(irow,irow)<0) { 
	  diagneg++;
	  cout << "Found negative diagonal element in row/parameter number  " << revfreeparlist.at(irow) << "\n";
	}
      }
      
      if (diagneg>0) cout << "OOPS, " << diagneg << "/" << _nN << " diagonal elements are negative!!\n";
      
      if (diagneg==0) {

	if (predicting ==0) {
	  Int_t ipar = 0;
	  for (UInt_t i = 0 ; i < nparam ; i++) {
	    if (parfixed[i]==0) {    
	      if (i<nfac) {
		Setparam(i, fabs(thisparam[ipar]));
		if (invhess(ipar,ipar)>=0) Setparam_err(i, sqrt(invhess(ipar,ipar)));
	      }
	      else if ((fac_nmix>1)&&(i>=f_nvariance+nfac+1)&&(i<f_nvariance+nfac+1+nfac)) {
		Setparam(i, fabs(thisparam[ipar]));
		if (invhess(ipar,ipar)>=0) Setparam_err(i, sqrt(invhess(ipar,ipar)));
	      }
	      else if ((fac_nmix>2)&&(i>=2*(f_nvariance+nfac+1))&&(i<2*(f_nvariance+nfac+1)+nfac)) {
		Setparam(i, fabs(thisparam[ipar]));
		if (invhess(ipar,ipar)>=0) Setparam_err(i, sqrt(invhess(ipar,ipar)));
	      }
	      // Force sign on means:
	      else if ((fac_nmix>1)&&(i>=f_nvariance)&&(i<f_nvariance+nfac)) {
		Setparam(i, fabs(thisparam[ipar]));
		if (invhess(ipar,ipar)>=0) Setparam_err(i, sqrt(invhess(ipar,ipar)));
	      }
	      else if ((fac_nmix>2)&&(i>=2*f_nvariance+nfac+1)&&(i<2*f_nvariance+nfac+1+nfac)) {
		Setparam(i, -1.0*fabs(thisparam[ipar]));
		//if (i==2*f_nvariance+nfac+1) param[i, -1.0*fabs(thisparam[ipar]));
		if (invhess(ipar,ipar)>=0) Setparam_err(i, sqrt(invhess(ipar,ipar)));
	      }
	      else {
		Setparam(i, thisparam[ipar]);
		if (invhess(ipar,ipar)>=0) Setparam_err(i, sqrt(invhess(ipar,ipar)));
	      }
	      ipar++;
	    }
	  }
	} // Predicting==0
	else {

	  Double_t * trialparam = new Double_t[_nN];
	  Double_t * Lratio = new Double_t[_nN];
	  Double_t * sigmai = new Double_t[_nN];
	  double baseline = fvalue;

	  
	  for (Int_t ifac = 0 ; ifac < _nN ; ifac++) {
	    double minlratio = 10.0;
	    sigmai[ifac] = sqrt(invhess(ifac,ifac));
	    
	    for (Int_t jfac = 0 ; jfac < _nN ; jfac++) trialparam[jfac] = thisparam[jfac];

	    while (minlratio>1.0) {
	      
	      trialparam[ifac] = thisparam[ifac] + sigmai[ifac];
	      LkhdFcn(_nN,grad,fvalue,trialparam,1,hess);
	      if (fvalue - baseline < minlratio) minlratio = fvalue - baseline;

	      trialparam[ifac] = thisparam[ifac] - sigmai[ifac];
	      LkhdFcn(_nN,grad,fvalue,trialparam,1,hess);
	      if (fvalue - baseline < minlratio) minlratio = fvalue - baseline;

	      if (minlratio>1.0) sigmai[ifac] = 0.8*sigmai[ifac];
	    }

	    //	    trialparam[ifac] = 0.0;
	    //	    LkhdFcn(_nN,grad,fvalue,trialparam,1,hess);
	    Lratio[ifac] = minlratio;
	    
	  }
	  
	  
	  // Save factor predictions
	  facprediction.resize(nfac,0.0);
	  for (UInt_t ifac = 0; ifac < nfac ; ifac++) {

	    // Reset extreme values that have scores>sqrt(2)*4.5*sigma
	    // This is the most extreme point used with Gauss-Hermite quadrature (n=16)
	    // I am guessing this will not happen as often once we include the prior (factor distribution)
	    // when estimating the factor scores. It is pretty rare now.
	    if ((fac_nmix==1) && (fabs(thisparam[ifac]) > 4.688738939305818364688*sqrt(2.0)*Getfvar(0,ifac)) ) {
	      facprediction.at(ifac) = 4.688738939305818364688*sqrt(2.0)*Getfvar(0,ifac);
	      if (thisparam[ifac]<0) facprediction.at(ifac) = -1.0*facprediction.at(ifac);
	    }
	    else {
	      facprediction.at(ifac) = thisparam[ifac];
	    }
	  }
	  
	  
	  // print factor predictions for this observation
	    FILE * pFile;
	    TString filename;

	    if ((nsubsamples==0)&&(nbootsamples==0)) {
	      filename = TString("factor_predictions.txt");
	    }
	    else {
	      std::stringstream out;
	      out << GetCurrentSample();
	      filename = TString("factor_predictions_").Append(out.str()).Append(".txt");
	    }

	    //	    filename = TString("factor_predictions.txt");
	    filename.Prepend(workingdir);
	    pFile = fopen(filename.Data(),"a");
	    if (indexvar==-1) {
	      fprintf(pFile,"%10d ", predictobs);
	    }
	    else {
	      fprintf(pFile,"%10d ", Int_t(data[predictobs*nvar+indexvar]));
	    }
	    // append the factor predictions
	    for (UInt_t ifac = 0; ifac < nfac ; ifac++) {
	      //	      fprintf(pFile,"%10.5f %10.5f %10.5f ",thisparam[ifac],sigmai[ifac], Lratio[ifac]);
	      fprintf(pFile,"%10.5f %10.5f ",facprediction[ifac],sigmai[ifac]);
	    }

	    // if (predictobs<10) {
	    //   fprintf(pFile,"\n Likelihood for this observation:\n");
	    //   for (float igrid = 0.0 ; igrid<40.0 ; igrid++) { 
	    // 	fprintf(pFile,"%10.5f ", -2.0 + igrid/10.0);
		
	    // 	for (Int_t ifac = 0 ; ifac < _nN ; ifac++) {
	    // 	  for (Int_t jfac = 0 ; jfac < _nN ; jfac++) trialparam[jfac] = thisparam[jfac];
	    // 	  //		trialparam[ifac] += sqrt(invhess(ifac,ifac));
	    // 	  trialparam[ifac] = -2.0 + igrid/10.0;
	    // 	  LkhdFcn(_nN,grad,fvalue,trialparam,1,hess);
	    // 	  fprintf(pFile,"%10.5f ", fvalue);
	    // 	}
	    // 	fprintf(pFile,"\n");
	    //   }
	    //	  }
	    
	    fprintf(pFile,"\n");
	    fclose (pFile);

	  delete [] Lratio;
	  delete [] trialparam;
	}
      }
    }

    if ( (((diagneg>0)&&(HessStdErr==2)) || (HessStdErr==0)) && (predicting==0) ) {

      //      if ((nsubsamples==0)&&(nbootsamples==0)) assert(0);

      Int_t ipar = 0;
      for (UInt_t i = 0 ; i < nparam ; i++) {
	if (parfixed[i]==0) {    
	  if (i<nfac) {
	    Setparam(i, fabs(thisparam[ipar]));
	  }
	  else if ((fac_nmix>1)&&(i>=f_nvariance+nfac+1)&&(i<f_nvariance+nfac+1+nfac)) {
	    Setparam(i, fabs(thisparam[ipar]));
	  }
	  else if ((fac_nmix>2)&&(i>=2*(f_nvariance+nfac+1))&&(i<2*(f_nvariance+nfac+1)+nfac)) {
	    Setparam(i, fabs(thisparam[ipar]));
	  }
	  // Force sign on means:
	  else if ((fac_nmix>1)&&(i>=f_nvariance)&&(i<f_nvariance+nfac)) {
	    Setparam(i, fabs(thisparam[ipar]));
	  }
	  else if ((fac_nmix>2)&&(i>=2*f_nvariance+nfac+1)&&(i<2*f_nvariance+nfac+1+nfac)) {
	    Setparam(i, -1.0*fabs(thisparam[ipar]));
	    //if (i==2*f_nvariance+nfac+1) param[i] = -1.0*fabs(thisparam[ipar]));
	  }
	  else {
	    Setparam(i, thisparam[ipar]);
	  }
	  //	  Setparam_err(i] = CalcStdErrLkhdRatio(i);

	  ipar++;
	}
      }
    }

  }
  else cout << "*******MINIMIZATION FAILED!!********" << endl;

  delete [] thisparam;
  delete [] grad;
  delete [] hess;
  
  return status;
}

// Double_t TMinLkhd::CalcStdErrLkhdRatio(UInt_t errparam) {

//   Int_t printlevel = 1;

//   Int_t currentparfixed[nparam];
//   Double_t optparam = param[errparam];
//   for (UInt_t i = 0 ; i < nparam; i++) {
//     currentparfixed[i] = parfixed.at(i);
//     if (i!=errparam) FixPar(i);
//   }
  
//   //Calculate number of fixed parameters and make a lookup table of freeparameters

//   vector<Int_t> revfreeparlist(nparam,-1);
//   Int_t nfixed = 0;
//   UInt_t ifreepar = 0;

//   for (UInt_t i = 0 ; i < nparam; i++) {
//     if (parfixed.at(i)==0) {
//       revfreeparlist.at(i) = ifreepar;
//       ifreepar++;
//     }
//     nfixed += parfixed.at(i);
//   }


//   vector<UInt_t> freeconstraints;
//   for (UInt_t i = 0 ; i < nparam; i++) {

//     // Is this parameter constrained?
//     if (parconstrained.at(i)!= -1) {

//       // Check to see if both of the constrained parameters are free in this estimation
//       if ((parfixed.at(i)==0)&&(parfixed.at(parconstrained.at(i))==0)) {
// 	  // We have a constraint!
// 	  freeconstraints.push_back(i);
//       }
//     }
//   }


//   if (printlevel>0) printf("Calculating standard errors via likelihood ratio: Minimizing %d parameters with %d constraints and while fixing %d.\n",nparam-nfixed,Int_t(freeconstraints.size()),nfixed);
//   //  cout << "Super-special counter of free parameters says there are: " << nfreeparam << ".\n";
//   Int_t nparam_min = nparam - nfixed;



//     /*---- CREATE A NEW KNITRO SOLVER INSTANCE. */

//   int status;
  
//   KTR_context_ptr  kc;
//   kc = KTR_new();
//   if (kc == NULL)
//     {
// 	cout << "*** KTR_new failed, maybe a license issue?\n";
// 	exit( EXIT_FAILURE );
//     }
  
//   status = KTR_load_param_file (kc, "knitro.opt");

//   if (printlevel==0) KTR_set_int_param (kc, KTR_PARAM_OUTLEV, KTR_OUTLEV_NONE);
//   //  KTR_set_int_param (kc, KTR_PARAM_GRADOPT, KTR_GRADOPT_CENTRAL);
//   KTR_set_int_param (kc, KTR_PARAM_HESSOPT, KTR_HESSOPT_LBFGS);

//   if (status != 0)
//     {
//       cout << "*** KTR_load_param_file() returned " << status << "\n";
//       return( false );
//     }
//   else  if (printlevel>0) cout << "**knitro options set\n";

//     int       _nN;
//     double *  _daXInit = NULL;
//     double *  _daXLo = NULL;
//     double *  _daXUp = NULL;
//     int       _nM;
//     int    *  _naCType = NULL;
//     double *  _daCLo = NULL;
//     double *  _daCUp = NULL;
//     int       _nNnzJ;
//     int    *  _naJacIndexVars = NULL;
//     int    *  _naJacIndexCons = NULL;
//     int       _nNnzH;
//     int    *  _naHessCols = NULL;
//     int    *  _naHessRows = NULL;


//     _nN = nparam_min; // number of parameters
//     _nM = freeconstraints.size(); // number of constraints
//     _nNnzJ = 2*_nM;

//     //---- VARIABLE BOUNDS
//     _daXLo  = new double[_nN];
//     _daXUp  = new double[_nN];

//     //---- PUT THE CONSTANT TERM IN THE RIGHT-HAND SIDE.
//     if (_nM>0) {
//       _naCType  = new int[_nM];
//       _daCLo    = new double[_nM];
//       _daCUp    = new double[_nM];
//       //---- SPECIFY THE CONSTRAINT JACOBIAN SPARSITY STRUCTURE.
//       _naJacIndexVars = new int[_nNnzJ];
//       _naJacIndexCons = new int[_nNnzJ];
      
//       for (int i = 0 ; i < _nM; i++) {
// 	_daCLo[i] = 0.0;
// 	_daCUp[i] = 0.0;
// 	_naCType[i] = KTR_CONTYPE_LINEAR;
	
// 	_naJacIndexCons[2*i] = i;  
// 	_naJacIndexVars[2*i] = revfreeparlist.at(freeconstraints.at(i));
	
// 	_naJacIndexCons[2*i+1] = i;  
// 	_naJacIndexVars[2*i+1] = revfreeparlist.at(parconstrained.at(freeconstraints.at(i)));
//       }
//     }

//     //---- SPECIFY THE HESSIAN OF THE LAGRANGIAN SPARSITY STRUCTURE.
//     //    _nNnzH = 0;
//     _nNnzH = nparam_min*(nparam_min+1)/2;
//     _naHessRows = new int[_nNnzH];
//     _naHessCols = new int[_nNnzH];
//     Int_t counter = 0;
//     for (Int_t irow = 0 ; irow < nparam_min ; irow++) {
//       for (Int_t icol = irow ; icol < nparam_min ; icol++) {
// 	_naHessRows[counter] = irow;  _naHessCols[counter] = icol;
// 	counter++;
//       }
//     }

//     //---- INITIAL GUESS FOR x AND lambda.
//     _daXInit = new double[_nN];
//     double * daLambdaInit = NULL;
// //     double *  daLambdaInit = new double[_nM + _nN];
// //     for (int i = 0; i < _nM + _nN; i++)
// //         daLambdaInit[i] = 0.0;


//   // Tell knitro about the parameter values and bounds:
//   Int_t ipar = 0;
//   for (UInt_t i = 0 ; i < nparam ; i++) {
//     if (parfixed[i]==0) {

//      if ((fac_corr!=0) && (i>=nfac) && (i<f_nvariance)) {
// 	_daXLo[ipar] = -1.0;
// 	_daXUp[ipar] = 1.0;
//         _daXInit[ipar] = 0.95*param[i];
//      }
//      // set limits for covariances for second mixture
//      else if ((fac_nmix>1) && (fac_corr!=0) && (i>=f_nvariance+nfac+1+nfac) && (i<f_nvariance+nfac+1+f_nvariance)) {
// 	_daXLo[ipar] = -1.0;
// 	_daXUp[ipar] = 1.0;
//         _daXInit[ipar] = 0.95*param[i];
//      }
//      // set limits for covariances for third mixture
//      else if ((fac_nmix>2) && (fac_corr!=0) && (i>=2*(f_nvariance+nfac+1)+nfac) && (i<2*(f_nvariance+nfac+1)+f_nvariance)) {
// 	_daXLo[ipar] = -1.0;
// 	_daXUp[ipar] = 1.0;
//         _daXInit[ipar] = 0.95*param[i];
//      }
//      else {
//         _daXLo[ipar] = -KTR_INFBOUND;
//         _daXUp[ipar] = KTR_INFBOUND;
// 	_daXInit[ipar] = 0.95*param[i];
//      }
//       ipar++;
//     }
//   }


//   double * jac = NULL;
//   double * c = NULL;

//   if (_nM>0) {
//     c = new double[_nM];
//     jac = new double[_nNnzJ];

//     for (int i = 0 ; i < _nM; i++) {
//       jac[2*i] = 1.0;
//       jac[2*i+1] = -1.0;
//       printf("(%d, %d): %f\n",_naJacIndexCons[2*i], _naJacIndexVars[2*i], jac[2*i]);
//       printf("(%d, %d): %f\n",_naJacIndexCons[2*i+1], _naJacIndexVars[2*i+1], jac[2*i+1]);
//     }
//   }
//   double * daLambda = new double[_nN+_nM];

//   Double_t fvalue = 0;
//   Double_t * thisparam = new Double_t[_nN];
//   Double_t * grad = new Double_t[_nN];
//   Double_t * hess = new Double_t[_nNnzH];


//   ifreepar = 0;
//   for (UInt_t ipar = 0 ; ipar < nparam ;ipar++) {
//     if (parfixed[ipar]==0) {
//       thisparam[ifreepar] = param[ipar];
//       ifreepar++;
//     }
//   }

 
//   //  cout << "ok evaluating likelihood now!\n";
//   //  for (UInt_t ipar = 0 ; ipar < nparam ;ipar++) thisparam[ipar] = param[ipar];

//   status = KTR_init_problem (kc, _nN,
// 			     KTR_OBJGOAL_MINIMIZE, KTR_OBJTYPE_GENERAL,
// 			     _daXLo, _daXUp,
// 			     _nM, _naCType, _daCLo, _daCUp,
// 			     _nNnzJ, _naJacIndexVars, _naJacIndexCons,
// 			     _nNnzH, _naHessRows, _naHessCols,
// 			     _daXInit, daLambdaInit);
//   if (status != 0)
//     {
//       cout << "*** KTR_init_problem() returned " << status << "\n";
//       return( status);
//     }
//   else if (printlevel>0) cout << "**knitro initialized\n";


//   delete [] _daXLo;
//   delete [] _daXUp;
//   if (_nM>0) {
//     delete [] _naCType;
//     delete [] _daCLo;
//     delete [] _daCUp;
//     delete [] _naJacIndexVars;
//     delete [] _naJacIndexCons;
//   }
//   delete [] _naHessRows;
//   delete [] _naHessCols;
//   //  delete [] daLambdaInit;
//   delete [] _daXInit;

//   LkhdFcn(_nN,grad,fvalue,thisparam,1,hess);
//   double optlikelihood = fvalue; 

//   printf("Calculating stderr for param %d=%f\n",errparam,thisparam[0]);

//   status =1 ;

//   LkhdFcn(_nN,grad,fvalue,thisparam,status,hess);
//   printf("lkhd=%f for param %d=%f\n",fvalue,errparam,thisparam[0]);

//   counter = 0;
//   Int_t badstatuscounter = 0;
//   while ((status==1)||(status==2)||(status==3)) {
//     if (counter>0) {

//       LkhdFcn(_nN,grad,fvalue,thisparam,status,hess);

//       for (int i = 0 ; i < _nM ; i++) {
// 	c[i] = thisparam[revfreeparlist.at(freeconstraints.at(i))] 
// 	  - thisparam[revfreeparlist.at(parconstrained.at(freeconstraints.at(i)))];
//       }
//     }

//     double lkhd = fvalue;
//     double lkhdgrad = grad[0];    

//     fvalue = (lkhd - optlikelihood - 0.5)*(lkhd - optlikelihood - 0.5);
//     if (status>=2) grad[0] = 2.0*(lkhd - optlikelihood - 0.5)*lkhdgrad;
//     if (status==3) hess[0] = 2.0*lkhdgrad*lkhdgrad + 2.0*(lkhd - optlikelihood - 0.5)*hess[0];

//     printf("iter %d: lkhd=%f for param %d=%f\n",counter,(lkhd - optlikelihood - 0.5),errparam,thisparam[0]);

//     status = KTR_solve (kc, thisparam , daLambda, 0, &fvalue,
//   			c, grad, jac, hess, NULL, NULL);    


// //     if (counter>0) LkhdFcn(_nN,grad,fvalue,thisparam,status,NULL);
// //     status = KTR_solve (kc, _daXInit , daLambda, 0, &fvalue,
// //  			c, grad, jac, NULL, NULL, NULL);    
//     counter++;
//     if (status==-100) {
//       printf("KNITRO FAILED WITH STATUS %d, restart %d...\n",status,badstatuscounter);
//       Double_t nudge = 0.001;
//       if (badstatuscounter==1) nudge = -0.001;
//       else if (badstatuscounter==2) nudge = 0.01;
//       else if (badstatuscounter==3) nudge = -0.01;
//       else if (badstatuscounter==4) nudge =  0.1;
      
//       if (badstatuscounter<5) {
// 	for (Int_t ipar = 0 ; ipar < _nN ;ipar++) thisparam[ipar] += nudge;
// 	status= KTR_restart(kc,thisparam,daLambda);
// 	status = 1;
// 	badstatuscounter++;
//       }
//     }
//   }

//   Double_t stderr = 0.0;
//   if (status==0) {
//     if (printlevel>0) cout << "Minimization done. Lkhd= " << fvalue << "\n";
//     stderr = fabs(thisparam[0] - optparam);
//   }
//   else cout << "*******MINIMIZATION FAILED!!********" << endl;

//   param[errparam] = optparam;

//   KTR_free (&kc);
  
//   delete [] jac;
//   delete [] daLambda;
//   delete [] c;
//   delete [] thisparam;
//   delete [] grad;
//   delete [] hess;
  
//   SetFixPar(currentparfixed);

//   return stderr;

// }


void TMinLkhd::EvalLkhd(Double_t & logLkhd, Double_t * gradL, Double_t * hessL, Double_t *gamma, Int_t iflag, Int_t rank, Int_t np) {

  //Calculate hessian using factor scores
  if ((HessFactorScores)&&(iflag>2)) {

    // First calculate "exact" likelihood and gradient
    CalcLkhd(logLkhd, gradL, hessL, gamma, 2, rank, np);

    //Store values
    Double_t tmplogLkhd = logLkhd;
    vector<Double_t> tmpgradL(nparam,0.0);
    for (UInt_t ipar = 0 ; ipar<nparam; ipar++) tmpgradL[ipar] = gradL[ipar];

    //calculate approximate Hessian
    Double_t current_threshold = adapt_int_thresh;
    SetAdaptIntThresh(10000.0);
    CalcLkhd(logLkhd, gradL, hessL, gamma, iflag, rank, np);

    // Set threshold, likelihood and gradient back to "exact" levels
    SetAdaptIntThresh(current_threshold);
    logLkhd = tmplogLkhd;
    for (UInt_t ipar = 0 ; ipar<nparam; ipar++) gradL[ipar] = tmpgradL[ipar]; 
  }
  else {
    CalcLkhd(logLkhd, gradL, hessL, gamma, iflag, rank, np);
  }
}

void TMinLkhd::CalcLkhd(Double_t & logLkhd, Double_t * gradL, Double_t * hessL, Double_t *gamma, Int_t iflag, Int_t rank, Int_t np) {

  
  //   cout << "point 1\n";
  // **********************
  // Do we really want to set param for every evaluation??
  // **********************
  // First force the sign of certain parameters
  // (means of mixtures and variances of factors)

  vector<Int_t> freeparlist;
  freeparlist.resize(nfreeparam);
  
  // If predicting then there are no free parameters
  UInt_t ifreepar = 0;
  for (UInt_t i = 0 ; i < nparam ; i++) {
    //    if (rank>0) cout << "processing parameter #" << i << endl;
    if (parfixed[i]==0) {
      //first force variances positive:
      if (i<nfac) {
	Setparam(i, fabs(gamma[ifreepar]));
      }
      else if ((fac_nmix>1)&&(i>=f_nvariance+nfac+1)&&(i<f_nvariance+nfac+1+nfac)) {
	Setparam(i, fabs(gamma[ifreepar]));
      }
      else if ((fac_nmix>2)&&(i>=2*(f_nvariance+nfac+1))&&(i<2*(f_nvariance+nfac+1)+nfac)) {
	Setparam(i, fabs(gamma[ifreepar]));
      }
      // Force sign on means:
      else if ((fac_nmix>1)&&(i>=f_nvariance)&&(i<f_nvariance+nfac)) {
	Setparam(i, fabs(gamma[ifreepar]));
                                           // 	if (i==f_nvariance) param[i] = fabs(gamma[ifreepar]);
	// 	else param[i] = -1.0*fabs(gamma[ifreepar]);
      }
      else if ((fac_nmix>2)&&(i>=2*f_nvariance+nfac+1)&&(i<2*f_nvariance+nfac+1+nfac)) {
	Setparam(i, -1.0*fabs(gamma[ifreepar]));
	// 	if (i==(2*f_nvariance+nfac+1)) param[i] = -1.0*fabs(gamma[ifreepar]);
	// 	else fabs(gamma[ifreepar]);
      }
      else Setparam(i, gamma[ifreepar]);

      freeparlist[ifreepar] = i;

      ifreepar++;
    }
  }


  //  PrintParam();
//   if ((nfac>1)&&(fac_corr!=0)) {
//     for (Int_t imix = 0 ; imix < fac_nmix ; imix++) {
//       Int_t first_corr = imix*f_nvariance+nfac;
//       for (UInt_t i = nfac ; i < f_nvariance ; i++) {
// 	if (param[i]>0.99) {
// 	  cout << "Correlation close to 1.0! Likelihood found a fake min!\n";
// 	  cout << "parameter " << i << "=" << param[i] << endl;
// 	  assert(0);
// 	}
//       }
//     }
//   }

// if (rank>0)   cout << "point 2, rank=" << rank << endl;
  // cout << "point 2\n";

  //Skip integration over factor distribution if predicting or initializing
  UInt_t nfac_eff = nfac;
  UInt_t ntyp_eff = ntyp;
  if ((initializing)||(predicting)) {
    nfac_eff = 0;
    ntyp_eff = 0;
  }
  // if (rank>0)   cout << "point 3, rank=" << rank << endl;
  // cout << "point 3\n";

  vector <Double_t> w_mix;
  //  cout << "Resizing vector from " <<w_mix.size() << " to " << fac_nmix << ".\n";
  w_mix.resize(est_nmix);

  vector <Double_t> fac_mean;
  //  fac_mean.resize(nfac*est_nmix);  // mean, nasty, evil BUG!!

  // cout << "point 3.02\n";

  //Calculate the denominator for the mixture weights
  Double_t w_mix_den = 1.0;
  if (est_nmix>1) {
    for (UInt_t imix = 0 ; imix < est_nmix-1 ; imix++) {
      w_mix_den += exp(Getfw(imix));
    }
  }
  // cout << "point 3.1\n";
 
  // Calculate actual weights and means for different mixtures:
  for (UInt_t imix = 0 ; imix < est_nmix*(nfac_eff>0)+(nfac_eff==0) ; imix++) {
    //    cout << "point 3.1" << imix << "\n";
    //Calculate weight for this mixture (logit-like scheme)
    w_mix[imix] = 1.0;
    if ((est_nmix>1)&&(nfac_eff>0)) { 
      w_mix[imix] = 1.0/w_mix_den;
      if (imix<est_nmix-1) w_mix[imix] = exp(Getfw(imix))*w_mix[imix];
    }
    
    //Get means for this mixture
    if ((est_nmix==1)||(nfac_eff==0)) {
      for (UInt_t ifac = 0 ; ifac < nfac; ifac++) fac_mean.push_back(0.0);
    }
    else if (imix < est_nmix-1) {
      for (UInt_t ifac = 0 ; ifac < nfac; ifac++) fac_mean.push_back(Getfmean(imix,ifac));
    }
    else {
      // Calculate mean of first mixture such that mean of all mixtures = 0.0:
      for (UInt_t ifac = 0 ; ifac < nfac; ifac++) {
	Double_t weighted_sum = 0.0;
	for (UInt_t imix_mean = 0 ; imix_mean < est_nmix-1 ; imix_mean++) {
	  weighted_sum += -exp(Getfw(imix_mean))*Getfmean(imix_mean,ifac);
	}
	fac_mean.push_back(weighted_sum);
      }
    }
  }
//    cout << "finished calculating means and weights:\n";
//    cout << "fac_mean size=" << fac_mean.size() << " and w_mix size=" << w_mix.size() << endl;
//    for (UInt_t imix = 0 ; imix < est_nmix*(nfac_eff>0)+(nfac_eff==0) ; imix++) {
//      cout << "imix=" << imix << ", mu1=" << fac_mean[imix*nfac] 
//  	 << ", mu2=" << fac_mean[imix*nfac+1]
//  	 << ", weight=" << w_mix[imix] << ".\n";

//    }
  // cout << "point 3.2\n";


// Find the integration points for each mixture for bivariate normal case with correlation
//  const Double_t pi = TMath::Pi();
  logLkhd = 0.0;
  vector <Double_t> x2sc, dx2sc_drho, x2;

  if ((nfac_eff==2)&&(fac_corr!=0)) {
    x2sc.resize(est_nmix*nquad_points*nquad_points);
    x2.resize(est_nmix*nquad_points*nquad_points);
    dx2sc_drho.resize(est_nmix*nquad_points*nquad_points);

    for (UInt_t imix = 0 ; imix < est_nmix ; imix++) {
      for (UInt_t i = 0; i <  nquad_points; i++) {

	for (UInt_t j = 0; j <  nquad_points; j++) {
	  x2sc[Getx2i(imix,i,j)] = Getfvar(imix,1)*
	    (Getfvar(imix,2)*GetHGQx(nquad_points,i)+sqrt(1.0-Getfvar(imix,2)*Getfvar(imix,2))*GetHGQx(nquad_points,j));
	  
	  x2[Getx2i(imix,i,j)] = Getfvar(imix,2)*GetHGQx(nquad_points,i) + sqrt(1.0-Getfvar(imix,2)*Getfvar(imix,2))*GetHGQx(nquad_points,j);
	  
	  dx2sc_drho[Getx2i(imix,i,j)] = Getfvar(imix,1)*
	    (GetHGQx(nquad_points,i)-Getfvar(imix,2)*GetHGQx(nquad_points,j)/(sqrt(1.0-Getfvar(imix,2)*Getfvar(imix,2))));
	}
      }
    }
  }

  // if (rank>0)   cout << "point 4, rank=" << rank << endl;
  // cout << "point 4\n";
  // Reset gradient and hessian vectors
  vector <Double_t> totalgrad, gradilk, totalhess, hessilk;

  if (iflag>=2) {
    UInt_t numberparam = nfreeparam;
    if (predicting==1) numberparam = nfac;

    for (ifreepar = 0 ; ifreepar < numberparam; ifreepar++) {
      gradL[ifreepar] = 0;
    }
    
    totalgrad.resize(nparam);
    gradilk.resize(nparam);
    
    if (iflag==3) {
      for (UInt_t ihess = 0; ihess < (numberparam*(numberparam+1)/2) ; ihess++) {
	hessL[ihess] = 0.0;
      }

      totalhess.resize(nparam*nparam);
      hessilk.resize(nparam*nparam);
    }

//     ifreepar = 0;
//     for (UInt_t ipar = 0 ; ipar < nparam; ipar++) {
//       if (parfixed[ipar]==0) {
// 	gradL[ifreepar] = 0;
// 	ifreepar++;
//       }
//     }

//     ifreepar = 0;
//       for (UInt_t ipar1 = 0 ; ipar1 < nparam; ipar1++) {
// 	for (UInt_t ipar2 = ipar1 ; ipar2 < nparam; ipar2++) {
// 	  if ((parfixed[ipar1]==0)&&(parfixed[ipar2]==0)) {
// 	    hessL[ifreepar] = 0;
// 	    ifreepar++;
// 	  }
// 	}
//       }

  }

// This is the calculation for dividing up the observations between processors:
  UInt_t div,start,end;
  if (predicting==0) {
    div = nobs/np;
    start = rank*div;
    end = (rank+1)*div;
    if (rank==np-1) {
      end = nobs;
    }
  }
  else {
    // If predicting only minimize for one observation
    start = predictobs;
    end = start+1;
  }
  //  cout << start <<", " << end << "\n";


  // Find first outcome model
  UInt_t first_outmodel=9999;
  UInt_t estimating_outcomes = 0;
  for (UInt_t imod = 0 ; imod < models.size() ; imod++) {
    // Find first outcome model
    if ((models[imod].GetGroup()==1)&&(first_outmodel==9999)) first_outmodel = imod;

    // Use the fact that all outcome models are ignored when estimating measurement system
    if ((models[imod].GetGroup()>0)&&(models[imod].GetIgnore()==0)) estimating_outcomes = 1;
  }

  //comment this out later:
  //  cout << "Estimating_outcomes=" << estimating_outcomes << endl;
  
  // if (rank>0)   cout << "point 5, rank=" << rank << endl;
  //  cout << "point 5\n";
  //Declare some variables to save time:
  vector<Double_t> modhess;
  vector<Double_t> modelEval;

  // Loop over observations
//   for (UInt_t iobs = start ; iobs < end ; iobs++) {
//     UInt_t i = bootstrapobs[iobs];
  TRandom *r3 = new TRandom3();
  Int_t countstochflag1 = 0;

  for (UInt_t i = start ; i < end ; i++) {

    Int_t stochflag = iflag;
    if ((iflag>2)&&(stoch_deriv_frac>0.0)) {
      Double_t draw = r3->Uniform();
      if (draw>stoch_deriv_frac) {
	stochflag = 2;
	countstochflag1++;
      }
    }
    
    //Check to see if observation is missing for relevant outcome models
    UInt_t estimate_thisobs = 0;
    if (estimating_outcomes==0) estimate_thisobs = 1;
    else {
      for (UInt_t imod = first_outmodel ; imod < models.size() ; imod++) {
	if (models[imod].GetIgnore()==0) {
	  if (models[imod].GetMissing()>-1) {
	    if (data[i*nvar + models[imod].GetMissing()] != 0) estimate_thisobs = 1;
	  }
	  else estimate_thisobs = 1;
	}
      }
    }
    if ( (estimate_thisobs==1) && (skipobs[i]==0) &&  ((bootstrapobs[i]>0) | (predicting==1)) ) {
      //   for (UInt_t i = 0 ; i < 1 ; i++) {
      //Reset i quantities
      Double_t totalprob = 0.0;
      if (stochflag>=2) {
	for (UInt_t ipar = 0 ; ipar < nparam; ipar++) totalgrad[ipar] = 0.0;
	//      if (stochflag==3) for (UInt_t ipar = 0 ; ipar < nparam*nparam; ipar++) totalhess[ipar] = 0.0;
	if (stochflag==3) {
	  
	  if (predicting==0) {
	    for (ifreepar = 0 ; ifreepar < nfreeparam; ifreepar++) {
	      for (UInt_t jfreepar = ifreepar ; jfreepar < nfreeparam; jfreepar++) {
		Int_t fullHessindex = freeparlist[ifreepar]*nparam+freeparlist[jfreepar];
		totalhess[fullHessindex] = 0.0;
	      }
	    }
	  }
	  else {
	    for (UInt_t ifac = 0 ; ifac < nfac; ifac++) {
	      for (UInt_t jfac = ifac ; jfac < nfac; jfac++) {
		Int_t fullHessindex = ifac*nparam+jfac;
		totalhess[fullHessindex] = 0.0;
	      }
	    }
	  }
	}
      }

      //      cout << "point 6\n";
      for (UInt_t imix = 0 ; imix < est_nmix*(nfac_eff>0)+(nfac_eff==0) ; imix++) {

	// Store: npoints, factor_values and weights for each point
	UInt_t nint_points = 1;
	nint_points = pow(nquad_points,nfac_eff);
	//       std::vector<UInt_t> fac_npoints(nfac,nquad_points);

	if (fac_npoints.size()==0) fac_npoints.resize(nfac,nquad_points);

	std::vector<UInt_t> thisobs_nquadpoints(fac_npoints);

	//This stores the location in the array of the correct weights and integration points
	//       std::vector<Int_t> facint_offset(nfac,0);
	
	if ((adapt_int==1)&&(nfac_eff>0)) {
	  Double_t ave_std = 0.0;
	  Int_t nfac_noRE = 0;
	
	  for (UInt_t ifac = 0 ; ifac <nfac ; ifac++) {

	    //If random effect (fstderr<0) then use default

	    if (fstderr[i][ifac]>0.0) {
	      // Ad-hoc rule, probably needs to be adjusted depending on the model
	      // Odd number of integration points indicates importance sampling
	      thisobs_nquadpoints.at(ifac) = 1+2*floor(fstderr[i][ifac]/Getfvar(imix,ifac)/adapt_int_thresh);

	      //If more points than default are specified, use default
	      if (thisobs_nquadpoints.at(ifac) > fac_npoints.at(ifac)) thisobs_nquadpoints.at(ifac) = fac_npoints.at(ifac);

	      //keep track of total uncertainty in location of factor
	      ave_std += fstderr[i][ifac]/Getfvar(imix,ifac);
	      nfac_noRE++;
	    }
	  }

	  //If std is big on average then use default, otherwise use adaptive points
	  if ( ave_std/nfac_noRE > 1.5) {
	    for (UInt_t ifac = 0 ; ifac <nfac ; ifac++) {
	      thisobs_nquadpoints.at(ifac) = fac_npoints.at(ifac);
	      //	      if (thisobs_nquadpoints[ifac]<fac_npoints[ifac]) fac_npoints[ifac] = tmp_npoints[ifac];
	    }
	  }

	  // Get offset for each case:
	  // if ((fac_npoints[ifac]==1)||(fac_npoints[ifac]==nquad_points)) facint_offset[ifac] = 0;
	  // if (fac_npoints[ifac]==3) facint_offset[ifac] = 1;
	  // else if (fac_npoints[ifac]==5) facint_offset[ifac] = 4;
	  // else if (fac_npoints[ifac]==7) facint_offset[ifac] = 9;

	  // Calculate total number of integration points
	  nint_points = 1;
	  for (UInt_t ifac = 0 ; ifac <nfac ; ifac++) nint_points *= thisobs_nquadpoints[ifac];

	}

	
	//	cout << "*********Finished calculating adaptive integration numbers! " << nint_points << " integration points for obs " << Int_t(data[i*nvar+indexvar]) << ".\n";

	//loop over types (ntyp is number of types - 1)
	for (UInt_t itype = 0; itype <= ntyp_eff; itype++) {

	  //loop over integration points
	  vector<UInt_t> facint(nfac,0);
	  for (UInt_t intpt = 0; intpt < nint_points ; intpt++) {
	    
	    // get integration points for each factor from loop
	    if (intpt!=0) {
	      for (UInt_t ifac = 0 ; ifac <nfac_eff ; ifac++) {
		facint[ifac]++;
		if (facint[ifac]<thisobs_nquadpoints[ifac]) break;
		else facint[ifac] = 0;
	      }
	    }
	    //	cout << "****calculating element imix=" << imix << ", intpt=" << intpt << "\n";
	    // Reset ilk quantities (i=obs, l=mixture, k=integration point)
	    Double_t probilk = 1.0;

	  if (stochflag>=2) {

	    // USE ASSIGN TO RESET TO 0?????
	    for (UInt_t ipar = 0 ; ipar < nparam; ipar++) gradilk[ipar] = 0.0;

 	    if (stochflag==3) {

	      if (predicting==0) {
		for (ifreepar = 0 ; ifreepar < nfreeparam; ifreepar++) {
		  for (UInt_t jfreepar = ifreepar ; jfreepar < nfreeparam; jfreepar++) {
		    Int_t fullHessindex = freeparlist[ifreepar]*nparam+freeparlist[jfreepar];
		    hessilk[fullHessindex] = 0.0;
		  }
		}
	      }
 	      else {
 		for (UInt_t ifac = 0 ; ifac < nfac; ifac++) {
 		  for (UInt_t jfac = ifac ; jfac < nfac; jfac++) {
 		    Int_t fullHessindex = ifac*nparam+jfac;
 		    hessilk[fullHessindex] = 0.0;
 		  }
 		}
 	      }
	    }
	  }

	  //	  cout << "calculating distribution part\n";
	  //Factor distribution part of likelihood
	  if (nfac_eff>0) {

	    //weight from mixture
	    probilk = probilk*w_mix[imix];

	    //weight of quadrature approximation to integral
	    for (UInt_t ifac = 0 ; ifac < nfac_eff; ifac++) {
	      probilk = probilk * GetHGQw(thisobs_nquadpoints[ifac],facint[ifac]);
	      //	      printf("%8.5f; ",GetHGQw(thisobs_nquadpoints[ifac],facint[ifac]));
	    }
	    //	    printf("\n");
	    /*
	    if (adapt_int==0) {
	      for (UInt_t ifac = 0 ; ifac < nfac_eff; ifac++) {
		probilk = probilk *w[facint[ifac]];
	      }
	    }
	    else {
	      //	      cout << "Calculating quadrature weights!\n";
	      for (UInt_t ifac = 0 ; ifac < nfac_eff; ifac++) {
		probilk = probilk * GetHGQw(fac_npoints[ifac],facint[ifac]);
		  //		if (fac_npoints[ifac]%2==1) probilk = probilk * wadapt[facint[ifac] + facint_offset[ifac]];
		  //		else probilk = probilk *w[facint[ifac]];	  
	      }
	    }
	    */
	    
	    if ((stochflag>=2)&&(stage!=2)) {
	      //  for (UInt_t ifac = 0; ifac < nfac ; ifac++) gradilk[Getifvar(imix,ifac)] += 1.0/Getfvar(imix,ifac);
 	      // gradient of weight term:
 	      for (UInt_t imix_grad = 0; imix_grad < est_nmix-1 ; imix_grad++) {
 		if (imix_grad==imix) gradilk[Getifw(imix_grad)] += 1.0-w_mix[imix_grad];
 		else gradilk[Getifw(imix_grad)] += -w_mix[imix_grad];
 	      }
 	      //Hessian of weight terms:
 	      if (stochflag==3) {
 		for (UInt_t imix_grad1 = 0; imix_grad1 < est_nmix-1 ; imix_grad1++) {
 		  for (UInt_t imix_grad2 = imix_grad1; imix_grad2 < est_nmix-1 ; imix_grad2++) {
 		    if (imix_grad1==imix_grad2) {
 		      hessilk[Getifw(imix_grad1)*nparam+Getifw(imix_grad2)] += -w_mix[imix_grad1]*(1.0-w_mix[imix_grad2]);
 		    }
 		    else {
 		      hessilk[Getifw(imix_grad1)*nparam+Getifw(imix_grad2)] += w_mix[imix_grad1]*w_mix[imix_grad2];
 		    }
 		  }
 		}
 	      }
	    }
  	  }

	  //	  printf("After weights prob=%8.5f\n",probilk);
	  //	  cout << "calculating model part\n";
	  //Model part of likelihood
 	  UInt_t firstpar = nfac_param;

	  //Calculate value of factors at this integration point
 	  vector<Double_t> fac_val;
 	  if (initializing) {
	    if ((initEstOutcomeLoadings) && (estimating_outcomes==1)) {
	      for (UInt_t ifac = 0; ifac < nfac ; ifac++) fac_val.push_back(fscore[i][ifac]);
	    }
	    else {
	      for (UInt_t ifac = 0; ifac < nfac ; ifac++) fac_val.push_back(0.0);
	    }
 	  }
 	  else if (predicting) {
	    for (UInt_t ifac = 0; ifac < nfac ; ifac++) fac_val.push_back(gamma[ifac]);
	  }
 	  else if (adapt_int==0) {
	    if (fac_corr==0) {
	      for (UInt_t ifac = 0; ifac < nfac ; ifac++) {
		fac_val.push_back(Getfvar(imix,ifac)*GetHGQx(thisobs_nquadpoints[ifac],facint[ifac])+fac_mean[ifac+nfac*imix]);
 	      }
	    }
	    else {
	      if (nfac>0) fac_val.push_back(Getfvar(imix,0)*GetHGQx(nquad_points,facint[0])+fac_mean[0+nfac*imix]);
	      if (nfac>1) fac_val.push_back(x2sc[Getx2i(imix,facint[0],facint[1])]+fac_mean[1+nfac*imix]);
	    }
 	  }
	  else if (adapt_int==1) {
	    for (UInt_t ifac = 0; ifac < nfac ; ifac++) {

	      //If npoints is odd, then we use importance sampling
	      if (thisobs_nquadpoints[ifac]%2==1) {
		  fac_val.push_back(fscore[i][ifac] + fstderr[i][ifac]*GetHGQx(thisobs_nquadpoints[ifac],facint[ifac]));

		  // Additional terms from importance sampling:

		  probilk = probilk * exp(GetHGQx(thisobs_nquadpoints[ifac],facint[ifac])*GetHGQx(thisobs_nquadpoints[ifac],facint[ifac])/2.0);
		  probilk = probilk * exp(-1.0*fac_val[ifac]*fac_val[ifac]/(2*Getfvar(imix,ifac)*Getfvar(imix,ifac)));
		  
		}
		else {
		  fac_val.push_back(Getfvar(imix,ifac)*GetHGQx(fac_npoints[ifac],facint[ifac])+fac_mean[ifac+nfac*imix]);
		}
	    }
	  }

	  //Add type to factor array
	  vector<Double_t> type(ntyp,0.0);
	  if ((ntyp_eff>0)&&(itype>0))  {
	    type[itype-1] = 1.0;
	  }
	  fac_val.insert(fac_val.end(),type.begin(),type.end());	

	  
	  //	  printf("After importance: prob=%8.5f\n",probilk);
 	  // // printf("%5d %10d, %3d/%3d: ", i,Int_t(data[i*nvar+indexvar]), intpt, nint_points);
	  //	  for (UInt_t ifac = 0; ifac < nfac ; ifac++) printf("%4.5f; ", fac_val.at(ifac));
	  //	  printf("\n");

	  // printf("\n");
	  //	  if (i>100) assert(0);
	  //	  cout << "now calculating actual models\n";
 	  for (UInt_t imod = 0 ; imod < models.size() ; imod++) {
	    //ADD IF MISSING HERE???
	    if (models[imod].GetIgnore()==0) {

	      //	      if ((!initializing)&&(intpt==0)&&(i==0)&&(imix==0)) cout << "Processing model " << models[imod].GetName() << " of type " <<  models[imod].GetType() << " ignore=" << models[imod].GetIgnore() << "\n";

	      //Do not estimate gradient and hessian for measurement system if evaluating outcomes
	      Int_t thisflag = stochflag;
	      if ((imod<first_outmodel)&&(estimating_outcomes==1)) thisflag = 1;
	      
	      models[imod].Eval(i*nvar,data,param,firstpar,fac_val,modelEval,modhess,thisflag);
//	      models[imod].Eval(i*nvar,data,param,firstpar,fac_val,modelEval,modhess,2);
 	      probilk = probilk*modelEval[0];

//	      if (imod==type_model) {
//		printf("Evalulated Type model with prob=%8.4f, where types are:",modelEval[0]);
//		for (UInt_t it = 0; it < ntyp; it++) {
//		  printf("%2.1f, ",fac_val[nfac+it]);
//		}
//		printf("\n");
//	      }
	      
	      if (thisflag>=2) {
 		//gradients for factor-specific parameters
		if (nfac_eff>0) {
		  if ((fac_corr!=0)&&(nfac>1)) {
		    gradilk[Getifvar(imix,2)] += modelEval[2]*dx2sc_drho[Getx2i(imix,facint[0],facint[1])];
		  }
		  for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {
		    // gradient of factor variance (sigma_theta)
		    //   		  if (ifac==0) gradilk[Getifvar(imix,ifac)] += modelEval[ifac+1]*x[int1];
		    //   		  if (ifac==1) gradilk[Getifvar(imix,ifac)] += modelEval[ifac+1]*x[facint[1]];
		    if ((fac_corr!=0)&&(ifac==1)) {
		      gradilk[Getifvar(imix,ifac)] += modelEval[ifac+1]*x2sc[Getx2i(imix,facint[0],facint[1])]/Getfvar(imix,ifac);
		    }
		    else {
		      gradilk[Getifvar(imix,ifac)] += modelEval[ifac+1]*Getfvar(imix,ifac)*x[facint[ifac]]/Getfvar(imix,ifac);
		    }

		    //gradient for means and weights
		    if (est_nmix>1) {
		      if (imix== est_nmix-1) {
			// last mixture includes all means:
			for (UInt_t imix_grad = 0; imix_grad < est_nmix-1 ; imix_grad++) {
			  gradilk[Getifmean(imix_grad,ifac)] += modelEval[ifac+1]*(-exp(Getfw(imix_grad)));
			  
			  // need to include gradients of weights in last mixture
			  // since they go into the calculation of the last mean
			  gradilk[Getifw(imix_grad)] += modelEval[ifac+1]*(-exp(Getfw(imix_grad))*Getfmean(imix_grad,ifac));
			}
		      }
		      else {
			gradilk[Getifmean(imix,ifac)] += modelEval[ifac+1];
		      }
		    }
		  }
		}
		else if ((predicting==1)&&(nfac_eff==0)) {
		  for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) gradilk[ifac] += modelEval[ifac+1];
		}
	      
		if (predicting==0) {
		  // gradients for model-specific parameters:
		  for (UInt_t iparam = 0 ; iparam < nparam_models[imod] ; iparam++) {
		    gradilk[iparam+firstpar] += modelEval[iparam+1+nfac];
		  }
		  
// 		  for (UInt_t iparam = 0 ; iparam < (nparam_models[imod] + 1 + nfac) ; iparam++) {
// 		    if (isnan(modelEval[iparam])) 
// 		      printf("****EvalLkhd: ISNAN: obs %5d, model %5d,par %5d/%5d. Lkhd=%8.4e\n",i,imod,iparam,(nparam_models[imod] + 1 + nfac),modelEval[0]);
// 		  }
		}

		// Ok, now calculate Hessian if not missing data (modhess.size=0)
		if ((thisflag==3)&&(modhess.size()>0)) {
		  Double_t nDimModHess = nparam_models[imod] + nfac;
		  if (nDimModHess*nDimModHess!=modhess.size()) {
		    cout << "Shit Got Hessian dimensions wrong! Expected " << nDimModHess*nDimModHess
			 << " and got " << modhess.size() << "\n";
		    assert(0);
		  }
		  if (predicting==1) {
		    for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {        //loop over each fac
		      for (UInt_t jfac = ifac ; jfac < nfac ; jfac++) {        //loop over each fac
			Int_t modhessindex = (nDimModHess)*(ifac) + jfac; 
			Int_t fullHessindex = nfac*ifac + jfac;
			hessilk[fullHessindex] += modhess[modhessindex];
		      }
		    }
		  }

		  if (nfac_eff>0) {
		    
		    if ((est_nmix>1)&&(imix== est_nmix-1)) {
		      Int_t fullHessindex;
		      
		      // d/dwdw
		      for (UInt_t imixhes1= 0 ; imixhes1 < est_nmix-1 ; imixhes1++) {    //loop over each mixture
			for (UInt_t imixhes2 = imixhes1 ; imixhes2 < est_nmix-1 ; imixhes2++) {        //loop over each mixture
			  
			  fullHessindex = (nparam)*Getifw(imixhes1)+Getifw(imixhes2);
			  
			  for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {        //loop over each fac
			    for (UInt_t jfac = 0 ; jfac < nfac ; jfac++) {        //loop over each fac
			      // store d/dthetai dthetaj for this model
			      Int_t modhessindex = (nDimModHess)*(ifac) + jfac; 
			      if (jfac<ifac) modhessindex = (nDimModHess)*(jfac) + ifac; // only upper triangle defined
			      Double_t ddthetasq = modhess[modhessindex];
			      
			      Double_t multiplier = (-exp(Getfw(imixhes1))*Getfmean(imixhes1,ifac))*(-exp(Getfw(imixhes2))*Getfmean(imixhes2,jfac));
			      hessilk[fullHessindex] += ddthetasq*multiplier;
			    }
			    
			    if (imixhes1==imixhes2) {
			      hessilk[fullHessindex] += modelEval[ifac+1]*(-exp(Getfw(imixhes1))*Getfmean(imixhes1,ifac));
			    }
			    
			  }
			  
			}
		      }

		      //d/dmudw
		      for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {        //ifac of mu (f')
			for (UInt_t imixhes1= 0 ; imixhes1 < est_nmix-1 ; imixhes1++) {    //imix of mu (n)
			  for (UInt_t imixhes2 = 0 ; imixhes2 < est_nmix-1 ; imixhes2++) {        //imix of omega (m)
			    
			    if (Getifmean(imixhes1,ifac)<Getifw(imixhes2) ) {
			      fullHessindex = (nparam)*Getifmean(imixhes1,ifac)+Getifw(imixhes2);
			    }
			    else fullHessindex = (nparam)*Getifw(imixhes2)+Getifmean(imixhes1,ifac);
			    
			    for (UInt_t jfac = 0 ; jfac < nfac ; jfac++) {        //loop over each fac (f)
			      // store d/dthetai dthetaj for this model
			      Int_t modhessindex = (nDimModHess)*(ifac) + jfac; 
			      if (jfac<ifac) modhessindex = (nDimModHess)*(jfac) + ifac; // only upper triangle defined
			      Double_t ddthetasq = modhess[modhessindex];
			      
			      Double_t multiplier = (-exp(Getfw(imixhes1)))*(-exp(Getfw(imixhes2)));
			      hessilk[fullHessindex] += ddthetasq*multiplier*Getfmean(imixhes2,jfac);
			    }
			    
			    if (imixhes1==imixhes2) {
			      hessilk[fullHessindex] +=  modelEval[ifac+1]*(-exp(Getfw(imixhes2)));
			    }
			  }
			}
		      }
		      
		    }
		  
		    for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {    //loop over each fac
		      
		      //distribution parameters cross distribution parameters:
		      for (UInt_t jfac = 0 ; jfac < nfac ; jfac++) {        //loop over each fac
			
			// store d/dthetai dthetaj for this model
			Int_t modhessindex = (nDimModHess)*(ifac) + jfac; 
			if (jfac<ifac) modhessindex = (nDimModHess)*(jfac) + ifac;
			Double_t ddthetasq = modhess[modhessindex];
			
			// d/dsigmadsigma
			if (jfac >= ifac) {
			  Int_t fullHessindex = (nparam)*Getifvar(imix,ifac)+Getifvar(imix,jfac);

			  if ((fac_corr!=0)&&(jfac==1))  {
			    if ((ifac==0)&&(jfac==1)) hessilk[fullHessindex] += ddthetasq*x[facint[0]]*x2[Getx2i(imix,facint[0],facint[1])];
			    else if ((ifac==1)&&(jfac==1)) hessilk[fullHessindex] += ddthetasq*x2[Getx2i(imix,facint[0],facint[1])]*x2[Getx2i(imix,facint[0],facint[1])];;
			  }
			  else {
			    hessilk[fullHessindex] += ddthetasq*x[facint[ifac]]*x[facint[jfac]];
			  }
			}
			
			if ((fac_corr!=0)&&(jfac==1)) {
			  //d/dsigmadrho
			  Int_t fullHessindex = (nparam)*Getifvar(imix,ifac)+Getifvar(imix,2);
			  if (ifac==0) hessilk[fullHessindex] += ddthetasq*x[facint[0]]*dx2sc_drho[Getx2i(imix,facint[0],facint[1])];
			  if (ifac==1) {
			    hessilk[fullHessindex] += ddthetasq*x2[Getx2i(imix,facint[0],facint[1])]*dx2sc_drho[Getx2i(imix,facint[0],facint[1])];
			    hessilk[fullHessindex] += modelEval[2]*dx2sc_drho[Getx2i(imix,facint[0],facint[1])]/Getfvar(imix,1);
			    
			    //d/drho/drho
			    fullHessindex = (nparam)*Getifvar(imix,2)+Getifvar(imix,2);
			    hessilk[fullHessindex] += ddthetasq*dx2sc_drho[Getx2i(imix,facint[0],facint[1])]*dx2sc_drho[Getx2i(imix,facint[0],facint[1])];
			    hessilk[fullHessindex] += modelEval[2]*Getfvar(imix,1)*(-1.0)*x[facint[1]]/(1-Getfvar(imix,2)*Getfvar(imix,2))/sqrt(1-Getfvar(imix,2)*Getfvar(imix,2));
			  }
			}

			if (est_nmix>1) {
			  // d/dsigmadmu
			  
			  Double_t multiplier = 1.0;
			  if (imix== est_nmix-1) {
			    for (UInt_t imix_hess = 0; imix_hess < est_nmix-1 ; imix_hess++) {
			      Int_t fullHessindex = (nparam)*Getifvar(imix,ifac)+Getifmean(imix_hess,jfac);
			      if (Getifvar(imix,ifac)>Getifmean(imix_hess,jfac) ) {
				fullHessindex = (nparam)*Getifmean(imix_hess,jfac)+Getifvar(imix,ifac);
			      }
			      if ((fac_corr!=0)&&(ifac==1)) hessilk[fullHessindex] += ddthetasq*x2[Getx2i(imix,facint[0],facint[1])]*(-exp(Getfw(imix_hess)));
			      else hessilk[fullHessindex] += ddthetasq*x[facint[ifac]]*(-exp(Getfw(imix_hess)));

			    }
			  }
			  else {
			    //			  fullHessindex = (nparam)*Getifvar(imix,ifac)+Getifmean(imix,jfac);
			    Int_t fullHessindex = (nparam)*Getifvar(imix,ifac)+Getifmean(imix,jfac);
			    if (Getifvar(imix,ifac)>Getifmean(imix,jfac) ) {
			      fullHessindex = (nparam)*Getifmean(imix,jfac)+Getifvar(imix,ifac);
			    }
			    if ((fac_corr!=0)&&(ifac==1)) hessilk[fullHessindex] += ddthetasq*x2[Getx2i(imix,facint[0],facint[1])];
			    else hessilk[fullHessindex] += ddthetasq*x[facint[ifac]];
			  }
			  
			  //d/dmudrho
			  if ((fac_corr!=0)&&(ifac==1) ) {
			    if (imix== est_nmix-1) {
			      for (UInt_t imix_hess = 0; imix_hess < est_nmix-1 ; imix_hess++) {
				Int_t fullHessindex = (nparam)*Getifvar(imix,2)+Getifmean(imix_hess,jfac);
				if (Getifvar(imix,2)>Getifmean(imix_hess,jfac) ) {
				  fullHessindex = (nparam)*Getifmean(imix_hess,jfac)+Getifvar(imix,2);
				}
				//			      if (ifac==0) hessilk[fullHessindex] += ddthetasq*x[facint[0]]*(-exp(Getfw(imix_hess)));
				hessilk[fullHessindex] += ddthetasq*dx2sc_drho[Getx2i(imix,facint[0],facint[1])]*(-exp(Getfw(imix_hess)));
			      }
			    } 
			    else {
			      //			  fullHessindex = (nparam)*Getifvar(imix,ifac)+Getifmean(imix,jfac);
			      Int_t fullHessindex = (nparam)*Getifvar(imix,2)+Getifmean(imix,jfac);
			      if (Getifvar(imix,2)>Getifmean(imix,jfac) ) {
				fullHessindex = (nparam)*Getifmean(imix,jfac)+Getifvar(imix,2);
			      }
			      hessilk[fullHessindex] += ddthetasq*dx2sc_drho[Getx2i(imix,facint[0],facint[1])];
			    }
			  }
			  
			  // d/dsigmadw
			  if (imix== est_nmix-1) {
			    for (UInt_t imix_grad = 0; imix_grad < est_nmix-1 ; imix_grad++) {
			      Int_t fullHessindex = (nparam)*Getifvar(imix,ifac)+Getifw(imix_grad);
			      
			      if (Getifvar(imix,ifac)>Getifw(imix_grad) ) {
				fullHessindex = (nparam)*Getifw(imix_grad)+Getifvar(imix,ifac);
			      }
			      multiplier = (-exp(Getfw(imix_grad))*Getfmean(imix_grad,jfac));
			      if ((fac_corr!=0)&&(ifac==1)) hessilk[fullHessindex] += ddthetasq*x2[Getx2i(imix,facint[0],facint[1])]*multiplier;
			      else hessilk[fullHessindex] += ddthetasq*x[facint[ifac]]*multiplier;

			    }
			  }
			  
			  //d/dwdrho
			  if ((fac_corr!=0)&&(ifac==1) ) {
			    if (imix== est_nmix-1) {
			      for (UInt_t imix_grad = 0; imix_grad < est_nmix-1 ; imix_grad++) {
				Int_t fullHessindex = (nparam)*Getifvar(imix,2)+Getifw(imix_grad);
				
				if (Getifvar(imix,ifac)>Getifw(imix_grad) ) {
				  fullHessindex = (nparam)*Getifw(imix_grad)+Getifvar(imix,2);
				}
				multiplier = (-exp(Getfw(imix_grad))*Getfmean(imix_grad,jfac));
				hessilk[fullHessindex] += ddthetasq*dx2sc_drho[Getx2i(imix,facint[0],facint[1])]*multiplier;
			      }
			    }
			  }		



			  //			if (jfac>=ifac) {
			  // d/dmudmu
			  if (imix== est_nmix-1) {
			    for (UInt_t imixhes1= 0 ; imixhes1 < est_nmix-1 ; imixhes1++) {    //loop over each mixture
			      for (UInt_t imixhes2 = imixhes1 ; imixhes2 < est_nmix-1 ; imixhes2++) {        //loop over each mixture
				
				if (Getifmean(imixhes1,ifac) <= Getifmean(imixhes2,jfac)) {
				  Int_t fullHessindex = (nparam)*Getifmean(imixhes1,ifac)+Getifmean(imixhes2,jfac);
				  hessilk[fullHessindex] += ddthetasq*(-exp(Getfw(imixhes1)))*(-exp(Getfw(imixhes2)));
				}
			      }
			    }
			  }
			  else {
			    if (jfac>=ifac) {
			      Int_t fullHessindex = (nparam)*Getifmean(imix,ifac)+Getifmean(imix,jfac);
			      hessilk[fullHessindex] += ddthetasq;
			    }
			  }
			  //			}
			}
		      }  // 2nd loop over factors


		      //distribution parameters cross model parameters:
		      //sigmax,sigmay just sqrt(2)* x or y
		      //mux, muy just -exp(omegam) if l=lmax, otherwise no multiplier
		      //omega just -exp(omegam)mu_m if l=lmax, other zero
		      
		      for (UInt_t jparam = 0 ; jparam < nparam_models[imod] ; jparam++) {      //loop over model parameters
			
			Int_t modhessindex = (nDimModHess)*(ifac) + nfac + jparam; // gives us d/dtheta dmodel
			
			//Now we just need to add dtheta/dparameter to d/dtheta dmodel
			
			//d/dsigma dmodel
			Int_t fullHessindex = (nparam)*Getifvar(imix,ifac)+(firstpar+jparam);
			if ((fac_corr!=0)&&(ifac==1)) hessilk[fullHessindex] += modhess[modhessindex]*x2[Getx2i(imix,facint[0],facint[1])];
			else hessilk[fullHessindex] += modhess[modhessindex]*x[facint[ifac]];
			
			// d/drho dmodel
			if ((fac_corr!=0)&&(ifac==1) ) {
			  fullHessindex = (nparam)*Getifvar(imix,2)+(firstpar+jparam);
			  hessilk[fullHessindex] += modhess[modhessindex]*dx2sc_drho[Getx2i(imix,facint[0],facint[1])];
			}
			if (est_nmix>1) {
			  // d/dmu dmodel
			  modhessindex = (nDimModHess)*(ifac) + nfac + jparam;
			  
			  if (imix== est_nmix-1) {
			    for (UInt_t imix_grad = 0; imix_grad < est_nmix-1 ; imix_grad++) {
			      fullHessindex = (nparam)*Getifmean(imix_grad,ifac)+(firstpar+jparam);
			      hessilk[fullHessindex] += modhess[modhessindex]*(-exp(Getfw(imix_grad)));
			    }
			  }
			  else {
			    fullHessindex = (nparam)*Getifmean(imix,ifac)+(firstpar+jparam);
			    hessilk[fullHessindex] += modhess[modhessindex];
			  }
			  // d/domega dmodel
			  if (imix== est_nmix-1) {
			    for (UInt_t imix_grad = 0; imix_grad < est_nmix-1 ; imix_grad++) {
			      fullHessindex = (nparam)*Getifw(imix_grad)+(firstpar+jparam);
			      hessilk[fullHessindex] += modhess[modhessindex]*(-exp(Getfw(imix_grad))*Getfmean(imix_grad,ifac));
			    }
			  }
			}
		      }
		    }
		  }

		  //Model specific parameters: d/dmodel dmodel 
		  if (predicting==0) {
		    for (UInt_t iparam = 0 ; iparam < nparam_models[imod] ; iparam++) {
		      for (UInt_t jparam = iparam ; jparam < nparam_models[imod] ; jparam++) {
			Int_t fullHessindex = (nparam)*(firstpar+iparam)+(firstpar+jparam);
			Int_t modhessindex = (nDimModHess)*(nfac+iparam) + nfac + jparam;
			hessilk[fullHessindex] += modhess[modhessindex];
		      }
		    }
		  }
		} // Calc hessian
 	      } // Calc Gradient

 	    } // Check if model should be ignored
 	    firstpar += nparam_models[imod];
 	  } // loop over models
	  //	  cout << "Done with model part\n";
	  
	  // sum likelihood over integration points and mixtures
	  totalprob = totalprob + probilk;
	  
	  // sum gradiant and hessian over integration points and mixtures
	  if (stochflag>=2) {

 	    if (stochflag==3) {

	      if (predicting==0) {
		for (ifreepar = 0 ; ifreepar < nfreeparam; ifreepar++) {
		  for (UInt_t jfreepar = ifreepar ; jfreepar < nfreeparam; jfreepar++) {
		    Int_t fullHessindex = freeparlist[ifreepar]*nparam+freeparlist[jfreepar];
		    totalhess[fullHessindex] += (hessilk[fullHessindex] + gradilk[freeparlist[ifreepar]]*gradilk[freeparlist[jfreepar]])*probilk;
		  }
		}
	      }
	      else {
		for (UInt_t ifac = 0 ; ifac < nfac; ifac++) {
		  for (UInt_t jfac = ifac ; jfac < nfac; jfac++) {
		    Int_t fullHessindex = ifac*nfac+jfac;
		    totalhess[fullHessindex] += (hessilk[fullHessindex] + gradilk[ifac]*gradilk[jfac])*probilk;
		  }
		}
	      }

// 	      for (UInt_t ipar = 0 ; ipar < nparam; ipar++) {
// 		for (UInt_t jpar = ipar ; jpar < nparam; jpar++) {
// 		  Int_t fullHessindex = (nparam)*(ipar)+(jpar);
// 		  totalhess[fullHessindex] += (hessilk[fullHessindex] + gradilk[ipar]*gradilk[jpar])*probilk;
// 		}
// 	      }
 	    }
	    
	    UInt_t numberparam = nparam;
	    if (predicting==1) numberparam = nfac;
	    for (UInt_t ipar = 0 ; ipar < numberparam; ipar++) {
	      gradilk[ipar] = gradilk[ipar]*probilk;
	      totalgrad[ipar] = totalgrad[ipar] + gradilk[ipar];
	    }
	  }
	  //	} // loop over second set of int points
	  } // loop over first set of int points
      } // loop over types
    } // loop over mixture integrals
    // Calculate logLkhd and gradL for this observation
    //    if (i==1) cout << "Obs " << i <<": "<< totalprob << endl;
    Double_t repeatobs = 1.0;
    if (bootstrapobs[i]>1) repeatobs *= bootstrapobs[i];
    if (weightvar>-1) repeatobs *= data[i*nvar+weightvar];

    logLkhd += -log(totalprob)*repeatobs;
    //    printf("Finished calculating likelihoood! %16.12f for obs %d\n.", (-log(totalprob)), Int_t(data[i*nvar+indexvar]));
    //    if (i==5) assert(0);	

    //Reweight hessian for stochastic derivatives
    //    if ((stoch_deriv_frac>0.0)&&(iflag>2)) repeatobs /= stoch_deriv_frac;

    if (stochflag>=2) {
      if (predicting==0) {
	ifreepar = 0;
	for (UInt_t ipar = 0 ; ipar < nparam; ipar++) {
	  if (parfixed[ipar]==0) {
	    gradL[ifreepar] += (-totalgrad[ipar]/totalprob)*repeatobs;
	    ifreepar++;
	  }
	}
      }      else {
	// If predicting factors:
	for (UInt_t ifac = 0 ; ifac < nfac; ifac++) {
	  gradL[ifac] += (-totalgrad[ifac]/totalprob)*repeatobs;
	}
      }

      
       if (iflag==3) {
	 if (predicting==0) {
	   Int_t hessindex = 0;
	   for (ifreepar = 0 ; ifreepar < nfreeparam; ifreepar++) {
	     for (UInt_t jfreepar = ifreepar ; jfreepar < nfreeparam; jfreepar++) {
	       Int_t fullHessindex = freeparlist[ifreepar]*nparam+freeparlist[jfreepar];
	       hessL[hessindex] += (totalgrad[freeparlist[ifreepar]]*totalgrad[freeparlist[jfreepar]]/(totalprob*totalprob))*repeatobs;

	       if (stoch_deriv_frac<0.0) hessL[hessindex] += (-totalhess[fullHessindex]/totalprob)*repeatobs;
	       //	       else if (stochflag==3) hessL[hessindex] += (-totalhess[fullHessindex]/totalprob/stoch_deriv_frac)*repeatobs;
	       hessindex++;
	     }
	   }
	 }
	 else {
	   Int_t hessindex = 0;
	   for (UInt_t ifac = 0 ; ifac < nfac; ifac++) {
	     for (UInt_t jfac = ifac ; jfac < nfac; jfac++) {
	       Int_t fullHessindex = ifac*nfac+jfac;
	       hessL[hessindex] += (-totalhess[fullHessindex]/totalprob + totalgrad[ifac]*totalgrad[jfac]/(totalprob*totalprob))*repeatobs;
	       hessindex++;
	     }
	   }
	 }
// 	ifreepar = 0;
// 	for (UInt_t ipar1 = 0 ; ipar1 < nparam; ipar1++) {
// 	  for (UInt_t ipar2 = ipar1 ; ipar2 < nparam; ipar2++) {
// 	    if ((parfixed[ipar1]==0)&&(parfixed[ipar2]==0)) {
// 	      Int_t fullHessindex = (nparam)*(ipar1)+(ipar2);
// 	      hessL[ifreepar] += -totalhess[fullHessindex]/totalprob + totalgrad[ipar1]*totalgrad[ipar2]/(totalprob*totalprob);
// 	      ifreepar++;
// 	    }
// 	  }
// 	}
       }

    }

    } //check skipobs
  } // loop over observations

  delete r3;

  //  printf("Found %d skips out of %d and frac=%f; %f\n",countstochflag1,end-start, float(countstochflag1)/(float(end-start)), stoch_deriv_frac);
  
   counter++;
   return;
}

// void LkhdFcn(Int_t &npar, Double_t *gin, Double_t &f,
//  	       Double_t *par, Int_t iflag) {
//    f = 0;
//    if (!gMinuit->GetObjectFit() ||
//        gMinuit->GetObjectFit()->IsA() != TMinLkhd::Class()) return;

//    TMinLkhd *sts = (TMinLkhd*)gMinuit->GetObjectFit();
//    sts->EvalLkhd(f,gin,par,iflag,0,1);
//    if ((sts->GetCounter())%100==0) printf("Iteration %5d Found logLkhd: %f \n",sts->GetCounter(),f);
//    return;
// }

void TMinLkhd::PrintParamTab(int pmod) {
  UInt_t ipar = 0;
  if (pmod==0) {

    cout << "Models with normalizations are:\n";
    for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {
      Int_t imod = norm_models.at(ifac);
      //      cout << "fac#" << ifac << " normalized in model " << imod;
      if (imod>-1) {
	vector <Double_t> fnorm = models.at(imod).GetNorm();
	cout << "fac #" << ifac << ": " << models.at(imod).GetTitle() << " with normalization=" << fnorm.at(ifac) << "\n";
      }
    }


    cout << endl << "****Factors*************" << endl;

    //     for ( ; ipar < nfac_param ; ipar++) {
//       if (ipar < nfac) {
// 	printf("Factor%1d sigma: %8.3f +/- %8.3f\n",ipar+1, param[ipar],param_err[ipar]);
//       }
//       else {
// 	// 2-factors only:
// 	printf("Factors rho  : %8.3f +/- %8.3f\n", param[ipar],param_err[ipar]);
//       }
//     }

    for (UInt_t imix = 0 ; imix < fac_nmix ; imix++) {

      if (fac_nmix>1) printf("\n**Mixture** %1d\n",imix);
      //variances
      for (UInt_t ivar = 0; ivar < f_nvariance ; ivar++) {
	if (ivar < nfac) {
	  printf("%3d. Factor%1d sigma: %8.3f +/- %8.3f\n",ipar,ivar+1, param[ipar],param_err[ipar]);
	  ipar++;
	}
	else {
	  // 2-factors only:
	  printf("%3d. Factors rho  : %8.3f +/- %8.3f\n", ipar,param[ipar],param_err[ipar]);
	  ipar++;
	} 
     }
      if (imix<fac_nmix-1) {
	// means
	for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {
	  printf("%3d. Factor%1d mean : %8.3f +/- %8.3f\n",ipar,ifac+1, param[ipar],param_err[ipar]);
	  ipar++;
	}
	
	// weights: =0.0 means all equally weighted
	printf("%3d. Mix weight   : %8.3f +/- %8.3f\n", ipar,param[ipar],param_err[ipar]);
	ipar++;
      }
    }
    if (ipar != nfac_param) {
      cout << "ERROR ERROR: ipar!= nfac_para\n";
      assert(0);
    } 
  }
  else ipar = nfac_param;
    
//   UInt_t start = 0;
//   UInt_t end = models.size();
//   if (pmod!=0) {
//     start = pmod-1;
//     end = pmod;
//     for (Int_t imod = 0 ; imod < pmod-1 ; imod++) ipar += nparam_models[imod];
//   }
  
//   for (UInt_t imod = start ; imod < end ; imod++) {
    
//     cout << endl << "****" << models[imod].GetTitle() << " Model*************" << endl;
//     vector <Int_t> regs = models[imod].GetReg();

//     for (UInt_t ireg = 0 ; ireg < regs.size(); ireg++) {
//       printf("%3d. %13s: %8.3f +/- %8.3f\n",ipar,var_table.at(regs[ireg]).Data(),param[ipar],param_err[ipar]);
//       ipar++;
//     }

//     vector <Double_t> fnorm = models[imod].GetNorm();
//     for (UInt_t i = 0 ; i < nfac ; i++) {
//       if (fnorm.size()==0) {
// 	printf("%3d. Fac%1d loading : %8.3f +/- %8.3f\n",ipar,i+1, param[ipar],param_err[ipar]);
// 	ipar++;
//       }
//       else {
// 	if (fnorm[i]>-9998) 	printf("     Fac%1d loading : %8.3f\n",i+1,fnorm[i]);
// 	else {
// 	  printf("%3d. Fac%1d loading : %8.3f +/- %8.3f\n",ipar,i+1, param[ipar],param_err[ipar]);
// 	  ipar++;
// 	}
//       }
//     }
//     if (models[imod].GetType()==1) {
//       printf("%3d. 1/Precision  : %8.3f +/- %8.3f\n", ipar,fabs(param[ipar]),param_err[ipar]);
//       ipar++;
//     }

//   }

  // Ok, we are just going to do it the dumb way (i.e. slow/inefficient)

  printf("\nPrinting tables!!\n\n");

  //  Double_t	StudentI(Double_t T, Double_t ndf)

  //Loop through tables
  UInt_t firstmod = 0;
  UInt_t lastmod = 0;
  FILE * pFile;
  TString filename;

  for (Int_t iprgrp = 0; iprgrp <= models.back().GetPrintGroup() ; iprgrp++) {

    firstmod = lastmod;

    filename = TString("EstimationTable_").Append(models[firstmod].GetName()).Append(".tex");
    filename.Prepend(workingdir);
    pFile = fopen ( ((char *)filename.Data()),"w");

    //    printf("\n\n*****Printing table group %d\n",iprgrp);
    //    fprintf (pFile,"Hello world %s\n",models[firstmod].GetName());

    printf("\n\\begin{sidewaystable}\n");
    printf("\\begin{center}\n");
    printf("	\\small\n");
    printf("\\caption{\\label{t_est%d} Estimates for Group %d Models}\n",iprgrp,iprgrp);

    // Get first and last model of this printgroup
    UInt_t done = 0;
    while (done==0) {
      lastmod++;
      if (lastmod < models.size()) {
	if (models.at(lastmod).GetPrintGroup()!=iprgrp) done=1;
      }
      else done=1;
    }

    vector<Double_t> margeffmult(lastmod-firstmod, 1.0);
    
    //Calculate marginal effects
    if (printmargeffect==1) {
      for (UInt_t imod = firstmod ; imod < lastmod ; imod++) {
	if (models[imod].GetType()==2) {
	  Double_t count = 0.0;
	  //	  Double_t aveprob = 0.0;
	  margeffmult.at(imod-firstmod) = 0.0;
	  printf("***calculating probabilities\n");
	  for (UInt_t iobs = 0; iobs < nobs; iobs++) {
	    
	    vector<Double_t> fac_val(nfac+ntyp,0.0);
	    
	    Double_t prob = models[imod].GetPdf(iobs*nvar,data,param,fparam_models[imod],fac_val);
	    // prob is -1 if this model's indicator is zero
	    if (prob>0.0) {
	      Double_t weight = 1.0;
	      if (weightvar>-1) weight = data[iobs*nvar+weightvar];
	    
	      //	      aveprob += weight*prob;
	      //	      margeffmult.at(imod) += weight*prob*(1.0-prob);
	      margeffmult.at(imod-firstmod) += weight*prob;
	      count += weight;
	    }
	  }
	  // Get ave marginal effect
	  margeffmult.at(imod-firstmod) = margeffmult.at(imod-firstmod)/count;
	  // aveprob = aveprob/count;
	  // printf("Average prob = %f\n\n",aveprob);
	}
      }
    }

    printf("\\begin{tabular}{l");
    for (UInt_t imod = firstmod ; imod < lastmod ; imod++) {
      int nchoice = 2;
      if (models[imod].GetType()==3) nchoice = models[imod].GetNchoice();
      for (int ichoice = 2 ; ichoice <= nchoice ; ichoice++) {
	printf("|cc");
      }
    }
    printf("} \\hline \\hline\n"); 
    printf("	Variable	");
    for (UInt_t imod = firstmod ; imod < lastmod ; imod++) {
      int nchoice = 2;
      if (models[imod].GetType()==3) nchoice = models[imod].GetNchoice();
      for (int ichoice = 2 ; ichoice <= nchoice ; ichoice++) {
	printf("&	\\multicolumn{2}{c}{%19s",models[imod].GetTitle());
	if (nchoice>2) printf("_ch%1d} ",ichoice);
	else printf("} ");
      }
    }
    printf("\\\\ \\hline\n");
    printf("		    ");
    for (UInt_t imod = firstmod ; imod < lastmod ; imod++) {
      int nchoice = 2;
      if (models[imod].GetType()==3) nchoice = models[imod].GetNchoice();
      for (int ichoice = 2 ; ichoice <= nchoice ; ichoice++) {
	printf("&	$\\beta$	&	StdEr.	");
      }
    }
    printf("\\\\ \\hline\n");
	  
  // make vector of parameters used in any of the models
    vector<Int_t> varuse;
    UInt_t linmodel = 0;
    int oprob_nchoice = 0;
    // print models
    //    printf("%18s: ","Models");
    for (UInt_t imod = firstmod ; imod < lastmod ; imod++) {
      //      printf("& %19s ",models[imod].GetTitle());
    
      vector <Int_t> regs = models[imod].GetReg();
      // make list of variables used
      for (UInt_t ireg = 0 ; ireg < regs.size(); ireg++) {	
	Int_t varincluded = 0;
	for (UInt_t i = 0 ; i < varuse.size() ; i++) if (varuse.at(i)==regs.at(ireg)) varincluded=1;
	if (varincluded==0) varuse.push_back(regs.at(ireg));
      }
      if (models[imod].GetType()==1) linmodel = 1;
      if (models[imod].GetType()==4) {
	if (oprob_nchoice < models[imod].GetNchoice()) {
	  oprob_nchoice = models[imod].GetNchoice();
	}
      }
    }
    //    printf("\\\\ \n");
    //Print regressors
    for (UInt_t ivar = 0; ivar < varuse.size() ; ivar++) {
      printf("%18s ",var_table.at(varuse[ivar]).Data());
      fprintf(pFile,"%18s ",var_table.at(varuse[ivar]).Data());

      for (UInt_t imod = firstmod ; imod < lastmod ; imod++) {
	int nchoice = 1;
	if (models[imod].GetType()==3) nchoice = models[imod].GetNchoice()-1;
	for (int ichoice = 0 ; ichoice < nchoice ; ichoice++) {

	  vector <Int_t> regs = models[imod].GetReg();
	  // check to see if model used this variable
	  Int_t varused = 0;
	  for (UInt_t ireg = 0 ; ireg < regs.size(); ireg++) {	
	    if (varuse.at(ivar)==regs.at(ireg)) {
	      int iparam = fparam_models[imod] + ireg + ichoice*(models[imod].GetNreg()+nfac+ntyp*(imod!=type_model));
	      if (param_fixval.at(iparam)<-9998.0) { 
		printf("& %8.3f & %8.3f ",param[iparam]*margeffmult.at(imod-firstmod),param_err[iparam]*margeffmult.at(imod-firstmod));
		fprintf(pFile,"& %8.3f & %8.3f ",param[iparam]*margeffmult.at(imod-firstmod),param_err[iparam]*margeffmult.at(imod-firstmod));
	      }
	      else {
		printf("& %8.3f & %8.3f ",param[iparam]*margeffmult.at(imod-firstmod),0.0);
		fprintf(pFile,"& %8.3f & %8.3f ",param[iparam]*margeffmult.at(imod-firstmod),0.0);
	      }
	      varused=1;
	    }
	  }
	  if (varused==0) {
	    printf("& %8s & %8s ","","");
	    fprintf(pFile,"& %8s & %8s ","","");
	  }
	}
      }
      printf("\\\\ \n");
      fprintf(pFile,"\\\\ \n");
    }
    //Print factor loadings
    vector<Int_t> nfacparam(lastmod-firstmod,0);

    for (UInt_t ifac = 0 ; ifac < nfac + ntyp; ifac++) {
      if (ifac<nfac) {
	printf("Fac%1d loading       ",ifac);
	fprintf(pFile,"Fac%1d_loading       ",ifac);
      }
      else {
	printf("Typ%1d loading       ",ifac-nfac+1);
	fprintf(pFile,"Typ%1d_loading       ",ifac-nfac+1);
      }	
      for (UInt_t imod = firstmod ; imod < lastmod ; imod++) {
	int nchoice = 1;
	if (models[imod].GetType()==3) nchoice = models[imod].GetNchoice()-1;
	for (int ichoice = 0 ; ichoice < nchoice ; ichoice++) {

	  vector <Double_t> fnorm = models[imod].GetNorm();
	  if ( (ifac<nfac) || (imod!=type_model)) {
	    if (fnorm.size()==0) {
	      nfacparam.at(imod-firstmod) +=1;
	      Int_t iparam = fparam_models[imod] + ichoice*(models[imod].GetNreg()+nfac+ntyp*(imod!=type_model)) + models[imod].GetNreg() + ifac;
	      printf("& %8.3f & %8.3f ",param[iparam]*margeffmult.at(imod-firstmod),param_err[iparam]*margeffmult.at(imod-firstmod));
	      fprintf(pFile,"& %8.3f & %8.3f ",param[iparam]*margeffmult.at(imod-firstmod),param_err[iparam]*margeffmult.at(imod-firstmod));
	    }
	    else {
	      if (fnorm[ifac]>-9998) {
		printf("& %8.3f & %8s ",fnorm[ifac]*margeffmult.at(imod-firstmod),"");
		fprintf(pFile,"& %8.3f & %8s ",fnorm[ifac]*margeffmult.at(imod-firstmod),"");
	      }
	      else {
		nfacparam.at(imod-firstmod) +=1;
		Int_t iparam = fparam_models[imod] + ichoice*(models[imod].GetNreg()+nfac+ntyp*(imod!=type_model)) + models[imod].GetNreg() + ifac;
		printf("& %8.3f & %8.3f ",param[iparam]*margeffmult.at(imod-firstmod),param_err[iparam]*margeffmult.at(imod-firstmod));
		fprintf(pFile,"& %8.3f & %8.3f ",param[iparam]*margeffmult.at(imod-firstmod),param_err[iparam]*margeffmult.at(imod-firstmod));
	      }
	    }
	  }
	  else {
	    printf("& %8s & %8s ","","");
	    fprintf(pFile,"& %8s & %8s ","", "");
	  }
	} // loop through choices
      } // loop over models
      printf("\\\\ \n");
      fprintf(pFile,"\\\\ \n");
    }

    //Print ordered probit thresholds
    if (oprob_nchoice>0) {
      for (int ithresh = 1 ; ithresh < oprob_nchoice ; ithresh++) {
	printf("Threshold %2d       ",ithresh);
	fprintf(pFile,"Threshold %2d       ",ithresh);

	for (UInt_t imod = firstmod ; imod < lastmod ; imod++) {
	  
	  int nchoice = 1;
	  if (models[imod].GetType()==3) nchoice = models[imod].GetNchoice()-1;
	  for (int ichoice = 0 ; ichoice < nchoice ; ichoice++) {
	    if ((models[imod].GetType()==4)&& (models[imod].GetNchoice()>ithresh) ) {
	      Int_t parnum = fparam_models[imod]+models[imod].GetReg().size() + nfacparam.at(imod-firstmod)+ithresh-1;
	      printf("& %8.3f & %8.3f ",param[parnum],param_err[parnum]);
	      fprintf(pFile,"& %8.3f & %8.3f ",param[parnum],param_err[parnum]);
	    }
	    else {
	      printf("& %8s & %8s ","","");
	      fprintf(pFile,"& %8s & %8s ","","");
	    }
	  }
	}
	printf("\\\\ \n");
	fprintf(pFile,"\\\\ \n");
      }
    }

    //Print precision
    if (linmodel==1) {
      printf("%18s ","1/Precision");
      fprintf(pFile,"%18s ","1/Precision");
      
      for (UInt_t imod = firstmod ; imod < lastmod ; imod++) {
	int nchoice = 1;
	if (models[imod].GetType()==3) nchoice = models[imod].GetNchoice()-1;
	for (int ichoice = 0 ; ichoice < nchoice ; ichoice++) {
	  if (models[imod].GetType()==1) {
	    Int_t parnum = fparam_models[imod]+models[imod].GetReg().size() + nfacparam.at(imod-firstmod);
	    printf("& %8.3f & %8.3f ",param[parnum],param_err[parnum]);
	    fprintf(pFile,"& %8.3f & %8.3f ",param[parnum],param_err[parnum]);
	  }
	  else {
	    printf("& %8s & %8s ","","");
	    fprintf(pFile,"& %8s & %8s ","","");
	  }
	}
      }
      printf("\\\\ \n");
      fprintf(pFile,"\\\\ \n");
    }

      printf("%18s ","N");
      fprintf(pFile,"%18s ","N");
      
      for (UInt_t imod = firstmod ; imod < lastmod ; imod++) {
	int nchoice = 1;
	if (models[imod].GetType()==3) nchoice = models[imod].GetNchoice()-1;
	for (int ichoice = 0 ; ichoice < nchoice ; ichoice++) {
	  printf("& %8d & %8s ",nobs_models.at(imod),"");
	  fprintf(pFile,"& %8d & %8s ",nobs_models.at(imod),"");
	}
      }
      printf("\\\\ \n");
      fprintf(pFile,"\\\\ \n");

    printf("\\hline \\hline \\end{tabular}\n");
    fprintf(pFile,"\\hline \\hline \\end{tabular}\n");
    printf("\\end{center}\n");
    printf("Notes: \n");
    printf("\\end{sidewaystable}\n\n");

    fclose(pFile);
  }
}



void TMinLkhd::PrintParam(int pmod) {
  UInt_t ipar = 0;
  char constrained = ' ';

  if (pmod<1) {

    cout << "Models with normalizations are:\n";
    for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {
      Int_t imod = norm_models.at(ifac);
      //      cout << "fac#" << ifac << " normalized in model " << imod;
      if (imod>-1) {
	vector <Double_t> fnorm = models.at(imod).GetNorm();
	cout << "fac #" << ifac << ": " << models.at(imod).GetTitle() << " with normalization=" << fnorm.at(ifac) << "\n";
      }
      else cout << "fac #" << ifac << ": Variance fixed to 1.0\n";
    }

    cout << endl << "****Factors*************" << endl;

    //     for ( ; ipar < nfac_param ; ipar++) {
//       if (ipar < nfac) {
// 	printf("Factor%1d sigma: %8.3f +/- %8.3f\n",ipar+1, param[ipar],param_err[ipar]);
//       }
//       else {
// 	// 2-factors only:
// 	printf("Factors rho  : %8.3f +/- %8.3f\n", param[ipar],param_err[ipar]);
//       }
//     }

    for (UInt_t imix = 0 ; imix < fac_nmix ; imix++) {

      if (fac_nmix>1) printf("\n**Mixture** %1d\n",imix);
      //variances
      for (UInt_t ivar = 0; ivar < f_nvariance ; ivar++) {
	if (ivar < nfac) {
	  if (parconstrained.at(ipar)>-1) constrained = '*';
	  else constrained=' ';
	  printf("%3d. Factor%1d sigma: %8.3f +/- %8.3f %c \n",ipar,ivar+1, param[ipar],param_err[ipar], constrained);
	  ipar++;
	}
	else {
	  // 2-factors only:
	  if (parconstrained.at(ipar)>-1) constrained = '*';
	  else constrained=' ';
	  printf("%3d. Factors rho  : %8.3f +/- %8.3f %c \n", ipar,param[ipar],param_err[ipar], constrained);
	  ipar++;
	} 
     }
      if (imix<fac_nmix-1) {
	// means
	for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {
	  if (parconstrained.at(ipar)>-1) constrained = '*';
	  else constrained=' ';
	  printf("%3d. Factor%1d mean : %8.3f +/- %8.3f %c \n",ipar,ifac+1, param[ipar],param_err[ipar], constrained);
	  ipar++;
	}
	
	// weights: =0.0 means all equally weighted
	  if (parconstrained.at(ipar)>-1) constrained = '*';
	  else constrained=' ';
	printf("%3d. Mix weight   : %8.3f +/- %8.3f %c \n", ipar,param[ipar],param_err[ipar], constrained);
	ipar++;
      }
    }
    if (ipar != nfac_param) {
      cout << "ERROR ERROR: ipar!= nfac_para\n";
      assert(0);
    } 
  }
  else ipar = nfac_param;
    
  UInt_t start = 0;
  UInt_t end = models.size();
  if (pmod>0) {
    start = pmod-1;
    end = pmod;
    for (Int_t imod = 0 ; imod < pmod-1 ; imod++) ipar += nparam_models[imod];
  }
  
  for (UInt_t imod = start ; imod < end ; imod++) {
    int nchoice = 1;
    if (models[imod].GetType()==3) nchoice = models[imod].GetNchoice()-1;
    for (int ichoice = 0 ; ichoice < nchoice ; ichoice++) {
    
      // Print out only measurement system if pmod==0
      if ((pmod!=0)||(models[imod].GetGroup()==0)) {
	cout << endl << "****" << models[imod].GetTitle();
	if ( (models[imod].GetType()==3)&&(nchoice>1) ) cout << "(Choice " << (ichoice+2) << ") ";
	cout << " Model*************" << endl;
	vector <Int_t> regs = models[imod].GetReg();
      
	for (UInt_t ireg = 0 ; ireg < regs.size(); ireg++) {
	  if (parconstrained.at(ipar)>-1) constrained = '*';
	  else constrained=' ';
	  printf("%3d. %13s: %8.3f +/- %8.3f %c \n",ipar,var_table.at(regs[ireg]).Data(),param[ipar],param_err[ipar], constrained);
	  ipar++;
	}
	
	vector <Double_t> fnorm = models[imod].GetNorm();
	for (UInt_t i = 0 ; i < nfac+ntyp*(imod!=type_model) ; i++) {
	  if (fnorm.size()==0) {
	    if (parconstrained.at(ipar)>-1) constrained = '*';
	    else constrained=' ';
	    if (i<nfac) {
	      printf("%3d. Fac%1d loading : %8.3f +/- %8.3f %c \n",ipar,i+1, param[ipar],param_err[ipar], constrained);
	    }
	    else {
	      printf("%3d. Typ%1d loading : %8.3f +/- %8.3f %c \n",ipar,i+1-nfac, param[ipar],param_err[ipar], constrained);
	    }
	    ipar++;
	  }
	  else {
	    if (fnorm[i]>-9998) 	printf("     Fac%1d loading : %8.3f\n",i+1,fnorm[i]);
	    else {
	      if (parconstrained.at(ipar)>-1) constrained = '*';
	      else constrained=' ';
	      if (i<nfac) {
		printf("%3d. Fac%1d loading : %8.3f +/- %8.3f %c \n",ipar,i+1, param[ipar],param_err[ipar], constrained);
	      }
	      else {
		printf("%3d. Typ%1d loading : %8.3f +/- %8.3f %c \n",ipar,i+1-nfac, param[ipar],param_err[ipar], constrained);
	      }
	      ipar++;
	    }
	  }
	}
	if (models[imod].GetType()==1) {
	  if (parconstrained.at(ipar)>-1) constrained = '*';
	  else constrained=' ';
	  printf("%3d. 1/Precision  : %8.3f +/- %8.3f %c \n", ipar,fabs(param[ipar]),param_err[ipar], constrained);
	  ipar++;
	}
	else if (models[imod].GetType()==4) {
	  for (int ithresh = 1 ; ithresh < models[imod].GetNchoice() ; ithresh++) {
	    if (parconstrained.at(ipar)>-1) constrained = '*';
	    else constrained=' ';
	    if (ithresh==1) printf("%3d. Threshold %3d: %8.3f +/- %8.3f %c \n", ipar,ithresh,param[ipar],param_err[ipar], constrained);
	    else printf("%3d. Threshold %3d: %8.3f +/- %8.3f %c \n", ipar,ithresh,fabs(param[ipar]),param_err[ipar], constrained);
	    ipar++;
	  }
	}
      }
      printf("     N            : %8d \n", nobs_models.at(imod));
    }
  }
}


void TMinLkhd::PrintParam_Varlist() {


  FILE * pFile;
  TString filename;
  filename = TString("fullmodel_par_varlist.txt");
  filename.Prepend(workingdir);
  pFile = fopen ( ((char *)filename.Data()),"w");

  UInt_t ipar = 0;

  for (UInt_t imix = 0 ; imix < fac_nmix ; imix++) {

      //variances
      for (UInt_t ivar = 0; ivar < f_nvariance ; ivar++) {
	if (ivar < nfac) {
	  fprintf(pFile,"%4d factor%1d_mix%1d sigma \n",ipar,ivar,imix);

	  ipar++;
	}
	else {
	  // 2-factors only:
	  fprintf(pFile,"%4d factor_mix%1d rho \n",ipar, imix);
	  ipar++;
	} 
      }
      if (imix<fac_nmix-1) {
	// means
	for (UInt_t ifac = 0 ; ifac < nfac ; ifac++) {
	  fprintf(pFile,"%4d factor%1d_mix%1d mean \n",ipar,ifac,imix);
	  ipar++;
	}
	
	// weights: =0.0 means all equally weighted
	fprintf(pFile,"%4d factor_mix%1d weight \n",ipar,imix);
	ipar++;
      }
  }
  if (ipar != nfac_param) {
    cout << "ERROR ERROR: ipar!= nfac_para\n";
    assert(0);
  } 
    
  UInt_t start = 0;
  UInt_t end = models.size();
  
  for (UInt_t imod = start ; imod < end ; imod++) {
    int nchoice = 1;
    if (models[imod].GetType()==3) nchoice = models[imod].GetNchoice()-1;
    for (int ichoice = 0 ; ichoice < nchoice ; ichoice++) {
    
      vector <Int_t> regs = models[imod].GetReg();
      
      for (UInt_t ireg = 0 ; ireg < regs.size(); ireg++) {
	
	if  ( (models[imod].GetType()==3)&&(nchoice>2) ) {
	  fprintf(pFile,"%4d %s_ch%d %20s \n",ipar,models[imod].GetName(),ichoice+2,var_table.at(regs[ireg]).Data());
	}
	else{
	  fprintf(pFile,"%4d %20s %20s \n",ipar,models[imod].GetName(),var_table.at(regs[ireg]).Data());
	}
	ipar++;
      }
      
      vector <Double_t> fnorm = models[imod].GetNorm();
      for (UInt_t i = 0 ; i < nfac+ntyp*(imod!=type_model) ; i++) {
	if (fnorm.size()==0) {
	  
	  if  ( (models[imod].GetType()==3)&&(nchoice>2) ) {
	    if (i<nfac) fprintf(pFile,"%4d %s_ch%d factor%d \n",ipar,models[imod].GetName(),ichoice+2,i+1);
	    else fprintf(pFile,"%4d %s_ch%d type%d \n",ipar,models[imod].GetName(),ichoice+2,i+1-nfac);
	  }
	  else{
	    if (i<nfac) fprintf(pFile,"%4d %20s factor%d \n",ipar,models[imod].GetName(),i+1);
	    else fprintf(pFile,"%4d %20s type%d \n",ipar,models[imod].GetName(),i+1-nfac);
	  }
	  
	  ipar++;
	}
	else {
	  if (fnorm[i]<-9998) {
	    if  ( (models[imod].GetType()==3)&&(nchoice>2) ) {
	      if (i<nfac) fprintf(pFile,"%4d %s_ch%d factor%d \n",ipar,models[imod].GetName(),ichoice+2,i+1);
	      else fprintf(pFile,"%4d %s_ch%d type%d \n",ipar,models[imod].GetName(),ichoice+2,i+1-nfac);
	    }
	    else{
	      if (i<nfac) fprintf(pFile,"%4d %20s factor%d \n",ipar,models[imod].GetName(),i+1);
	      else fprintf(pFile,"%4d %20s type%d \n",ipar,models[imod].GetName(),i+1-nfac);
	    }
	    ipar++;
	  }
	}
      }
      if (models[imod].GetType()==1) {
	fprintf(pFile,"%4d %20s sigma \n",ipar,models[imod].GetName());
	ipar++;
      }
      else if (models[imod].GetType()==4) {
	for (int ithresh = 1 ; ithresh < models[imod].GetNchoice() ; ithresh++) {
	  fprintf(pFile,"%4d %20s threshold%d \n",ipar,models[imod].GetName(),ithresh);
	  ipar++;
	}
      }
    }
  }
  
  fclose(pFile);

}


int fexist(const char *filename ) {
  int status = 0;
  ifstream in1;
  in1.open(filename);
  if (in1.good()) status = 1;
  in1.close();
  return status;
}

int dexist(const char *dirname ) {
  //  string strPath = dirname;
  //  if ( access( strPath.c_str(), 0 ) == 0 )

//   if ( access( dirname, 0 ) == 0 )
//     {
  struct stat status;
  //      stat( strPath.c_str(), &status );
  stat( dirname, &status );
  
  if ( status.st_mode & S_IFDIR )
    {
      //	  cout << "The directory exists." << endl;
      return 0;
    }
  else
    {
      //	  cout << "Path doesn't exist." << endl;
      return 1;
    }
}

