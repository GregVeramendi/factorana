
#include "TModel.hh"

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
#include "TMatrix.h"
#include "TRandom.h"
#include "TMinuit.h"
#include "TPaveText.h"
#include "TFile.h"
#include "Math/DistFunc.h"

#include <cmath>
#include <iostream>
#include <assert.h>
#include <random>

using namespace std;

static std::mt19937_64 mt(7);

TModel::TModel(const char *name, const char *title, Int_t modeltype, Int_t modelgroup, Int_t prntgroup,
	       std::vector<Int_t> & moddata, Int_t nfac, Double_t * thisnormfac,
	       UInt_t nchoice, UInt_t nrank)
  : TNamed(name,title)
{
  ignore = 0;
  splitsim = -1;
  detailsim = 0;
  modtype = modeltype;
  modgroup = modelgroup;
  printgroup = prntgroup;
  numchoice = nchoice;
  numrank = nrank;
  rankshare = -9999.0;
  outcome = moddata[0];
  missing = moddata[1];
  nregressors = moddata.size()-2;
  regressors.reserve(nregressors);
  for (int i = 0 ; i < nregressors ; i++) regressors.push_back(moddata[i+2]);

  simresult = -9999.0;
  //This list records the outcome variable names
  endogVarList.clear();

  //These lists will be used for simulation
  endogRegList.clear();
  endogModelList.clear();
  endogChoiceList.clear();

//  printf("Regressors: ");
//  for (int i = 0 ; i < nregressors ; i++) {
//    if (i!=0) cout << ", "; 
//    cout << regressors[i];
//    //    cout << vartab[regressors[i]].Data();
//  }
//  printf("\n");

  if (outcome==-1 && modtype>1) {
      cout << "ERROR (TModel::TModel): Outcome \"none\" found for discrete choice model!"
	   << endl;
      assert(0);
  }

  numfac=nfac;
  if (thisnormfac) {
    facnorm.reserve(nfac);
    for (int i = 0 ; i < nfac ; i++) {
      facnorm.push_back(thisnormfac[i]);
    }
  }
}

TModel::~TModel() {
  regressors.clear();
  facnorm.clear();

}

void TModel::SplitSim(UInt_t ivar) {
  splitsim = (Int_t)ivar;
}


void SetRankShareVar(UInt_t rankshare) {
  if ((modtype==3) && (numrank>1) && (rankshare >= 0)) {
    ranksharevar = rankshare;
  }
  else {
    cout << "ERROR (TModel::SetRankShareVar): Either incorrect model type or negative varnumber!"
	 << "\n\t Model type " << modetype << " nrank=" << numrank 
	 << endl;
    assert(0);
  }

}


//void TModel::SetEndogenousReg(UInt_t endModel, std::vector<TModel> & models, std::vector<TString> & vartab) {
void TModel::SetEndogenousRegs(std::vector<TModel> & models, std::vector<TString> & vartab) {


  //These lists will be used for simulation
  endogRegList.clear();
  endogModelList.clear();
  endogChoiceList.clear();
  
  
  //Loop over list of endogenous variables
  for (UInt_t iendogvarnum = 0; iendogvarnum < endogVarList.size() ; iendogvarnum++) {

    UInt_t endModel = 9999;
    //Find model of endogenous variable:
    for (UInt_t imod = 0 ; imod < models.size() ; imod++) {
      if (models.at(imod).GetOutcome()==endogVarList.at(iendogvarnum)) {
	endModel = imod;
	break;
      }
    }
    if (endModel==9999) {
      cout << "ERROR (TModel::SetEndogenousReg): Endogenous regressor is not an outcome in any model!"
	   << "\n\t We were looking for var " << vartab.at(endogVarList.at(iendogvarnum))
	   << endl;
      assert(0);
    }
    
    Int_t foundreg = 0;

    TString outcomename;
    if (models.at(endModel).GetOutcome() > -1) outcomename = vartab.at(models.at(endModel).GetOutcome());
    else outcomename = "none";

    
    for (int ireg = 0 ; ireg < nregressors ; ireg++) {

      // Look for endogenous regressors from linear or probit models
      if (( models.at(endModel).GetType()==1 && models.at(endModel).ModelHasXtiles()==0 )
	  || (models.at(endModel).GetType()==2) ) {
	if (models.at(endModel).GetOutcome()==regressors[ireg]) {
	  foundreg++;
	  endogRegList.push_back(models.at(endModel).GetOutcome());
	  endogModelList.push_back(endModel); 
	  endogChoiceList.push_back(-1);
	  
	  cout << "Model "<< models.at(endModel).GetName() 
	       << " being set as endogenous covariate (" 
	       << models.at(endModel).GetOutcome() << ") in model " 
	       << GetName() << endl;
	}
      }
      else if ((models.at(endModel).GetType()==3)||(models.at(endModel).GetType()==4)||
	       ( models.at(endModel).GetType()==1 && models.at(endModel).ModelHasXtiles()==1 ) ) {
	
	//Look for endogenous regressors from multinomial models
	

	if (vartab.at(regressors[ireg]).BeginsWith(outcomename)) {
	  
	  TString thisendogreg = vartab.at(regressors[ireg]);
	  
	  int choicelength = thisendogreg.Length() - outcomename.Length();
	  if ((choicelength<1)||(choicelength>2)) {
	    cout << "Found regressor " << thisendogreg;
	    cout << "TModel::SetEndogenousVar(): choicelength can only be 1-2 characters!!" << endl;
	    assert(0);
	  }
	  
	  char charchoice = thisendogreg[thisendogreg.Length()-1];
	  int choice = -1;
	  //Use the fact that in ascii 48 is '0', 49 is '1', etc.
	  if ((charchoice>47)&&(charchoice<58)) {
	    choice = charchoice-48;
	    
	    if (choicelength==2) {
	      charchoice = thisendogreg[thisendogreg.Length()-2];
	      if ((charchoice>47)&&(charchoice<58)) {
		choice += 10*(charchoice-48);
	      }
	      else choice=-1;
	    }
	  }
	  if ((choice>0)&&(choice<100)) {
	    foundreg++;
	    endogRegList.push_back(regressors[ireg]);
	    endogModelList.push_back(endModel);
	    endogChoiceList.push_back(choice);

	  cout << "Model "<< models.at(endModel).GetName() 
	       << " being set as endogenous covariate (" 
	       << models.at(endModel).GetOutcome() << ") in model " 
	       << GetName() << ", where choice=" << choice << endl;
	  }
	  else {
	    cout << endl << "ERROR (TModel::SetEndogenousReg): Last character has to be a number between 1-99!!!"
		 << "\n\t We were looking for var " << outcomename
	      //		 << "\n\t We were looking for var " << vartab.at(models.at(endModel).GetOutcome())
		 << " and found " << thisendogreg
		 << " choicelength=" << choicelength
		 << " in model " << GetName() << endl;
	    
	    assert(0);
	  }
	}
      }
      

      //Look for missing indicator
      //We will assume that there is one indicator for multinomial choice models
      //TString varname_miss = vartab.at(models.at(endModel).GetOutcome()) + "_miss";
      TString varname_miss = outcomename + "_miss";
      
      if (vartab.at(regressors[ireg])==varname_miss) {
      	foundreg++;
	cout << "TModel::SetEndogenousVar(): Found missing indicator " << varname_miss << endl;
	endogRegList.push_back(regressors[ireg]);
	endogModelList.push_back(-1); 
	endogChoiceList.push_back(-1);
      }
      
    }
    
  
    if (foundreg==0) {
      cout << "ERROR (TModel::SetEndogenousReg): Endogenous regressor is not in list of covariates!"
	//	   << "\n\t We were looking for var " << vartab.at(models.at(endModel).GetOutcome())
	   << "\n\t We were looking for var " << outcomename
	   << " in model " << GetName()
	   << endl;
      assert(0);
    }
  }
}

//void TModel::SetEndogenousMajor(UInt_t endModel, std::vector<TModel> & models, std::vector<TString> & vartab) {
//
//  Int_t foundreg = 0;
//  TString outcomename = vartab.at(models.at(endModel).GetOutcome());
//  cout << "Looking for major regressors that start with " << outcomename << endl;
//
//  for (int ireg = 0 ; ireg < nregressors ; ireg++) {
//    if (vartab.at(regressors[ireg]).BeginsWith(outcomename)) {
//
//      TString thisendogreg = vartab.at(regressors[ireg]);
//
//      cout << "Found regressor " << thisendogreg;
//
//      char cmaj = thisendogreg[thisendogreg.Length()-1];
//      int maj = -1;
//      //Use the fact that in ascii 48 is '0', 49 is '1', etc.
//      if ((cmaj>48)&&(cmaj<55)) {
//	maj = cmaj-48;
//	foundreg++;
//	cout << ", where major=" << maj << endl;
//      }
//      else {
//	cout << "TModel::SendEndogenousMajor(): Last character has to be a number between 1-6!!" << endl;
//	assert(0);
//      }
//      endogRegList.push_back(regressors[ireg]);
//      endogModelList.push_back(endModel);
//      endogChoiceList.push_back(maj);
//    }
//  }
//  cout << "Found " << foundreg << " endogenous major indicators!\n";
//
//  if (foundreg==0) {
//    cout << "ERROR (TModel::SetEndogenousReg): Endogenous regressor is not in list of covariates!"
//         << "\n\t We were looking for var#" << models.at(endModel).GetOutcome()
//	 << " in model " << models.back().GetName()
//	 << endl;
//    assert(0);
//  }
//}


void TModel::PrintModel(std::vector<TString> vartab)
//void TModel::Print(TString * vartab)
{

  cout << "****";
  if (modtype==1) cout << "Linear";
  else if (modtype==2) cout << "Probit";
  else if (modtype==3) cout << "Logit (nchoice=" << numchoice << ") ";
  else if (modtype==4) cout << "Ordered Probit (nchoice=" << numchoice << ") ";

  cout << " model for " << this->GetTitle() << "\n";
  printf(" Model groups (estimation, printing): %4d, %4d \n", modgroup,printgroup);

  if (outcome==-1) printf(" Outcome:         none\n");
  else printf(" Outcome: %12s\n",vartab.at(outcome).Data());
  if (missing==-1) printf(" Missing:         none\n");
  else printf(" Missing: %12s\n",vartab.at(missing).Data());

  printf(" Regressors: ");
  for (int i = 0 ; i < nregressors ; i++) {
    if (i!=0) cout << ", "; 
    //    cout << regressors[i];
    cout << vartab.at(regressors[i]).Data();
  }
  printf("\n");

//  printf("Reg pointer: ");
//  cout << regressors << endl;
  
  printf(" Normalizations for factors: ");
  if (facnorm.size()) {
    for (uint i = 0 ; i < facnorm.size() ; i++) {
      if (i!=0) cout << ", " ;
      if (facnorm[i]>-9998) cout << facnorm[i];
       else cout << "none"; 
    }
  }
  else cout << "none";
  
  printf("\n\n");
}


void TModel::Eval(UInt_t iobs_offset,const std::vector<Double_t> & data,const std::vector<Double_t> & param, UInt_t firstpar, const std::vector<Double_t> & fac, std::vector<Double_t> & modEval,std::vector<Double_t> & hess, Int_t flag)
{
  //  cout << "Entering model " << this->GetTitle() << ", flag=" << flag <<"\n";

  //  std::vector<Double_t> modEval;
  if (flag >= 2) {

    //General parameters in gradient
    // 1 (likelihood) 2*numfac (df/dtheta and df/dalpha) df/dbeta
    Int_t ngrad = 1+2*numfac+nregressors;

    //Model specific parameters:
    if (modtype==1) ngrad += 1; // variance of error term
    if ((modtype==3)&&(numchoice>2)) ngrad = 1 + numfac + (numchoice-1)*(numfac+nregressors); // multinomial logit
    if (modtype==4) ngrad += (numchoice-1); // ordered probit intercepts
    modEval.resize(ngrad);
    for (int i = 1 ; i < ngrad ; i++) modEval[i] = 0.0;
    if (flag==3) hess.clear();
  }
  else modEval.resize(1);

  modEval[0]=1.0;

  if (missing>-1) {
    if (data[iobs_offset+missing]==0) return;
  }
  if (ignore) return;

  Int_t numlogitchoice = 2;
  vector<double> expres(1,0.0);

  if ((modtype==3)&&(numchoice>2)) {
    expres.resize(numchoice-1,0.0);
    numlogitchoice = numchoice;
  }

  UInt_t ifreefac = 0;

  for (int ichoice = 0; ichoice < numlogitchoice-1; ichoice++) {
    ifreefac = 0;
    UInt_t nparamchoice = nregressors + numfac;

    for (int i = 0; i < nregressors  ; i++) {
      expres[ichoice] += param[i+firstpar+ichoice*nparamchoice]*data[iobs_offset+regressors[i]];
    }
    
    for (int i = 0 ; i < numfac ; i++) {
      if (facnorm.size()==0) {
	expres[ichoice] += param[ifreefac+firstpar+ichoice*nparamchoice+nregressors]*fac[i];
	ifreefac++;
      }
      else {
	if (facnorm[i]>-9998) expres[ichoice] += facnorm[i]*fac[i];
	else {
	  expres[ichoice] += param[ifreefac+firstpar+ichoice*nparamchoice+nregressors]*fac[i];
	  ifreefac++;
	}
      }
    }
  } // loop over choices

  // nparameters: d/dtheta, d/dbeta, d/dalpha
  Int_t npar = numfac+nregressors+ifreefac;
  if ((modtype==3)&&(numchoice>2)) npar = numfac + (numchoice-1)*(nregressors+numfac);


  if (modtype==1) {
    Double_t Z = 0.0;
    if (outcome>-1) Z = data[iobs_offset+outcome]-expres[0];
    Double_t sigma = fabs(param[firstpar+nregressors+ifreefac]);

    //Now find the density:
    modEval[0] =ROOT::Math::normal_pdf(Z,sigma);
    //    return ROOT::Math::normal_pdf(data[iobs_offset+outcome]-Z[0],fabs(param[firstpar+nregressors+ifreefac]));
    //Now find the derivatives:
    if (flag>=2) {
      //      Int_t npar = numfac+ifreefac+nregressors+1;
      npar++; // add 1 for sigma
      if (flag==3) {
	hess.resize(npar*npar);
	for (int i = 0; i < npar ; i++) {
	  for (int j = i ; j < npar ; j++) {
	    hess[i*npar+j] = -1.0/(sigma*sigma);
	  }
	}
      }

      ifreefac = 0;
      for (int i = 0 ; i < numfac ; i++) {
	if (facnorm.size()==0) {
	  // gradient of factor-specific parameters (variance, mean, weights) d/dtheta
	  //	  modEval[i+1] = Z*(fac[i]*param[ifreefac+firstpar+nregressors]/param[i])/(sigma*sigma);
	  modEval[i+1] = Z*(param[ifreefac+firstpar+nregressors])/(sigma*sigma);
	  //gradient of factor loading (d/dalpha)
	  modEval[1+numfac+nregressors+ifreefac] = Z*fac[i]/(sigma*sigma);

	  if (flag==3) {
	    // d/dtheta_i row
	    for (int j = i; j < npar ; j++) hess[i*npar+j] *= param[ifreefac+firstpar+nregressors];
	    // d/dtheta_i col
	    for (int j = 0; j <= i ; j++) hess[j*npar+i] *= param[ifreefac+firstpar+nregressors];

	    // alpha_i index
	    Int_t index = numfac+nregressors+ifreefac;
	    // d/dalpha_i row
	    for (int j = index; j < npar ; j++) hess[index*npar+j] *= fac[i];
	    // d/dalpha_i col
	    for (int j = 0; j <= index ; j++) hess[j*npar+index] *= fac[i];

	    hess[i*npar+index] += Z/(sigma*sigma);
	  }
	  ifreefac++;
	}
	else {
	  if (facnorm[i]>-9998.0) {
	    // gradient of factor-specific parameters (variance, mean, weights)
	    //	    modEval[i+1] = Z*(fac[i]*facnorm[i]/param[i])/(sigma*sigma);
	    modEval[i+1] = Z*facnorm[i]/(sigma*sigma);

	    if (flag==3) {
	      // d/dtheta_i row
	      for (int j = i; j < npar ; j++) hess[i*npar+j] *= facnorm[i];
	      // d/dtheta_i col
	      for (int j = 0; j <= i ; j++) hess[j*npar+i] *= facnorm[i];
	    }
	  }
	  else {
	    // gradient of factor-specific parameters (variance, mean, weights)
	    //	    modEval[i+1] = Z*(fac[i]*param[ifreefac+firstpar+nregressors]/param[i])/(sigma*sigma);
	    modEval[i+1] = Z*param[ifreefac+firstpar+nregressors]/(sigma*sigma);
	    //gradient of factor loading (alpha)
	    modEval[1+numfac+nregressors+ifreefac] = Z*fac[i]/(sigma*sigma);

	    if (flag==3) {
	      // d/dtheta_i row
	      for (int j = i; j < npar ; j++) hess[i*npar+j] *= param[ifreefac+firstpar+nregressors];
	      // d/dtheta_i col
	      for (int j = 0; j <= i ; j++) hess[j*npar+i] *= param[ifreefac+firstpar+nregressors];
	      
	      // alpha_i index
	      Int_t index = numfac+nregressors+ifreefac;
	      // d/dalpha_i row
	      for (int j = index; j < npar ; j++) hess[index*npar+j] *= fac[i];
	      // d/dalpha_i col
	      for (int j = 0; j <= index ; j++) hess[j*npar+index] *= fac[i];
	      
	      hess[i*npar+index] += Z/(sigma*sigma);
	    }

	    ifreefac++;
	  }
	}
      }
      //Gradients of betas:
      for (int i = 0; i < nregressors  ; i++) {
	modEval[1+numfac+i] = Z*data[iobs_offset+regressors[i]]/(sigma*sigma);

	if (flag==3) {
	  // alpha_i index
	  Int_t index = numfac+i;
	  // d/dalpha_i row
	  for (int j = index; j < npar ; j++) hess[index*npar+j] *= data[iobs_offset+regressors[i]];
	  // d/dalpha_i col
	  for (int j = 0; j <= index ; j++) hess[j*npar+index] *= data[iobs_offset+regressors[i]];
	}
      }

      // gradient of precision:
      modEval[1+numfac+nregressors+ifreefac] = (Z*Z/(sigma*sigma*sigma)) - (1.0/sigma);
      if (flag==3) {
	// sigma_i index
	Int_t index = numfac+nregressors+ifreefac;
	// d/dsigma_i row
	for (int j = index; j < npar ; j++) hess[index*npar+j] *= 2.0*Z/sigma;
	// d/dsigma_i col
	for (int j = 0; j < index ; j++) hess[j*npar+index] *= 2.0*Z/sigma;

	// separately handle diagonal element
	hess[index*npar+index] = -3.0*Z*Z/pow(sigma,4) + 1.0/(sigma*sigma);
      }


    }
    //    cout << "Returning hessian size=" << hess.size() << "\n";


    return;
    //    return ROOT::Math::normal_pdf(data[iobs_offset+outcome]-Z[0],fabs(param[firstpar+nregressors+ifreefac]));
  }
  // Probit model
  else if (modtype==2) {
//     if (( int(data[iobs_offset+outcome])!=1)&&( int(data[iobs_offset+outcome])!=0)){
//       cout << "ERROR (TModel::Eval): Found non-binary number for probit outcome!"
// 	   << " In model " <<  this->GetTitle() << "\n"
// 	   << "Looking at outcome #" << outcome 
// 	   << ", Missing (" << missing << ")=" << data[iobs_offset+missing]
// 	   << ", obs#:" << iobs_offset/379.0
// 	   << ", Value=" << data[iobs_offset+outcome] << endl;
//       assert(0);
//     }

//     if ( int(data[iobs_offset+outcome])==1) modEval[0] = ROOT::Math::normal_cdf(Z[0]);
//     else if ( int(data[iobs_offset+outcome])==0) modEval[0] = ROOT::Math::normal_cdf(-Z[0]);

    // Was state observed?
    Double_t obsSign = 1.0;
    if ( int(data[iobs_offset+outcome])==0) obsSign = -1.0;

    //Now find the density:
    modEval[0] = ROOT::Math::normal_cdf(obsSign*expres[0]);
//     if ((obsSign*Z[0]<-35.0)&&(flag!=2)) {
//       modEval[0] = 1.0e-50/fabs(Z[0]);
//       //      cout << "TModel: Z[0]s too small, artificially fixing it\n";
//     }
    
    // Now find the derivatives: 
    if (flag>=2) {
      Double_t Z = obsSign*expres[0];
      Double_t pdf = ROOT::Math::normal_pdf(obsSign*expres[0]);   //(obsSign here is not necessary)
      //      if (fabs(expres[0])>35.0) pdf = 1.0e-30/fabs(expres[0]);
      Double_t cdf = modEval[0];
      if (obsSign*expres[0]<-35.0) cdf = 1.0e-50;

      if (flag==3) {
	hess.resize(npar*npar,0.0);
	for (int i = 0; i < npar ; i++) {
	  for (int j = i ; j < npar ; j++) {
	    hess[i*npar+j] = 1.0;
	  }
	}
      }

      ifreefac = 0;
      for (int i = 0 ; i < numfac ; i++) {
	if (facnorm.size()==0) {
	  // gradient of factor-specific parameters (variance, mean, weights)
	  //	  modEval[i+1] = pdf*(obsSign*fac[i]*param[ifreefac+firstpar+nregressors]/param[i])/cdf; 
	  modEval[i+1] = pdf*obsSign*param[ifreefac+firstpar+nregressors]/cdf; 
	  //gradient of factor loading (alpha)
	  modEval[1+numfac+nregressors+ifreefac] =  pdf*(obsSign*fac[i])/cdf;
// 	  if (obsSign*expres[0]<-35.0) modEval[1+numfac+nregressors+ifreefac] = -0.1*param[ifreefac+firstpar+nregressors]*param[ifreefac+firstpar+nregressors]/fabs(param[ifreefac+firstpar+nregressors]);

	  if (flag==3) {
	    // lambda^L(theta) - lambda^Prob(theta) (row)  
	    for (int j = i; j < npar ; j++) hess[i*npar+j] *= -Z*obsSign*param[ifreefac+firstpar+nregressors] - modEval[i+1];
	    // dZ/dtheta_i (col)
	    for (int j = 0; j <= i ; j++) hess[j*npar+i] *= obsSign*param[ifreefac+firstpar+nregressors];

	    // alpha_i index
	    Int_t index = numfac+nregressors+ifreefac;
	    // lambda^L(alpha) - lambda^Prob(alpha) (row)
	    for (int j = index; j < npar ; j++) hess[index*npar+j] *= -Z*obsSign*fac[i] - modEval[1+numfac+nregressors+ifreefac];
	    // dZ/dalpha (col)
	    for (int j = 0; j <= index ; j++) hess[j*npar+index] *= obsSign*fac[i];
	  }

	  ifreefac++;
	}
	else {
	  if (facnorm[i]>-9998.0) {
	    // gradient of factor-specific parameters (variance, mean, weights)
	    //	    modEval[i+1] = pdf*(obsSign*fac[i]*facnorm[i]/param[i])/cdf;
	    modEval[i+1] = pdf*obsSign*facnorm[i]/cdf;

	    if (flag==3) {
	      // lambda^L(theta) - lambda^Prob(theta) (row)  
	      for (int j = i; j < npar ; j++) hess[i*npar+j] *= -Z*obsSign*facnorm[i] - modEval[i+1];
	      // dZ/dtheta_i (col)
	      for (int j = 0; j <= i ; j++) hess[j*npar+i] *= obsSign*facnorm[i];
	    }
	  }
	  else {
	  // gradient of factor-specific parameters (variance, mean, weights)
	    //	    modEval[i+1] = pdf*(obsSign*fac[i]*param[ifreefac+firstpar+nregressors]/param[i])/cdf; 
	    modEval[i+1] = pdf*obsSign*param[ifreefac+firstpar+nregressors]/cdf; 
	  //gradient of factor loading (alpha)
	    modEval[1+numfac+nregressors+ifreefac] = pdf*(obsSign*fac[i])/cdf;
// 	    if (obsSign*expres[0]<-35.0) modEval[1+numfac+nregressors+ifreefac] = -0.1*param[ifreefac+firstpar+nregressors]*param[ifreefac+firstpar+nregressors]/fabs(param[ifreefac+firstpar+nregressors]);

	    if (flag==3) {
	      // lambda^L(theta) - lambda^Prob(theta) (row)  
	      for (int j = i; j < npar ; j++) hess[i*npar+j] *= -Z*obsSign*param[ifreefac+firstpar+nregressors] - modEval[i+1];
	      // dZ/dtheta_i (col)
	      for (int j = 0; j <= i ; j++) hess[j*npar+i] *= obsSign*param[ifreefac+firstpar+nregressors];
	      
	      // alpha_i index
	      Int_t index = numfac+nregressors+ifreefac;
	      // lambda^L(alpha) - lambda^Prob(alpha) (row)
	      for (int j = index; j < npar ; j++) hess[index*npar+j] *= -Z*obsSign*fac[i] - modEval[1+numfac+nregressors+ifreefac];
	      // dZ/dalpha (col)
	      for (int j = 0; j <= index ; j++) hess[j*npar+index] *= obsSign*fac[i];
	    }
	    
	    ifreefac++;
	  }
	}
      }

      for (int ireg = 0; ireg < nregressors  ; ireg++) {
	modEval[ireg+numfac+1] = pdf*(obsSign*data[iobs_offset+regressors[ireg]])/cdf;
	if (flag==3) {
	  Int_t index = numfac+ireg;
	  // lambda^L(theta) - lambda^Prob(theta) (row)  
	  for (int j = index; j < npar ; j++) hess[index*npar+j] *= -Z*obsSign*data[iobs_offset+regressors[ireg]] - modEval[ireg+numfac+1];
	  // dZ/dtheta_i (col)
	  for (int j = 0; j <= index ; j++) hess[j*npar+index] *= obsSign*data[iobs_offset+regressors[ireg]];
	} 
      }

      if (flag==3) {
	//Need to add dZ/dtheta dalpha
	ifreefac = 0;
	for (int i = 0 ; i < numfac ; i++) {
	  if (facnorm.size()==0) {
	    Int_t index = numfac+nregressors+ifreefac;
	    hess[i*npar+index] += obsSign;
	    ifreefac++;
	  }
	  else if (facnorm[i]<-9998.0) {
	    Int_t index = numfac+nregressors+ifreefac;
	    hess[i*npar+index] += obsSign;
	  ifreefac++;
	  }
	}

	// Finally multiply by pdf/cdf
	for (int i = 0; i < npar ; i++) {
	  for (int j = i ; j < npar ; j++) {
	    hess[i*npar+j] *= pdf/cdf;
	  }
	}
      }

    } 
    //    cout << "Returning hessian size=" << hess.size() << "\n";
    // check for NaN
//     for (UInt_t i = 0; i <  modEval.size() ; i++) 
//       if ( std::isnan( modEval[i])) cout << "Found the NaN!! model type 2, element " << i << " modEval[0] = " << modEval[0] << " expess=" << expres[0] << "\n";

    return;
    
    //     if ( int(data[iobs_offset+outcome])==1) return ROOT::Math::normal_cdf(expres[0]);
    //     else if ( int(data[iobs_offset+outcome])==0) return ROOT::Math::normal_cdf(-expres[0]);
  }

  // Ordered probit model:
  else if (modtype==4) {

    // Was state observed?
    Int_t obsCat = -1;
    for (int icat = 1 ; icat <= numchoice ; icat++) if (icat == int(data[iobs_offset+outcome])) obsCat = icat;
    
    if (obsCat==-1) {
      cout << "ERROR (TModel::Eval): Found invalid number for ordered probit outcome!"
 	   << " In model " <<  this->GetTitle() << "\n"
 	   << "Looking at outcome #" << outcome 
 	   << ", Missing (" << missing << ")=" << data[iobs_offset+missing]
 	   << ", obs#:" << iobs_offset/379.0
 	   << ", Value=" << data[iobs_offset+outcome] << endl;
      
      assert(0);
    }

    // if obsCat==1, then threshold[0] is ignored
    // if obsCat==numchoice, then threshold[1] is ignored
    double threshold[2];
    threshold[0] = param[firstpar+nregressors+ifreefac];
    threshold[1] = param[firstpar+nregressors+ifreefac];
    for (int icat = 2 ; icat <= obsCat ; icat++) {
      if (icat<obsCat) threshold[0] += abs(param[firstpar+nregressors+ifreefac+icat-1]);
      if (icat<numchoice) threshold[1] += abs(param[firstpar+nregressors+ifreefac+icat-1]);
    }


    double CDF[2] = {0.0, 1.0};
    if (obsCat>1) {
      CDF[0] = ROOT::Math::normal_cdf(threshold[0] - expres[0]);
    }
    if (obsCat<numchoice) {
      CDF[1] = ROOT::Math::normal_cdf(threshold[1] - expres[0]);
    }

    //Now find the density:
    modEval[0] = CDF[1] - CDF[0]; 

    //    printf("icat=%d, threshold0=%f, threshold1=%f, lkhd=%f\n",obsCat,threshold[0],threshold[1],modEval[0]);

    // Fix numerical problems (will divide by CDF later)
    if ((obsCat>1) && (CDF[0]<1.0e-50)) CDF[0] = 1.0e-50;
    if (CDF[1]<1.0e-50) CDF[1] = 1.0e-50;
    double diffCDF = modEval[0];
    if (diffCDF<1.0e-50) diffCDF = 1.0e-50;

    // Now find the derivatives: 
    if (flag>=2) {
      Double_t Z[2] = {-9999.0, -9999.0};
      Double_t PDF[2] = {0.0, 0.0};

      if (obsCat>1) {
	Z[0] = threshold[0] - expres[0];
	PDF[0] = ROOT::Math::normal_pdf(threshold[0] - expres[0]);
      }
      if (obsCat<numchoice) {
	Z[1] = threshold[1] - expres[0];
	PDF[1] = ROOT::Math::normal_pdf(threshold[1] - expres[0]);
      }


      npar += numchoice - 1; // add intercepts

      if (flag==3) {
	hess.resize(npar*npar,0.0);
	for (int i = 0; i < npar ; i++) {
	  for (int j = i ; j < npar ; j++) {
	    hess[i*npar+j] = 0.0;
	  }
	}
      }

      // ****************************
      // First calculate the gradients
      // ****************************

      // loop over two terms
      for (int iterm = 0; iterm < 2 ; iterm++) {
	// skip if they are end categories
	if ( ((obsCat>1)&&(iterm==0)) || ((obsCat<numchoice)&&(iterm==1)) ) {
	  
	  std::vector<double> tmpgrad(npar, 0.0);
	  
	  // First do the factor-specific gradient and hessian terms (theta and alpha)
	  ifreefac = 0;
	  for (int ifac = 0 ; ifac < numfac ; ifac++) {
	    //No normalizations:
	    if (facnorm.size()==0) {
	      // gradient of factor-specific parameters (variance, mean, weights)
	      tmpgrad[ifac] = (-1.0*param[ifreefac+firstpar+nregressors])*PDF[iterm]/diffCDF; 
	      //gradient of factor loading (alpha)
	      tmpgrad[numfac+nregressors+ifreefac] = (-1.0*fac[ifac])*PDF[iterm]/diffCDF;
	      
	      ifreefac++;
	    }
	    else {
	      if (facnorm[ifac]>-9998.0) {
		// gradient of factor-specific parameters (variance, mean, weights)
		tmpgrad[ifac] = (-1.0*facnorm[ifac])*PDF[iterm]/diffCDF;
	      }
	      else {
		// gradient of factor-specific parameters (variance, mean, weights)
		tmpgrad[ifac] = (-1.0*param[ifreefac+firstpar+nregressors])*PDF[iterm]/diffCDF; 
		// gradient of factor loading (alpha)
		tmpgrad[numfac+nregressors+ifreefac] = (-1.0*fac[ifac])*PDF[iterm]/diffCDF;
		
		ifreefac++;
	      }
	    }
	  }
	  
	  // Gradient for X's
	  for (int ireg = 0; ireg < nregressors  ; ireg++) {
	    tmpgrad[ireg+numfac] = (-1.0*data[iobs_offset+regressors[ireg]])*PDF[iterm]/diffCDF;
	  }

	  // Gradient for thresholds
	  // Example: four categories 0 1 2
	  //iterm = 1 --> 0-obschoice-1
	  //iterm = 0 --> 1-obschoice-2
	  int thres_offset = npar - (numchoice - 1);
	  int maxthresloop = obsCat;
	  if (iterm==0) maxthresloop--;
	  for (int ithres = 0 ; ithres < maxthresloop ; ithres++)  {
	    tmpgrad[thres_offset + ithres] = PDF[iterm]/diffCDF;
	  }
	  
	  Int_t obsSign = (iterm==0) ? -1.0 : 1.0;
	  //Add the tmp gradient and hessian to the totals
	  for (int i = 0; i < npar ; i++) {
	    modEval[i+1] += obsSign*tmpgrad[i];
	  }	  
	} // Check for end categories
      } // loop over two terms


      // *****************************
      // Calculate the Hessian
      if (flag==3) {
	for (int iterm = 0; iterm < 2 ; iterm++) {
	  // skip if they are end categories
	  if ( ((obsCat>1)&&(iterm==0)) || ((obsCat<numchoice)&&(iterm==1)) ) {
	    
	    std::vector<double> tmphess;
	    tmphess.clear();
	    tmphess.resize(npar*npar,1.0);
	    
	    Double_t obsSign = (iterm==0) ? -1.0 : 1.0;

	    // First do the factor-specific gradient and hessian terms (theta and alpha)
	    ifreefac = 0;
	    for (int ifac = 0 ; ifac < numfac ; ifac++) {
	      //No normalizations:
	      if (facnorm.size()==0) {
	      
		// lambda^L(theta) - lambda^Prob(theta) (row) 
		for (int j = ifac; j < npar ; j++) tmphess[ifac*npar+j] *= -Z[iterm]*(-1.0*param[ifreefac+firstpar+nregressors]) - modEval[1+ifac];
		// dZ/dtheta_i (col) 
		for (int j = 0; j <= ifac ; j++) tmphess[j*npar+ifac] *= -1.0*param[ifreefac+firstpar+nregressors];
		
		// alpha_i index
		Int_t index = numfac+nregressors+ifreefac;
		// lambda^L(alpha) - lambda^Prob(alpha) (row)
		for (int j = index; j < npar ; j++) tmphess[index*npar+j] *= -Z[iterm]*(-1.0*fac[ifac]) - modEval[1+numfac+nregressors+ifreefac];
		// dZ/dalpha (col)
		for (int j = 0; j <= index ; j++) tmphess[j*npar+index] *= -1.0*fac[ifac];
		
		ifreefac++;
	      }
	      else {
		if (facnorm[ifac]>-9998.0) {
		  // lambda^L(theta) - lambda^Prob(theta) (row)  
		  for (int j = ifac; j < npar ; j++) tmphess[ifac*npar+j] *= -Z[iterm]*(-1.0*facnorm[ifac]) - modEval[1+ifac];
		  // dZ/dtheta_i (col)
		  for (int j = 0; j <= ifac ; j++) tmphess[j*npar+ifac] *= -1.0*facnorm[ifac];
		}
		else {
		  // lambda^L(theta) - lambda^Prob(theta) (row)  
		  for (int j = ifac; j < npar ; j++) tmphess[ifac*npar+j] *= -Z[iterm]*(-1.0*param[ifreefac+firstpar+nregressors]) - modEval[1+ifac];
		  // dZ/dtheta_i (col)
		  for (int j = 0; j <= ifac ; j++) tmphess[j*npar+ifac] *= -1.0*param[ifreefac+firstpar+nregressors];
		  
		  // alpha_i index
		  Int_t index = numfac+nregressors+ifreefac;
		  // lambda^L(alpha) - lambda^Prob(alpha) (row)
		  for (int j = index; j < npar ; j++) tmphess[index*npar+j] *= -Z[iterm]*(-1.0*fac[ifac]) - modEval[1+numfac+nregressors+ifreefac];
		  // dZ/dalpha (col)
		  for (int j = 0; j <= index ; j++) tmphess[j*npar+index] *= -1.0*fac[ifac];
		
		  ifreefac++;
		}
	      }
	    }
	  
	    // Hessian for X's
	    for (int ireg = 0; ireg < nregressors  ; ireg++) {
	      

	      Int_t index = numfac+ireg;
	      // lambda^L(X) - lambda^Prob(X) (row)  
	      for (int j = index; j < npar ; j++) tmphess[index*npar+j] *= -Z[iterm]*(-1.0*data[iobs_offset+regressors[ireg]]) - modEval[1+ireg+numfac];
	      // dZ/dX_i (col)
	      for (int j = 0; j <= index ; j++) tmphess[j*npar+index] *= -1.0*data[iobs_offset+regressors[ireg]];
	    }

	    // Hessian for thresholds
	    // Example: four categories 0 1 2
	    //iterm = 1 --> 0-obschoice-1
	    //iterm = 0 --> 1-obschoice-2
	    int thres_offset = npar - (numchoice - 1);
	    int maxthresloop = obsCat;
	    if (iterm==0) maxthresloop--;
	    for (int ithres = 0 ; ithres < numchoice-1 ; ithres++)  {
	      //	    for (int ithres = 0 ; ithres < maxthresloop ; ithres++)  {
	      
	      Int_t index = thres_offset + ithres;
	      if (ithres<maxthresloop) {
		// lambda^L(chi) - lambda^Prob(chi) (row)   
		for (int j = index; j < npar ; j++) tmphess[index*npar+j] *= -Z[iterm] - modEval[1+thres_offset + ithres];
	      }
 	      else {
 		// lambda^L(chi) - lambda^Prob(chi) (row)   
 		for (int j = index; j < npar ; j++) tmphess[index*npar+j] *= - modEval[1+thres_offset + ithres];
 		// dZ/dchi_i (col) 
 		for (int j = 0; j <= index ; j++) tmphess[j*npar+index] = 0.0; 
 	      }
	    }
	  
	    //Need to add dZ/dtheta dalpha
	    ifreefac = 0;
	    for (int i = 0 ; i < numfac ; i++) {
	      if (facnorm.size()==0) {
		Int_t index = numfac+nregressors+ifreefac;
		tmphess[i*npar+index] += -1.0;
		ifreefac++;
	      }
	      else if (facnorm[i]<-9998.0) {
		Int_t index = numfac+nregressors+ifreefac;
		tmphess[i*npar+index] += -1.0;
		ifreefac++;
	      }
	    }

	    //Add the tmp gradient and hessian to the totals
	    for (int i = 0; i < npar ; i++) {
	      for (int j = i ; j < npar ; j++) { 
		// Finally multiply by obsSign*PDF[iterm]/diffCDF
		hess[i*npar+j] += obsSign*tmphess[i*npar+j]*PDF[iterm]/diffCDF;
		//if (iobs_offset<=52) printf("Putting the pieces together (%5d,%3d,%3d): obsSign=%d, hess=%8f, pdf=%8f, modEval=%8f\n",iobs_offset,i,j,int(obsSign),tmphess[i*npar+j],PDF[iterm],diffCDF);
	      }
	    }
	  
	  } // Check for end categories
	} // loop over two terms
      } // flag==3

    } // check (flag>=2)
    //    cout << "Returning hessian size=" << hess.size() << "\n";
    // check for NaN
//     for (UInt_t i = 0; i <  modEval.size() ; i++) 
//       if ( std::isnan( modEval[i])) cout << "Found the NaN!! model type 2, element " << i << " modEval[0] = " << modEval[0] << " expess=" << expres[0] << "\n";

    return;
    
    //     if ( int(data[iobs_offset+outcome])==1) return ROOT::Math::normal_cdf(expres[0]);
    //     else if ( int(data[iobs_offset+outcome])==0) return ROOT::Math::normal_cdf(-expres[0]);
  }
  else if (modtype==3) {

    // Was state observed? (1,2,...,numchoice)
    // We'll use (0,1,...,numchoice-1) for the variable obsCat
    // First set density to zero so we can sum up the possible ranked choices
    modEval[0]=0.0;

    for (int irank = 0 ; irank <= numrank ; irank++) {
      Int_t obsCat = -1;
      for (int icat = 1 ; icat <= numchoice ; icat++) if (icat == int(data[iobs_offset + outcome + irank])) obsCat = icat-1;
      
      //Individuals may not use all rankings, so we only check the first one:
      if ((obsCat==-1)&&(irank==0)) {
	  cout << "ERROR (TModel::Eval): Found invalid number for logit outcome!"
	       << " In model " <<  this->GetTitle() << "\n"
	       << "Looking at outcome #" << outcome 
	       << ", Missing (" << missing << ")=" << data[iobs_offset+missing]
	       << ", obs#:" << iobs_offset/379.0
	       << ", Value=" << data[iobs_offset+outcome] << endl;
	  
	  assert(0);
      }

      if (obsCat>-1) {
	UInt_t nparamchoice = nregressors + numfac;

	vector<double> rankedChoiceCorr(numchoice - 1, 0.0);
	if (ranksharevar >= 0 )  {
	  for (int icat = 0 ; icat < numchoice-1; icat++) {
	    if (data[iobs_offset + ranksharevar + (numchoice-1)*irank + icat]> -9998.0) {
	      rankedChoiceCorr[icat] = data[iobs_offset + ranksharevar + (numchoice-1)*irank + icat];
	    }
	  }
	}
    
	double logitdenom = 1.0;
	for (int icat = 1 ; icat < numchoice; icat++) {
	  logitdenom += exp(expres[icat-1] + rankedChoiceCorr[icat-1]);
	}

	double dens = 1.0/logitdenom;

	if (obsCat>0) {
	  dens *=  exp(expres[obsCat-1] + rankedChoiceCorr[icat-1]);
	}
    
	modEval[0] += dens;
    
	// Now find the derivatives: 
	if (flag>=2) {
	  vector<double> pdf(numchoice,1.0/logitdenom);
	  for (int icat = 1 ; icat < numchoice ; icat++) {
	    pdf[icat] *= exp(expres[icat-1] + rankedChoiceCorr[icat-1]);
	  }

	  vector<double> logitgrad;
	  if (flag==3) {
	    logitgrad.resize(numchoice*npar,0.0);
	  }

	  // gradient of factor parameters
	  for (int ifac = 0 ; ifac < numfac ; ifac++) {
	    
	    //***************
	    // gradient of factor-specific parameters (variance, mean, weights, etc)
	    //obsCat term:
	    if (obsCat>0) modEval[1+ifac] += param[firstpar+(obsCat-1)*nparamchoice+nregressors+ifac];
	    if (flag==3) {
	      for (int jcat = 1 ; jcat < numchoice ; jcat++) logitgrad[jcat*npar + ifac] += param[firstpar+(jcat-1)*nparamchoice+nregressors+ifac];
	    }
	    
	    //no parameters for choice=0 (i.e. Z(icat=0) = 0)
	    for (int icat = 1 ; icat < numchoice ; icat++) {
	      modEval[1+ifac] += -pdf[icat]*param[firstpar+(icat-1)*nparamchoice+nregressors+ifac];
	      if (flag==3) {
		for (int jcat = 0 ; jcat < numchoice ; jcat++) logitgrad[jcat*npar + ifac] += -pdf[icat]*param[firstpar+(icat-1)*nparamchoice+nregressors+ifac];
	      }
	    }

	    //****************
	    //gradient of factor loading (alpha) 
	    //obsCat term:
	    if (obsCat>0) modEval[1+numfac+(obsCat-1)*nparamchoice+nregressors+ifac] += fac[ifac];
	    if (flag==3) {
	      for (int jcat = 1 ; jcat < numchoice ; jcat++) logitgrad[jcat*npar + numfac + (jcat-1)*nparamchoice + nregressors + ifac] +=  fac[ifac];
	    }
	    
	    
	    //no parameters for choice=0 (i.e. Z(icat=0) = 0)
	    for (int icat = 1 ; icat < numchoice ; icat++) {
	      modEval[1+numfac+(icat-1)*nparamchoice+nregressors+ifac] += -pdf[icat]*fac[ifac];
	      if (flag==3) {
		for (int jcat = 0 ; jcat < numchoice ; jcat++) logitgrad[jcat*npar + numfac + (icat-1)*nparamchoice + nregressors + ifac] += -pdf[icat]*fac[ifac];
	      }
	    }
	  }
	  
	  //****************
	  //gradient of betas
	  for (int ireg = 0; ireg < nregressors  ; ireg++) {
	    
	    //obsCat term:
	    if (obsCat>0) modEval[1+numfac+(obsCat-1)*nparamchoice+ireg] += data[iobs_offset+regressors[ireg]];
	    if (flag==3) {
	      for (int jcat = 1 ; jcat < numchoice ; jcat++) logitgrad[jcat*npar + numfac + (jcat-1)*nparamchoice + ireg] +=  data[iobs_offset+regressors[ireg]];
	    }
	    
	    //no parameters for choice=0 (i.e. Z(icat=0) = 0)
	    for (int icat = 1 ; icat < numchoice ; icat++) {
	      modEval[1+numfac+(icat-1)*nparamchoice+ireg] += -pdf[icat]*data[iobs_offset+regressors[ireg]];
	      if (flag==3) {
		for (int jcat = 0 ; jcat < numchoice ; jcat++) logitgrad[jcat*npar + numfac + (icat-1)*nparamchoice + ireg] += -pdf[icat]*data[iobs_offset+regressors[ireg]];
	      }
	    }
	  }
	  
	  // Now Calculate Hessian
	  
	  if (flag==3) {
	    
	    
	    hess.resize(npar*npar,0.0);
	    //	for (int i = 0; i < npar ; i++) {
	    //	  for (int j = i ; j < npar ; j++) {
	    //	    hess[i*npar+j] = 1.0;
	    //	  }
	    //	}



	    // First do second order derivatives (for now only dZ/dalpha dtheta terms are non-zero)
	    //Need to add dZ/dtheta dalpha
	    if (obsCat>0) {
	      for (int ifac = 0 ; ifac < numfac ; ifac++) {
		Int_t index = numfac+(obsCat-1)*nparamchoice+nregressors+ifac;
		hess[ifac*npar+index] += 1.0;
	      }
	    }
	    
	    //second, second-order derivative term 
	    for (int icat = 1 ; icat < numchoice ; icat++) {
	      for (int ifac = 0 ; ifac < numfac ; ifac++) {
		Int_t index = numfac+(icat-1)*nparamchoice+nregressors+ifac;
		hess[ifac*npar+index] += -pdf[icat]*1.0;
	      }
	    }
	    
	    //Now the first-order derivative Hessian terms
	    // dtheta d....
	    for (int ifac = 0 ; ifac < numfac ; ifac++) {
	      
	      // loop over g categories (see notes)
	      for (int icat = 1 ; icat < numchoice ; icat++) {
		
		// dtheta dtheta
		for (int jfac = ifac ; jfac < numfac ; jfac++) {
		  hess[ifac*npar+jfac] += -pdf[icat]*logitgrad[icat*npar + jfac]*param[firstpar+(icat-1)*nparamchoice+nregressors+ifac];
		}
		
		// loop over parameters in each category
		for (int jcat = 1 ; jcat < numchoice ; jcat++) {
		  
		  //dtheta dalpha
		  for (int jfac = 0 ; jfac < numfac ; jfac++) {
		    Int_t index = numfac+(jcat-1)*nparamchoice+nregressors+jfac;
		    hess[ifac*npar+index] += -pdf[icat]*logitgrad[icat*npar + index]*param[firstpar+(icat-1)*nparamchoice+nregressors+ifac];
		  }
		  
		  //dtheta dbeta
		  for (int ireg = 0; ireg < nregressors  ; ireg++) {
		    Int_t index = numfac+(jcat-1)*nparamchoice+ireg;
		    hess[ifac*npar+index] += -pdf[icat]*logitgrad[icat*npar + index]*param[firstpar+(icat-1)*nparamchoice+nregressors+ifac];
		  }
		}
	      }
	    }
	    
	    // loop over row categories (see notes)
	    for (int icat = 1 ; icat < numchoice ; icat++) {
	      // loop over parameters in each category
	      for (int jcat = icat ; jcat < numchoice ; jcat++) {
		
		//dbeta dbeta
		for (int ireg = 0; ireg < nregressors  ; ireg++) {
		  for (int jreg = 0; jreg < nregressors  ; jreg++) {
		    if ( (jcat>icat) || (jreg>=ireg) ) {
		      Int_t index1 = numfac+(icat-1)*nparamchoice+ireg;	      
		      Int_t index2 = numfac+(jcat-1)*nparamchoice+jreg;	      
		      hess[index1*npar+index2] += -pdf[icat]*logitgrad[icat*npar + index2]*data[iobs_offset+regressors[ireg]];
		    }
		  }
		}
		//dbeta dalpha
		for (int ireg = 0; ireg < nregressors  ; ireg++) {
		  for (int jfac = 0 ; jfac < numfac ; jfac++) {
		    Int_t index1 = numfac+(icat-1)*nparamchoice+ireg;	      
		    Int_t index2 = numfac+(jcat-1)*nparamchoice+nregressors+jfac;	      
		    hess[index1*npar+index2] += -pdf[icat]*logitgrad[icat*npar + index2]*data[iobs_offset+regressors[ireg]];
		  }
		}
		
		//dalpha dbeta
		if (jcat>icat) {
		  for (int ifac = 0 ; ifac < numfac ; ifac++) {
		    for (int jreg = 0; jreg < nregressors  ; jreg++) {
		      Int_t index1 = numfac+(icat-1)*nparamchoice+nregressors+ifac;	      
		      Int_t index2 = numfac+(jcat-1)*nparamchoice+jreg;	      
		      hess[index1*npar+index2] += -pdf[icat]*logitgrad[icat*npar + index2]*fac[ifac];
		    }
		  }
		}
		
		//dalpha dalpha
		for (int ifac = 0 ; ifac < numfac ; ifac++) {
		  for (int jfac = 0; jfac < numfac  ; jfac++) {
		    if ( (jcat>icat) || (jfac>=ifac) ) {
		      Int_t index1 = numfac+(icat-1)*nparamchoice+nregressors+ifac;	      
		      Int_t index2 = numfac+(jcat-1)*nparamchoice+nregressors+jfac;	      
		      hess[index1*npar+index2] += -pdf[icat]*logitgrad[icat*npar + index2]*fac[ifac];
		    }
		  }
		}
	    
	      } //jcat
	    } //icat
	  } // flag==3
	} //flag>=2
      } // obsCat>-1
    } //irank
  } //modtype==3

    //    cout << "Returning hessian size=" << hess.size() << "\n";
    // check for NaN
//     for (UInt_t i = 0; i <  modEval.size() ; i++) 
//       if ( std::isnan( modEval[i])) cout << "Found the NaN!! model type 2, element " << i << " modEval[0] = " << modEval[0] << " expess=" << expres[0] << "\n";

    return;
    
    //     if ( int(data[iobs_offset+outcome])==1) return ROOT::Math::normal_cdf(expres[0]);
    //     else if ( int(data[iobs_offset+outcome])==0) return ROOT::Math::normal_cdf(-expres[0]);


//  cout << "ERROR (TModel::Eval): Non-supported model!!"
//       << "Looking at model #" << modtype << endl;
//  assert(0);

}

void TModel::Sim(UInt_t iobs_offset, const std::vector<Double_t> & data, std::vector<TModel> & models, const std::vector<Double_t> & param, UInt_t firstpar, const std::vector <Double_t> & fac, FILE * pFile, UInt_t gof)
{

  Int_t nout = 1;
  if (splitsim>-1) nout=2;

  for (Int_t isplit = 0 ; isplit < nout ; isplit++) {

    // for now just check if regressors are not -9999:
    Int_t this_missing = 0;

    if (missing>-1) {
      if (data[iobs_offset+missing]==0) this_missing = 1;
    }

    for (int ireg = 0; ireg < nregressors  ; ireg++) {

      Int_t thisendog = 0;
      if (gof==0) {
	for (UInt_t iendogvar = 0 ; iendogvar < endogRegList.size() ; iendogvar++) {
	  if (regressors[ireg] == endogRegList[iendogvar]) {
	    thisendog = 1;
	    if ( (endogModelList[iendogvar]>=0)&&(models.at(endogModelList[iendogvar]).GetSimResult()<-9998) ) this_missing = 2;
	  }
	}
      }
      if (thisendog==0) {
	if (data[iobs_offset+regressors[ireg]]<-9998) this_missing = 2;
      }
    }
    

    fprintf(pFile,", %10d",this_missing);
    
    Int_t numlogitchoice = 2;

    if ((modtype==3)&&(numchoice>2)) {
      numlogitchoice = numchoice;
    }


    if (this_missing >1) {

      //print out detailed variables (Vobs, Vend, Vfac0..VfacN, eps)
      if (detailsim) {
	// Vobs, eps, Vfac#
	int ndetailvar = 2 + numfac;
	if ((endogRegList.size()>0)&&(gof==0)) ndetailvar++;
	for (int ichoice = 2 ; ichoice <=numlogitchoice; ichoice++) {
	  for (int i = 0; i < ndetailvar; i++)  fprintf(pFile,", %10d",-9999);
	}
      }

      //print out prob for each choice
      if (modtype==2)  fprintf(pFile,", %10d",-9999);
      if ((modtype==3)||(modtype==4)) {
	for (int ichoice = 2 ; ichoice <=numchoice; ichoice++) fprintf(pFile,", %10d",-9999);
      }

      //print out predicted value
      fprintf(pFile,", %10d",-9999);
      simresult = -9999.0;
    }
    else {

      vector<double> expres(numlogitchoice-1,0.0);

      vector<double> vendog(numlogitchoice-1,0.0);

      UInt_t ifreefac = 0;

      for (int ichoice = 0; ichoice < numlogitchoice-1; ichoice++) {
	
	ifreefac = 0;
	UInt_t nparamchoice = nregressors + numfac;
	
	for (int ireg = 0; ireg < nregressors  ; ireg++) {
	  if (regressors[ireg]!=splitsim) {
	    
	    double thisterm = param[ireg+firstpar+ichoice*nparamchoice]*data[iobs_offset+regressors[ireg]];

	    // Don't use endogenous variables if doing goodness of fit
	    if (gof==0) {
	      for (UInt_t iendogvar = 0 ; iendogvar < endogRegList.size() ; iendogvar++) {
		if (regressors[ireg] == endogRegList[iendogvar]) {
		  
		  // process linear and probit models
		  if ( (endogChoiceList[iendogvar]==-1) && (endogModelList[iendogvar]>=0) ) {
		    thisterm = param[ireg+firstpar+ichoice*nparamchoice]*models.at(endogModelList[iendogvar]).GetSimResult();
		    vendog[ichoice] += thisterm;
		  }
		  // process multinomial choice models or Xtiles
		  else if ( (endogChoiceList[iendogvar]>-1) && (endogModelList[iendogvar]>-1) ) {
		    int simchoice = models.at(endogModelList[iendogvar]).GetSimResult();
		    if (endogChoiceList[iendogvar] == simchoice) {
		      thisterm = param[ireg+firstpar+ichoice*nparamchoice];
		      vendog[ichoice] += thisterm;
		    }
		    else thisterm = 0.0;
		  }
		  // don't do anything for missing indicators (set to zero)
		  else {
		    //require model >= 0 as model= -1 implies missing indicator which we ignore in simulation
		  // ignore missing indicators
		    thisterm = 0.0;
		}
		}
	      } //loop over endogenous variables and see if there is a match
	    }
	    expres[ichoice] += thisterm;
	      
	      //	     printf("Model %25s using %25s: obs=%6d reg=%4d beta=%8.3f endog=%8.3f vend=%8.3f\n",GetName(),models.at(endogModel).GetName(),iobs_offset,regressors[ireg],param[i+firstpar],models.at(endogModel).GetSimResult(),vendog);
	  }
	  else expres[ichoice] += param[ireg+firstpar+ichoice*nparamchoice]*isplit;
	}
	if (detailsim) {
	  if ( (endogRegList.size()==0)||(gof==1) ) {
	    fprintf(pFile,", %10.5f",expres[ichoice]);
	  }
	  else {
	    fprintf(pFile,", %10.5f",expres[ichoice]-vendog[ichoice]);
	    fprintf(pFile,", %10.5f",vendog[ichoice]);
	  }
	}      
	Double_t fac_comp = 0.0;
	for (int i = 0 ; i < numfac ; i++) {
	  if (facnorm.size()==0) {
	    fac_comp += param[ifreefac+firstpar+ichoice*nparamchoice+nregressors]*fac.at(i);
	    if (detailsim) fprintf(pFile,", %10.5f",param[ifreefac+firstpar+ichoice*nparamchoice+nregressors]*fac.at(i));
	    
	    ifreefac++;
	  }
	  else {
	    if (facnorm[i]>-9998) {
	      fac_comp += facnorm[i]*fac.at(i);
	      if (detailsim) fprintf(pFile,", %10.5f",facnorm[i]*fac.at(i));
	    }
	    else {
	      fac_comp += param[ifreefac+firstpar+ichoice*nparamchoice+nregressors]*fac.at(i);
	      if (detailsim) fprintf(pFile,", %10.5f",param[ifreefac+firstpar+ichoice*nparamchoice+nregressors]*fac.at(i));
	      ifreefac++;
	    }
	  }
	}
	//      if (detailsim) fprintf(pFile,", %10.5f",fac_comp);
	expres[ichoice] += fac_comp;
      }

      if (modtype==1) {
	Double_t sigma = fabs(param[firstpar+nregressors+ifreefac]);
	Double_t eps = gRandom->Gaus(0.0,sigma);
	if (detailsim) fprintf(pFile,", %10.5f",eps);
	fprintf(pFile,", %10.5f",expres[0]+eps);
	simresult = expres[0]+eps;

	//Save Xtile if Xtiles are used in other models
	if (endogXtiles.size()>0) {
	  simresult = 1;
	  for(UInt_t ixtile = 0 ; ixtile < endogXtiles.size() ; ixtile++) {
	    if (expres[0]+eps > endogXtiles[ixtile]) simresult = ixtile+2;
	  }
	}
	if (gof==1) {
	  simresult = - expres[0];
	  if (outcome>-1) simresult += data[iobs_offset+outcome];
	  //	  fprintf(pFile,", %10.5f",simresult);
	}
      }
      else if (modtype==2) {
	Double_t eps = gRandom->Gaus(0.0,1.0);
	Double_t prb = ROOT::Math::normal_cdf(expres[0]);
	//	if (detailsim) fprintf(pFile,", %10.5f",prb);
	if (detailsim) fprintf(pFile,", %10.5f",eps);
	fprintf(pFile,", %10.5f",prb);
	if ((expres[0]+eps)>0) { 
	  fprintf(pFile,", %10d",1);
	  simresult = 1.0;
	}
	else {
	  fprintf(pFile,", %10d",0);
	  simresult = 0.0;
	}
	if (gof==1) {
	  if ((expres[0])>0) {
	    simresult = (int(data[iobs_offset+outcome])==1);
	  }
	  else {
	    simresult = (int(data[iobs_offset+outcome])==0);
	  }
	  //	  fprintf(pFile,", %10d",int(simresult));
	}
      }
      else if (modtype==3) {

	std::extreme_value_distribution<double> EV1Dist(0,1.0);

	//Simulate choice
	int thischoice = 1;

	// Draw extreme error value for choice 1
	double eps1 = EV1Dist(mt);
        double maxval = eps1;

	// Get value for choice in data for GoF
	double obsval = 0.0;
        Int_t obsrank = 0;	
	if (int(data[iobs_offset+outcome])>1) obsval = expres.at(int(data[iobs_offset+outcome])-2);
	
	// See if any of the other choices are larger
	for (int icat = 1 ; icat < numchoice; icat++) {
	  double eps = EV1Dist(mt);

	  //print out difference in epsilons
	  if (detailsim) fprintf(pFile,", %10.5f",eps-eps1);

	  double thisval = expres[icat-1] + eps;

	  //Get rank of observed choice according to model
	  if ( (int(data[iobs_offset+outcome]) != icat+1) && (expres[icat-1] > obsval)) obsrank++;
	  
	  if (thisval > maxval ) {
	    maxval = thisval;
	    thischoice = icat+1;	    
	  }
	}

	//print out probabilities of choices...
	double logitdenom = 1.0;
	for (int icat = 1 ; icat < numchoice; icat++) {
	  logitdenom += exp(expres[icat-1]);
	}
	for (int icat = 1 ; icat < numchoice; icat++) fprintf(pFile,", %10.5f",exp(expres[icat-1])/logitdenom);

	fprintf(pFile,", %10d",thischoice);
	simresult = Double_t(thischoice);

	if (gof==1) {
	  simresult = Double_t(obsrank);
	  //	  fprintf(pFile,", %10d",obsrank);
	}
	  
      }
      else if (modtype==4) {

	Double_t eps = gRandom->Gaus(0.0,1.0);
	if (detailsim) fprintf(pFile,", %10.5f",eps);

	double threshold[2];
	for (int ichoice = 1 ; ichoice < numchoice; ichoice++) {
	  
	  threshold[0] = param[firstpar+nregressors+ifreefac];
	  threshold[1] = param[firstpar+nregressors+ifreefac];
	  for (int icat = 2 ; icat <= ichoice ; icat++) {
	    if (icat<ichoice) threshold[0] += abs(param[firstpar+nregressors+ifreefac+icat-1]);
	    if (icat<numchoice) threshold[1] += abs(param[firstpar+nregressors+ifreefac+icat-1]);
	  }
	  double CDF[2] = {0.0, 1.0};

	  if (ichoice>1) {
	    CDF[0] = ROOT::Math::normal_cdf(threshold[0] - expres[0]);
	  }
	  if (ichoice<numchoice) {
	    CDF[1] = ROOT::Math::normal_cdf(threshold[1] - expres[0]);
	  }

	  //Print out probability of this choice
	  fprintf(pFile,", %10.5f",CDF[1] - CDF[0]);
	}

	Int_t choice = 1;
	Double_t thisthreshold = param[firstpar+nregressors+ifreefac];
	for (int icat = 2 ; icat < numchoice ; icat++) { 
	  if ( (expres[0]+eps > thisthreshold) && (expres[0]+eps < thisthreshold + abs(param[firstpar+nregressors+ifreefac+icat-1])) ){ 
	    choice = icat;
	  }
	  thisthreshold += abs(param[firstpar+nregressors+ifreefac+icat-1]);
	}
	if (expres[0]+eps > thisthreshold) choice = numchoice;

	fprintf(pFile,", %10d",choice);
	simresult = Double_t(choice);
      }
      else {
	cout << "ERROR (TModel::Eval): Non-supported model!!"
	     << "Looking at model #" << modtype << endl;
	assert(0);
      }
    }
  }
}


Double_t TModel::GetPdf(UInt_t iobs_offset, const std::vector<Double_t> & data, const std::vector<Double_t> & param, UInt_t firstpar, const std::vector <Double_t> & fac)
{

  if (modtype!=2) {
    cout << "ERROR (TModel::GetProb): No need to estimate Marginal Effect for model other than probit!!\n";
    assert(0);
  }
  else {
    // for now just check if regressors are not -9999:
    Int_t this_missing = 0;
    
    if (missing>-1) {
      if (data[iobs_offset+missing]==0) this_missing = 1;
    }
    
    if (this_missing==0) {
      for (int i = 0; i < nregressors  ; i++) {
	if (data[iobs_offset+regressors[i]]<-9998) this_missing = 2;
      }
    }
    if (this_missing == 1) {
      return -1.0;
    }
    else if (this_missing==2) {
      cout << "ERROR (TModel::GetProb): Estimating marginal effect with missing covariate!!\n";
      assert(0);
    }
    else {
      double expres = 0.0;
      
      for (int i = 0; i < nregressors  ; i++) {
	expres += param[i+firstpar]*data[iobs_offset+regressors[i]];
      }
      
      UInt_t ifreefac = 0;
      Double_t fac_comp = 0.0;
      for (int i = 0 ; i < numfac ; i++) {
	if (facnorm.size()==0) {
	  fac_comp += param[ifreefac+firstpar+nregressors]*fac.at(i);	  
	  ifreefac++;
	}
	else {
	  if (facnorm[i]>-9998) {
	    fac_comp += facnorm[i]*fac.at(i);
	  }
	  else {
	    fac_comp += param[ifreefac+firstpar+nregressors]*fac.at(i);
	    ifreefac++;
	  }
	}
      }
      expres += fac_comp;

      Double_t prb = ROOT::Math::normal_pdf(expres);
      return prb;

    } //not missing

  } //Not probit type model

} 
