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

//ClassImp(TModel)

TModel::TModel(const char *name, const char *title, Int_t modeltype, Int_t modelgroup, Int_t prntgroup, std::vector<Int_t> & moddata, Int_t nfac, Double_t * thisnormfac, UInt_t nchoice)
  : TNamed(name,title)
{
  ignore = 0;
  splitsim = -1;
  detailsim = 0;
  modtype = modeltype;
  modgroup = modelgroup;
  printgroup = prntgroup;
  numchoice = nchoice;
  outcome = moddata[0];
  missing = moddata[1];
  nregressors = moddata.size()-2;
  regressors.reserve(nregressors);
  for (int i = 0 ; i < nregressors ; i++) regressors.push_back(moddata[i+2]);

//  printf("Regressors: ");
//  for (int i = 0 ; i < nregressors ; i++) {
//    if (i!=0) cout << ", "; 
//    cout << regressors[i];
//    //    cout << vartab[regressors[i]].Data();
//  }
//  printf("\n");


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

  printf(" Outcome: %12s\n",vartab.at(outcome).Data());
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
  //  cout << "Entering model, flag=" << flag <<"\n";

  //  std::vector<Double_t> modEval;
  if (flag >= 2) {

    //General parameters in gradient
    // 1 (likelihood) 2*numfac (df/dtheta and df/dalpha) df/dbeta
    Int_t ngrad = 1+2*numfac+nregressors;

    //Model specific parameters:
    if (modtype==1) ngrad += 1; // variance of error term
    if ((modtype==3)&&(numchoice>2)) ngrad += (numchoice-2)*(2*numfac+nregressors); // multinomial logit
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

  vector<double> Z(1,0.0);
  if ((modtype==3)&&(numchoice>2)) Z.resize(numchoice-1,0.0);

  for (int ichoice = 0; ichoice < numchoice-1; ichoice++) {
    for (int i = 0; i < nregressors  ; i++) {
      Z[ichoice] += param[i+firstpar]*data[iobs_offset+regressors[i]];
    }
    
    UInt_t ifreefac = 0;
    for (int i = 0 ; i < numfac ; i++) {
      if (facnorm.size()==0) {
	Z[ichoice] += param[ifreefac+firstpar+nregressors]*fac[i];
	ifreefac++;
      }
      else {
	if (facnorm[i]>-9998) Z[ichoice] += facnorm[i]*fac[i];
	else {
	  Z[ichoice] += param[ifreefac+firstpar+nregressors]*fac[i];
	  ifreefac++;
	}
      }
    }
  } // loop over choices

  // nparameters: d/dtheta, d/dbeta, d/dalpha
  Int_t npar = numfac+nregressors+ifreefac;

  if ((modtype==3)&&(numchoice>2)) npar *= (numchoice-2);


  if (modtype==1) {
    Double_t Z = data[iobs_offset+outcome]-Z[0];
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
	  // gradient of factor-specific parameters (variance, mean, weights)
	  //	  modEval[i+1] = Z*(fac[i]*param[ifreefac+firstpar+nregressors]/param[i])/(sigma*sigma);
	  modEval[i+1] = Z*(param[ifreefac+firstpar+nregressors])/(sigma*sigma);
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
    modEval[0] = ROOT::Math::normal_cdf(obsSign*Z[0]);
//     if ((obsSign*Z[0]<-35.0)&&(flag!=2)) {
//       modEval[0] = 1.0e-50/fabs(Z[0]);
//       //      cout << "TModel: Z[0]s too small, artificially fixing it\n";
//     }
    
    // Now find the derivatives: 
    if (flag>=2) {
      Double_t Z = obsSign*Z[0];
      Double_t pdf = ROOT::Math::normal_pdf(obsSign*Z[0]);   //(obsSign here is not necessary)
      //      if (fabs(Z[0])>35.0) pdf = 1.0e-30/fabs(Z[0]);
      Double_t cdf = modEval[0];
      if (obsSign*Z[0]<-35.0) cdf = 1.0e-50;

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
// 	  if (obsSign*Z[0]<-35.0) modEval[1+numfac+nregressors+ifreefac] = -0.1*param[ifreefac+firstpar+nregressors]*param[ifreefac+firstpar+nregressors]/fabs(param[ifreefac+firstpar+nregressors]);

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
// 	    if (obsSign*Z[0]<-35.0) modEval[1+numfac+nregressors+ifreefac] = -0.1*param[ifreefac+firstpar+nregressors]*param[ifreefac+firstpar+nregressors]/fabs(param[ifreefac+firstpar+nregressors]);

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
//       if ( std::isnan( modEval[i])) cout << "Found the NaN!! model type 2, element " << i << " modEval[0] = " << modEval[0] << " expess=" << Z[0] << "\n";

    return;
    
    //     if ( int(data[iobs_offset+outcome])==1) return ROOT::Math::normal_cdf(Z[0]);
    //     else if ( int(data[iobs_offset+outcome])==0) return ROOT::Math::normal_cdf(-Z[0]);
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
      CDF[0] = ROOT::Math::normal_cdf(threshold[0] - Z[0]);
    }
    if (obsCat<numchoice) {
      CDF[1] = ROOT::Math::normal_cdf(threshold[1] - Z[0]);
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
	Z[0] = threshold[0] - Z[0];
	PDF[0] = ROOT::Math::normal_pdf(threshold[0] - Z[0]);
      }
      if (obsCat<numchoice) {
	Z[1] = threshold[1] - Z[0];
	PDF[1] = ROOT::Math::normal_pdf(threshold[1] - Z[0]);
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
//       if ( std::isnan( modEval[i])) cout << "Found the NaN!! model type 2, element " << i << " modEval[0] = " << modEval[0] << " expess=" << Z[0] << "\n";

    return;
    
    //     if ( int(data[iobs_offset+outcome])==1) return ROOT::Math::normal_cdf(Z[0]);
    //     else if ( int(data[iobs_offset+outcome])==0) return ROOT::Math::normal_cdf(-Z[0]);
  }
  else if (modtype==3) {
    if ( int(data[iobs_offset+outcome])==1) {
      modEval[0] =  exp(Z[0])/(1.0+exp(Z[0]));
    }
    else if ( int(data[iobs_offset+outcome])==0) {
      modEval[0] =  1.0/(1.0+exp(Z[0]));
    }
    cout << "ERROR (TModel::Eval): Found non-binary number for logit outcome!"
	 << " In model, " <<  this->GetTitle() << "\n"
 	 << "Looking at outcome #" << outcome 
 	 << ", Value=" << data[iobs_offset+outcome] << endl;
    assert(0);
  }

  cout << "ERROR (TModel::Eval): Non-supported model!!"
       << "Looking at model #" << modtype << endl;
  assert(0);

}

void TModel::Sim(UInt_t iobs_offset, const std::vector<Double_t> & data, const std::vector<Double_t> & param, UInt_t firstpar, const std::vector <Double_t> & fac, FILE * pFile)
//std::vector<Double_t> TModel::Eval(UInt_t iobs_offset, Double_t * data, std::vector<Double_t> param, UInt_t firstpar, std::vector<Double_t> fac)
{

  Int_t nout = 1;
  if (splitsim>-1) nout=2;

  for (Int_t isplit = 0 ; isplit < nout ; isplit++) {

    // for now just check if regressors are not -9999:
    Int_t this_missing = 0;

    if (missing>-1) {
      if (data[iobs_offset+missing]==0) this_missing = 1;
    }

    for (int i = 0; i < nregressors  ; i++) {
      if (data[iobs_offset+regressors[i]]<-9998) this_missing = 2;
    }
    

    fprintf(pFile,", %10d",this_missing);
    
    if (this_missing >1) {
      if (detailsim) for (int i = 0; i < 4; i++)   fprintf(pFile,", %10d",-9999);
      else {
	fprintf(pFile,", %10d",-9999);
	if (modtype==2) fprintf(pFile,", %10d",-9999);
      }
    }
    else {
      Double_t Z[0] = 0.0;
      for (int i = 0; i < nregressors  ; i++) {
	if (regressors[i]!=splitsim) {
	  Z[0] += param[i+firstpar]*data[iobs_offset+regressors[i]];
	}
	else Z[0] += param[i+firstpar]*isplit;
      }
      if (detailsim) fprintf(pFile,", %10.5f",Z[0]);
      
      
      UInt_t ifreefac = 0;
      Double_t fac_comp = 0.0;
      for (int i = 0 ; i < numfac ; i++) {
	if (facnorm.size()==0) {
	  fac_comp += param[ifreefac+firstpar+nregressors]*fac.at(i);
	  if (detailsim) fprintf(pFile,", %10.5f",param[ifreefac+firstpar+nregressors]*fac.at(i));

	  ifreefac++;
	}
	else {
	  if (facnorm[i]>-9998) {
	    fac_comp += facnorm[i]*fac.at(i);
	    if (detailsim) fprintf(pFile,", %10.5f",facnorm[i]*fac.at(i));
	  }
	  else {
	    fac_comp += param[ifreefac+firstpar+nregressors]*fac.at(i);
	    if (detailsim) fprintf(pFile,", %10.5f",param[ifreefac+firstpar+nregressors]*fac.at(i));
	    ifreefac++;
	  }
	}
      }
      //      if (detailsim) fprintf(pFile,", %10.5f",fac_comp);
      Z[0] += fac_comp;
    
      if (modtype==1) {
	Double_t sigma = fabs(param[firstpar+nregressors+ifreefac]);
	Double_t eps = gRandom->Gaus(0.0,sigma);
	if (detailsim) fprintf(pFile,", %10.5f",eps);
	fprintf(pFile,", %10.5f",Z[0]+eps);
      }
      else if (modtype==2) {
	Double_t eps = gRandom->Gaus(0.0,1.0);
	Double_t prb = ROOT::Math::normal_cdf(Z[0]);
	//	if (detailsim) fprintf(pFile,", %10.5f",prb);
	fprintf(pFile,", %10.5f",prb);
	if ((Z[0]+eps)>0) fprintf(pFile,", %10d",1);
	else fprintf(pFile,", %10d",0);
      }
      else if (modtype==3) {
	cout << "ERROR (TModel::Sim): Simulation of Logit not supported yet!!\n";
	assert(0);
      }
      else if (modtype==4) {
	Double_t eps = gRandom->Gaus(0.0,1.0);
	//	if (detailsim) fprintf(pFile,", %10.5f",prb);
	fprintf(pFile,", %10.5f",eps);
	Int_t choice = 1;
	Double_t threshold = param[firstpar+nregressors+ifreefac];
	for (int icat = 2 ; icat < numchoice ; icat++) { 
	  if ( (Z[0]+eps > threshold) && (Z[0]+eps < threshold + abs(param[firstpar+nregressors+ifreefac+icat-1])) ){ 
	    choice = icat;
	  }
	  threshold += abs(param[firstpar+nregressors+ifreefac+icat-1]);
	}
	if (Z[0]+eps > threshold) choice = numchoice;

	fprintf(pFile,", %10d",choice);
      }
      else {
	cout << "ERROR (TModel::Eval): Non-supported model!!"
	     << "Looking at model #" << modtype << endl;
	assert(0);
      }
    }
  }

}
