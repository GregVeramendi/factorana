#ifndef __CINT__

#include "Riostream.h"
#include "Rtypes.h"
#include "TROOT.h"
#include "TFile.h"

//#include "nlopt.h"

#include "TMinLkhd.hh"

#include "mpi.h"

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <sstream>


//using namespace ROOT::Math; 
using namespace std;


#endif

extern TMinLkhd mymodel;


int fexist(const char *filename );
void init(int mpirank);
void MinuitLkhdFcn(Int_t &npar, Double_t *gin, Double_t &f,
	       Double_t *par, Int_t iflag);

void LkhdFcn(Int_t &npar, Double_t *gin, Double_t &f,
		   Double_t *par, Int_t iflag, Double_t *hess);

int main(int argc, char *argv[]) { 
  init(0);

  TString workingdir = mymodel.GetWorkingDir();

  if ((mymodel.GetNSubSample()==0)&&(mymodel.GetNBootSample()==0)) {
    mymodel.Minimize(0);
  }
  else if (mymodel.GetNBootSample()>0) {

    //Check to see if bootstrap sample was passed as command line argument
    UInt_t bootSampleStart = mymodel.GetBootStrapStartSample();
    UInt_t bootSampleEnd = mymodel.GetNBootSample();
    
    if (argc>=2) {
      std::istringstream ss(argv[1]);
      if (ss >> bootSampleStart)
	{
	  printf("Bootstrap sample set at command line to %d\n",bootSampleStart);
	  bootSampleEnd = bootSampleStart + 1;	  
	}
      else printf("*****Error in specifying bootstrap sample using command line. Unable to convert to integer! Argument: %s\n",argv[1]);
    }

    if (argc>=3) {
      std::istringstream ss(argv[2]);
      UInt_t sim_nobs;
      if (ss >> sim_nobs)
	{
	  printf("Simulation size set at command line to %d\n",sim_nobs);
	  mymodel.SetSimNobs(sim_nobs);
	}
      else printf("*****Error in setting simulation size using command line. Unable to convert to integer! Argument: %s\n",argv[2]);
    }
    
    for (UInt_t isample =  bootSampleStart; isample < bootSampleEnd ; isample++) {

      Int_t nobs = mymodel.GetNobs();
      Int_t * obslist = new Int_t[nobs];

      TString filename;
      std::stringstream out;
      out << isample;
      filename = TString("newbootstrapsample_").Append(out.str()).Append(".txt");
      filename.Prepend(workingdir);

      if (fexist((char *)filename.Data())) {
	cout << "\n\n\n" 
	     << "*********************************************************************\n";
	cout << "*********************************************************************\n";
	cout << "***Processing bootstrap sample #" << isample << "************************************\n";
	cout << "*********************************************************************\n";
	cout << "*********************************************************************\n\n\n";

	ifstream in1;
	in1.open(filename.Data());
	for (Int_t iobs = 0 ; iobs < nobs; iobs++) {
	  in1 >> obslist[iobs];
	  if (!in1.good()) {
	    cout << "Problem reading bootstrap sample values\n";
	    //	  assert(0);
	    return 0;
	  }
	}
	in1.close();
	mymodel.SetBootStrap(obslist);
	mymodel.SetCurrentSample(isample);
	delete [] obslist;
	
	mymodel.Minimize(0);
      }
    } 
  }
  else {
    for (Int_t isample = 0 ; isample < mymodel.GetNSubSample() ; isample++) {
      cout << "\n\n\n *********************************************************************\n";
      cout << "***Processing sample #" << isample << "\n";
      Int_t nobs = mymodel.GetNobs();
      Int_t * obslist = new Int_t[nobs];

      TString filename;
      std::stringstream out;
      out << isample;
      filename = TString("subsample_").Append(out.str()).Append(".txt");
      filename.Prepend(workingdir);

      ifstream in1;
      in1.open(filename.Data());
      for (Int_t iobs = 0 ; iobs < nobs; iobs++) {
	in1 >> obslist[iobs];
	if (!in1.good()) {
	  cout << "Problem reading subsample values\n";
	  //	  assert(0);
	  return 0;
	}
      }
      in1.close();
      mymodel.SetSkipObs(obslist);
      mymodel.SetCurrentSample(isample);
      delete [] obslist;

      mymodel.Minimize(0);
    } 
  }
  return 0;
}

void MinuitLkhdFcn(Int_t &npar, Double_t *gin, Double_t &f,
		   Double_t *par, Int_t iflag) {
  Int_t flag = iflag;
  if (flag !=2) flag = 1;
  LkhdFcn(npar,gin,f,par,flag,NULL);
}

void LkhdFcn(Int_t &npar, Double_t *gin, Double_t &f,
	     Double_t *par, Int_t iflag, Double_t *hess) {
  f = 0;
  //  mymodel.EvalLkhd(f,gin,par,2,0,1);
  
  //  cout << "The flag is " << iflag << "\n";
  // iflag=2;
  mymodel.EvalLkhd(f,gin,hess,par,iflag,0,1);
  if ((mymodel.GetCounter())%100==0) printf("Iteration %5d Found logLkhd: %f \n",mymodel.GetCounter(),f);
  //  printf("Iteration %5d, flag=%2d: Found logLkhd: %f \n",mymodel.GetCounter(),iflag,f);
//   if (mymodel.GetCounter()<10) {
//     cout << "The flag is " << iflag << "\n";
//     cout << "Gradiant:" << endl;
//     for (UInt_t i = 0 ; i < npar; i++) cout << i << " "<< gin[i] << endl;
//   }
  return;
}
