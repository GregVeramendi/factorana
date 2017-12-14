#ifndef __CINT__

#include "Riostream.h"
#include "Rtypes.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cmath>
#include "TROOT.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"

#include "TMinLkhd.hh"

#include "mpi.h"

using namespace ROOT::Math; 
using namespace std;


#endif

extern TMinLkhd mymodel;

int fexist(const char *filename );
void init(int mpirank);
void MinuitLkhdFcn(Int_t &npar, Double_t *gin, Double_t &f,
	       Double_t *par, Int_t iflag);

void LkhdFcn(Int_t &npar, Double_t *gin, Double_t &f,
		   Double_t *par, Int_t iflag, Double_t *hess);

// void LkhdFcn(Int_t &npar, Double_t *gin, Double_t &f,
// 	       Double_t *par, Int_t iflag);

int main(int argc,char *argv[]) { 

  int    my_rank, np;
  int    namelen;
  char   processor_name[MPI_MAX_PROCESSOR_NAME];
  

  int already_initialised;

  MPI_Initialized(&already_initialised);
  if (!already_initialised) {
    MPI_Init(&argc,&argv);
  }

  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Get_processor_name(processor_name,&namelen);
  
  fprintf(stdout,"Initializing process %d of %d is on %s\n",
	  my_rank, np, processor_name);
  fflush(stdout);

  init(my_rank);
  mymodel.ResetFitInfo();
  TString workingdir = mymodel.GetWorkingDir();

  // Int_t source, dest, tag, finished, flag;
  Int_t finished, flag;
  Double_t thisprob, totalprob;
  Int_t num_param = mymodel.GetNparam();

  Double_t * par = new Double_t[num_param];
  Double_t * thisgrad = new Double_t[num_param];
  Double_t * totalgrad = new Double_t[num_param];
  Double_t * thishess = new Double_t[num_param*(num_param+1)/2];
  Double_t * totalhess = new Double_t[num_param*(num_param+1)/2];

  UInt_t nobs = mymodel.GetNobs();
  UInt_t nmod = mymodel.GetNmodels();

  Int_t * ignoremod = new Int_t[nmod];
  Int_t * skipobs = new Int_t[nobs];
  Int_t * bootstrapobs = new Int_t[nobs];
  Int_t * fpar = new Int_t[num_param];

//   Int_t newbootstrap;
//   Int_t newskipobs;
//   Int_t newfixpar;
//   Int_t newmodignore;
  
  if (my_rank!=0) {    
//     source = 0;
//     dest = 0;
//     tag = 0;
    MPI_Bcast(&finished, 1, MPI_INT, 0, MPI_COMM_WORLD);

    while (!finished) {

      //Check to see if any settings have changed:
      UInt_t newflag = 0;
      MPI_Bcast(&newflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
      mymodel.SetNewFlag(newflag); 

      if (newflag !=0 ) { 

      // Check to see if sub-sampling has changed
//       MPI_Bcast(&newskipobs, 1, MPI_INT, 0, MPI_COMM_WORLD);
//       if (newskipobs==1) {
	
	if (mymodel.GetNewSkipObs()==1) {
	  MPI_Bcast(skipobs, nobs, MPI_INT, 0, MPI_COMM_WORLD);
	  mymodel.SetSkipObs(skipobs);
	  
	  Int_t nskip = 0;
	  Int_t firstskip = -1;
	  for (UInt_t iobs = 0; iobs < nobs ; iobs++) {
	    nskip += skipobs[iobs];
	  if ((firstskip==-1)&&(skipobs[iobs]==1)) firstskip = iobs;
	  }
// 	  fprintf(stdout,"Process %d of %d is on %s: skipping %d,first=%d\n",
// 		  my_rank, np, processor_name, nskip,firstskip);
	}

      // Check to see if bootstrap has changed
//       MPI_Bcast(&newbootstrap, 1, MPI_INT, 0, MPI_COMM_WORLD);
//       if (newbootstrap==1) {
	
	if (mymodel.GetNewBootStrap()==1) {
	  MPI_Bcast(bootstrapobs, nobs, MPI_INT, 0, MPI_COMM_WORLD);
	  mymodel.SetBootStrap(bootstrapobs);
	  
// 	  fprintf(stdout,"Process %d of %d is on %s: BootStrap Set,first obs=%d, last obs=%d\n",
// 		  my_rank, np, processor_name, bootstrapobs[0],bootstrapobs[nobs-1]);
	}


      //Check to see if fixing of parameters has changed
//       MPI_Bcast(&newfixpar, 1, MPI_INT, 0, MPI_COMM_WORLD);
//       if (newfixpar==1) {

	if (mymodel.GetNewFixPar()==1) {
	  MPI_Bcast(fpar, num_param, MPI_INT, 0, MPI_COMM_WORLD);
	  mymodel.SetFixPar(fpar);
	}

 	// Int_t nfreeparam = mymodel.GetNfreeparam();
 	// if (mymodel.GetNewFixPar()==1) fprintf(stdout,"Process %d of %d is on %s: has %d free parameters\n",
 	// 			  my_rank, np, processor_name, nfreeparam);

	//Check to see if ignoring models has changed
	//       MPI_Bcast(&newmodignore, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//       if (newmodignore==1) {

	if (mymodel.GetNewIgnoreMod()==1) {
	  MPI_Bcast(ignoremod, nmod, MPI_INT, 0, MPI_COMM_WORLD);
	  mymodel.SetAllIgnoreMod(ignoremod);

// 	  Int_t nignoremod = mymodel.GetNignoremod();
// 	  fprintf(stdout,"Process %d of %d is on %s: ignoring %d models\n",my_rank, np, processor_name, nignoremod);
	}

	if (mymodel.GetNewEstNMix()==1) {
	  Int_t newest_nmix = 0;
	  MPI_Bcast(&newest_nmix, 1, MPI_INT, 0, MPI_COMM_WORLD);
	  mymodel.SetEstNMix(newest_nmix);

// 	  fprintf(stdout,"Process %d of %d is on %s: Estimating %d mixtures.\n",my_rank, np, processor_name, newest_nmix);
	}
      } 

      Int_t nfreeparam = mymodel.GetNfreeparam();

      MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(par, nfreeparam, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      mymodel.EvalLkhd(thisprob,thisgrad,thishess,par,flag,my_rank,np);

      MPI_Reduce(&thisprob, &totalprob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);      
      if (flag>=2) MPI_Reduce(thisgrad, totalgrad, nfreeparam, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if (flag==3) MPI_Reduce(thishess, totalhess, nfreeparam*(nfreeparam+1)/2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      MPI_Bcast(&finished, 1, MPI_INT, 0, MPI_COMM_WORLD);

    }
  }
  else {
    // Master Node:

    if ((mymodel.GetNSubSample()==0)&&(mymodel.GetNBootSample()==0)) {
      mymodel.Minimize(0);
    }
    else if (mymodel.GetNBootSample()>0) {
      for (UInt_t isample = mymodel.GetBootStrapStartSample(); isample < mymodel.GetNBootSample() ; isample++) {

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
    
    finished = 1;
    //    tag = 0;
    MPI_Bcast(&finished, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }

  //  cout << "Process " << my_rank << " Finished\n";
  MPI_Finalize();
  return 0;
    
}

void MinuitLkhdFcn(Int_t &npar, Double_t *gin, Double_t &f,
		   Double_t *par, Int_t iflag) {
  Int_t flag = iflag;
  if (flag !=2) flag = 1;
  LkhdFcn(npar,gin,f,par,iflag,NULL);
}


// void LkhdFcn(Int_t &npar, Double_t *gin, Double_t &f,
// 	     Double_t *par, Int_t iflag) {
void LkhdFcn(Int_t &npar, Double_t *gin, Double_t &f,
	     Double_t *par, Int_t iflag, Double_t *hess) {
  
  if ((mymodel.GetCounter())%100==0) {
    TString filename = "checkpoint.txt";
    filename.Prepend(mymodel.GetWorkingDir());

    FILE * pFile;
    pFile = fopen (filename.Data(),"w");
    for (Int_t ipar=0 ; ipar< npar ; ipar++)
      {
	fprintf (pFile, "%5d %12.8f %12.8f\n",ipar,par[ipar],0.0);
      }
    fclose (pFile);
  }
  
  f = 0;

  if (mymodel.GetInitializing()==1) {
    mymodel.EvalLkhd(f,gin,hess,par,iflag,0,1);
    return;
  }
  else {

    Int_t np;
    Int_t finished = 0;
    int    namelen;
    char   processor_name[MPI_MAX_PROCESSOR_NAME];
    
    //    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&np);
    MPI_Get_processor_name(processor_name,&namelen);

    MPI_Bcast(&finished, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    Int_t newflag = mymodel.GetNewFlag(); 
    MPI_Bcast(&newflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (newflag !=0 ) {
      
      //       Int_t newskipobs = mymodel.GetNewSkipObs();      
      //       MPI_Bcast(&newskipobs, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
      if (mymodel.GetNewSkipObs()==1) {
	vector<UInt_t> thisskipobs = mymodel.GetSkipObs();
	UInt_t nobs = mymodel.GetNobs();
	Int_t * skipobs = new Int_t[nobs];
	Int_t firstskip = -1;
	for (UInt_t iobs = 0; iobs < nobs ; iobs++) {
	  skipobs[iobs] = thisskipobs[iobs];
	  if ((firstskip==-1)&&(skipobs[iobs]==1)) firstskip = iobs;
	}
	MPI_Bcast(skipobs, nobs, MPI_INT, 0, MPI_COMM_WORLD);
	
	Int_t nskip = 0;
	for (UInt_t iobs = 0; iobs < nobs ; iobs++) nskip += skipobs[iobs];
	//      fprintf(stdout,"Process %d of %d: skipping %d, first=%d\n", 0, np, nskip,firstskip);
	mymodel.ResetSkipObs();
	delete [] skipobs;
      }
      
      
      //       Int_t newbootstrap = mymodel.GetNewBootStrap();      
      //       MPI_Bcast(&newbootstrap, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
      if (mymodel.GetNewBootStrap()==1) {
	vector<UInt_t> thisbootobs = mymodel.GetBootStrapObs();
	UInt_t nobs = mymodel.GetNobs();
	Int_t * bootstrapobs = new Int_t[nobs];
	for (UInt_t iobs = 0; iobs < nobs ; iobs++) bootstrapobs[iobs] = thisbootobs[iobs];
	MPI_Bcast(bootstrapobs, nobs, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (mymodel.GetPrintLevel()>0) fprintf(stdout,"Process %d of %d is on %s: BootStrap Set,first obs=%d, last obs=%d\n",
					       0, np, processor_name, bootstrapobs[0],bootstrapobs[nobs-1]);
	mymodel.ResetBootStrap();
	delete [] bootstrapobs;
      }
      
//       Int_t newfixpar = mymodel.GetNewFixPar();      
//       MPI_Bcast(&newfixpar, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
      if (mymodel.GetNewFixPar()==1) {
	vector<UInt_t> thisfixpar = mymodel.GetFixPar();
	
	UInt_t ntotparam = mymodel.GetNparam();
	Int_t * fixpar = new Int_t[ntotparam];
	for (UInt_t ipar = 0; ipar < ntotparam ; ipar++) {
	  fixpar[ipar] = thisfixpar[ipar];
	  //	  fprintf(stdout,"Param %8d Fixed to %d\n",ipar,fixpar[ipar]);
	}
	MPI_Bcast(fixpar, ntotparam, MPI_INT, 0, MPI_COMM_WORLD);
	
	Int_t nfreeparam = mymodel.GetNfreeparam();
	if (mymodel.GetPrintLevel()>0) fprintf(stdout,"Process %d of %d is on %s: has %d free parameters\n",0, np, processor_name, nfreeparam);
	
	mymodel.ResetNewFixPar();
	delete [] fixpar;
      }
      
      
//       Int_t newmodignore = mymodel.GetNewIgnoreMod();      
//       MPI_Bcast(&newmodignore, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
      if (mymodel.GetNewIgnoreMod()==1) {
	vector<UInt_t> thisignoremod = mymodel.GetIgnoreMod();
	
	UInt_t nmod = mymodel.GetNmodels();
	Int_t * ignoremod = new Int_t[nmod];
	for (UInt_t imod = 0; imod < nmod ; imod++) {
	  ignoremod[imod] = thisignoremod[imod];
	  //	  fprintf(stdout,"Model %8d Fixed to %d\n",imod,ignoremod[imod]);
	}
	MPI_Bcast(ignoremod, nmod, MPI_INT, 0, MPI_COMM_WORLD);
	
	Int_t nignoremod = mymodel.GetNignoremod();
	if (mymodel.GetPrintLevel()>0) fprintf(stdout,"Process %d of %d is on %s: ignoring %d models\n",0, np, processor_name, nignoremod);
	
	mymodel.ResetNewIgnoreMod();
      }
      
      if (mymodel.GetNewEstNMix()==1) {
	
	Int_t newest_nmix = mymodel.GetEstNMix();
	MPI_Bcast(&newest_nmix, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (mymodel.GetPrintLevel()>0) fprintf(stdout,"Process %d of %d is on %s: Estimating %d mixtures.\n",0, np, processor_name, newest_nmix);
	
	mymodel.ResetNewEstNMix();
      }
      
    }
    
    //send flag
    MPI_Bcast(&iflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // send free parameters
    MPI_Bcast(par, npar, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
    Double_t thisprob;
    Double_t * thisgrad = new Double_t[npar];
    Double_t * thishess = new Double_t[npar*(npar+1)/2];
    
    mymodel.EvalLkhd(thisprob,thisgrad,thishess,par,iflag,0,np);
    
    MPI_Reduce(&thisprob, &f, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (iflag>=2) MPI_Reduce(thisgrad, gin, npar, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (iflag==3) MPI_Reduce(thishess, hess, npar*(npar+1)/2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if ((mymodel.GetCounter())%100==0) printf("Iteration %5d Found logLkhd: %f \n",mymodel.GetCounter(),f);
    if ((mymodel.GetCounter())%100000==0) {
      printf("Finished Iteration %5d : current parameters are \n",mymodel.GetCounter());
      mymodel.PrintParam();
    }
    
    delete [] thisgrad;
    delete [] thishess;
  }
  return;
}
