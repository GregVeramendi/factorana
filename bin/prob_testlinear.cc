#ifndef __CINT__

#include "Riostream.h"
#include "Rtypes.h"
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <vector>
#include "TString.h"

#include "TMinLkhd.hh"

#endif


TMinLkhd mymodel("testlinear analysis","testlinear",0,0,1,8,1);
void init(int mpirank) {
  // do not touch this line:
  mymodel.SetMPRank(mpirank);
  
  // ***********************************************************
  // Editing starts from here:
  // ****************************

  //  mymodel.SetWorkingDir("/home/mendi/data/factorana/");
  mymodel.SetPrintLevel(2);
  //   mymodel.SetSimNobs(10000);
  
  // Defines how initial parameters are set for measurement system estimation
  // When running bootstrap samples, fullmodel_par.txt will be used instead of init_par.txt
  // mymodel.SetInitLoadingMultiple(0.01);
  //  mymodel.SetInitFixedLoading();
  // mymodel.SetInitBetaOLS();

  //Define dataset
  mymodel.SetData(6,(char *)"test/testdata_linear.txt",(char *) "test/testdata_linear_varlist.txt");

  // How Measurement system is estimated
  //  mymodel.EstSequenttialMixtures();
  //  mymodel.ConstrainFactorCorrelations();

  //For Bootstrapping
  //  mymodel.BootStrap(200);
  // mymodel.SetBootStrapStartSample(0); // first sample is 0

  //Double_t norm[1];
  //  Double_t norm[2];
//    norm[0] = -9999.0;
//   norm[1] = 1.0;

//   vector<TString> y_model;
//   y_model.push_back("y");
//   y_model.push_back("none");
//   y_model.push_back("x1");
//   y_model.push_back("x2");
//   y_model.push_back("x3");
//   y_model.push_back("const");

//   mymodel.AddModel("ymod","Test Linear model","linear",y_model);

  vector<TString> prob_model;
  prob_model.push_back("D");
  prob_model.push_back("none");
  prob_model.push_back("x1");
  prob_model.push_back("x2");
  prob_model.push_back("x3");
  prob_model.push_back("const");

  mymodel.AddModel("pmod","Test probit model","probit",prob_model);

      // mymodel.AllModel_Detailsim();
      // mymodel.SimIncludeData();

  if ((mpirank==0)&&(mymodel.GetPrintLevel()==1)) mymodel.PrintModels();
}
