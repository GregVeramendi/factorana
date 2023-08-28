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


TMinLkhd mymodel("testlinear analysis","testlinear",3,0,1,16,1);
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
  mymodel.SetData(16,(char *)"test/testdata_meassys_f3.txt",(char *) "test/testdata_meassys_f3_varlist.txt");

  // How Measurement system is estimated
  //  mymodel.EstSequenttialMixtures();
  //  mymodel.ConstrainFactorCorrelations();

  //For Bootstrapping
  //  mymodel.BootStrap(200);
  // mymodel.SetBootStrapStartSample(0); // first sample is 0

  // Double_t norm[1];
  Double_t norm[3];
  norm[0] = 1.0;
  norm[1] = 0.0;
  norm[2] = 0.0;

  vector<TString> M1_model;
  M1_model.push_back("M1");
  M1_model.push_back("none");
  M1_model.push_back("x1");
  M1_model.push_back("x2");
  M1_model.push_back("x3");
  //  M1_model.push_back("f");
  M1_model.push_back("const");
  mymodel.AddModel("m1","M1 model","linear",M1_model,norm);

  norm[0] = -9999.0;

  vector<TString> M2_model;
  M2_model.push_back("M2");
  M2_model.push_back("none");
  M2_model.push_back("x1");
  M2_model.push_back("x2");
  M2_model.push_back("x3");
  //  M2_model.push_back("f");
  M2_model.push_back("const");
  mymodel.AddModel("m2","M2 model","linear",M2_model,norm);

  vector<TString> M3_model;
  M3_model.push_back("M3");
  M3_model.push_back("none");
  M3_model.push_back("x1");
  M3_model.push_back("x2");
  M3_model.push_back("x3");
  //  M3_model.push_back("f");
  M3_model.push_back("const");
  mymodel.AddModel("m3","M3 model","linear",M3_model,norm);


  //  norm[0] = -9999.0;
  norm[0] = 0.0;
  norm[1] = 1.0;
  norm[2] = 0.0;

  vector<TString> M4_model;
  M4_model.push_back("M4");
  M4_model.push_back("none");
  M4_model.push_back("x1");
  M4_model.push_back("x2");
  M4_model.push_back("x3");
  //  M4_model.push_back("f");
  M4_model.push_back("const");
  mymodel.AddModel("m4","M4 model","linear",M4_model,norm);

  //  Int_t ifirstgpamodel = mymodel.GetNmodels() - 1;

  norm[1] = -9999.0;

  vector<TString> M5_model;
  M5_model.push_back("M5");
  M5_model.push_back("none");
  M5_model.push_back("x1");
  M5_model.push_back("x2");
  M5_model.push_back("x3");
  //  M5_model.push_back("f");
  M5_model.push_back("const");
  mymodel.AddModel("m5","M5 model","linear",M5_model,norm);
  //  mymodel.ConstrainLastFactorLoadingToModel(ifirstgpamodel,1);

  vector<TString> M6_model;
  M6_model.push_back("M6");
  M6_model.push_back("none");
  M6_model.push_back("x1");
  M6_model.push_back("x2");
  M6_model.push_back("x3");
  //  M6_model.push_back("f");
  M6_model.push_back("const");
  mymodel.AddModel("m6","M6 model","linear",M6_model,norm);
  //  mymodel.ConstrainLastFactorLoadingToModel(ifirstgpamodel,1);


  //  norm[0] = -9999.0;
  norm[0] = 0.0;
  norm[1] = 0.0;
  norm[2] = 1.0;

  vector<TString> M7_model;
  M7_model.push_back("M7");
  M7_model.push_back("none");
  M7_model.push_back("x1");
  M7_model.push_back("x2");
  M7_model.push_back("x3");
  //  M7_model.push_back("f");
  M7_model.push_back("const");
  mymodel.AddModel("m7","M7 model","linear",M7_model,norm);

  //  Int_t ifirstgpamodel = mymodel.GetNmodels() - 1;

  norm[2] = -9999.0;

  vector<TString> M8_model;
  M8_model.push_back("M8");
  M8_model.push_back("none");
  M8_model.push_back("x1");
  M8_model.push_back("x2");
  M8_model.push_back("x3");
  //  M8_model.push_back("f");
  M8_model.push_back("const");
  mymodel.AddModel("m8","M8 model","linear",M8_model,norm);
  //  mymodel.ConstrainLastFactorLoadingToModel(ifirstgpamodel,1);

  vector<TString> M9_model;
  M9_model.push_back("M9");
  M9_model.push_back("none");
  M9_model.push_back("x1");
  M9_model.push_back("x2");
  M9_model.push_back("x3");
  //  M9_model.push_back("f");
  M9_model.push_back("const");
  mymodel.AddModel("m9","M9 model","linear",M9_model,norm);
  //  mymodel.ConstrainLastFactorLoadingToModel(ifirstgpamodel,1);

      // mymodel.AllModel_Detailsim();
      // mymodel.SimIncludeData();

  if ((mpirank==0)&&(mymodel.GetPrintLevel()==1)) mymodel.PrintModels();
}
