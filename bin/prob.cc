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

//  meaning of numbers:
// number of factors, correlated bivariate normal, number of mixtures, number of integration points, analysis stage
// stages 1:Estimate measurement system 2: Estimate outcomes 3: Estimate first measurement system and then outcomes
TMinLkhd mymodel("testlinear analysis","testlinear",3,0,1,16,1);
void init(int mpirank) {
  // do not touch this line:
  mymodel.SetMPRank(mpirank);
  
  // ***********************************************************
  // Editing starts from here:
  // ****************************

  //   mymodel.SetWorkingDir("/home/rlanderso/factorana/temp/");
  mymodel.SetPrintLevel(1);
 // mymodel.SetSimNobs(10000);
  //  mymodel.CalcHessStdErr(0);
  
  // Defines how initial parameters are set for measurement system estimation
  // When running bootstrap samples, fullmodel_par.txt will be used instead of init_par.txt
  // mymodel.SetInitLoadingMultiple(0.01);
  //  mymodel.SetInitFixedLoading();
  mymodel.SetInitBetaOLS();

//////////////////////////
// DEFINE DATA SET AND NAME OF VARIABLE LIST
//////////////////////////

  //Define dataset
  mymodel.SetData(73,(char *)"/home/rlanderso/factorana/data_2.raw",(char *) "/home/rlanderso/factorana/varlist_2.txt");

  // How Measurement system is estimated
  //  mymodel.EstSequenttialMixtures();
  //  mymodel.ConstrainFactorCorrelations();

  //For Bootstrapping
 // mymodel.BootStrap(500);
 //   mymodel.SetBootStrapStartSample(0); // first sample is 0

  vector<TString> bkgnd;

  bkgnd.push_back("mom_smoke_preg");
  bkgnd.push_back("m_faminc79_mi");  
  bkgnd.push_back("m_broken");
  bkgnd.push_back("m_afqt_rawcorr");
  bkgnd.push_back("magebir1");
  bkgnd.push_back("magebir3");
  bkgnd.push_back("magebir4");  
  bkgnd.push_back("black");  
  bkgnd.push_back("hispanic");
  bkgnd.push_back("boy"); 
  bkgnd.push_back("impute_m_faminc79_m");  
  bkgnd.push_back("impute_m_broken");
  bkgnd.push_back("impute_m_afqt_rawcorr");
  bkgnd.push_back("constant"); 

  vector<TString> bkgnd2;
  bkgnd2.push_back("constant");

  vector<TString> bkgnd3;
  bkgnd3.push_back("boy"); 
  bkgnd3.push_back("constant");

  vector<TString> bkgnd4;
  bkgnd4.push_back("boy"); 
  bkgnd4.push_back("black");  
  bkgnd4.push_back("hispanic");
  bkgnd4.push_back("constant");


/////////////////////////////
// DEFINE NO OF NORMALIZATIONS. 1 factor=1, 2 factors=2, etc.... CHANGE LOADINGS TO -9999 IF FACTORS ARE TO BE INDEPENDENT
/////////////////////////////

  //Double_t norm[1];
 // Double_t norm[2];
 Double_t norm[3];
  norm[0] = -9999;
  norm[1] = 0.0;
  norm[2] = 0.0;


////////////////////////////
// 1  DEFINE TEST SCORE, 2  DEFINE INDICATOR OF MISSING, 3 DEFINE X's + constant, 4 DEFINE ESTIMATION (linear or probit)
///////////////////////////

//1
  vector<TString> M1_model;
  M1_model.push_back("PIAT_comp");  //First string is dependent variable
//2 
  M1_model.push_back("ind_PIAT_comp"); //Second string is indicator for missing
//3 
  M1_model.push_back("mom_smoke_preg");  //Third and all subsequent strings are covariates
  M1_model.push_back("m_faminc79_mi");
  M1_model.push_back("m_broken");
  M1_model.push_back("m_afqt_rawcorr");
  M1_model.push_back("magebir1");
  M1_model.push_back("magebir3");
  M1_model.push_back("magebir4");
  M1_model.push_back("boy");
  M1_model.push_back("black");
  M1_model.push_back("hispanic");
  M1_model.push_back("impute_m_faminc79_m");
  M1_model.push_back("impute_m_broken");
  M1_model.push_back("impute_m_afqt_rawcorr");
  M1_model.push_back("PIAT_comp_age10");
  //  M1_model.push_back("f");
  M1_model.push_back("constant");
//4
 ///  mymodel.AddModel("m1","M1 model","probit",M1_model,norm);
  mymodel.AddModel("m1","M1 model","linear",M1_model,norm);
  mymodel.LastModel_Detailsim();

  vector<TString> M2_model;
  M2_model.push_back("PIAT_math");
  M2_model.push_back("ind_PIAT_math");
  M2_model.push_back("mom_smoke_preg");
  M2_model.push_back("m_faminc79_mi");
  M2_model.push_back("m_broken");
  M2_model.push_back("m_afqt_rawcorr");
  M2_model.push_back("magebir1");
  M2_model.push_back("magebir3");
  M2_model.push_back("magebir4");
  M2_model.push_back("boy");
  M2_model.push_back("black");
  M2_model.push_back("hispanic");
  M2_model.push_back("impute_m_faminc79_m");
  M2_model.push_back("impute_m_broken");
  M2_model.push_back("impute_m_afqt_rawcorr");
  M2_model.push_back("PIAT_math_age10");
  //  M2_model.push_back("f");
  M2_model.push_back("constant");
  mymodel.AddModel("m2","M2 model","linear",M2_model,norm);
  mymodel.LastModel_Detailsim();

  norm[0]=1.0;

  vector<TString> M3_model;
  M3_model.push_back("PIAT_recog");
  M3_model.push_back("ind_PIAT_recog");
  M3_model.push_back("mom_smoke_preg");
  M3_model.push_back("m_faminc79_mi");
  M3_model.push_back("m_broken");
  M3_model.push_back("m_afqt_rawcorr");
  M3_model.push_back("magebir1");
  M3_model.push_back("magebir3");
  M3_model.push_back("magebir4");
  M3_model.push_back("boy");
  M3_model.push_back("black");
  M3_model.push_back("hispanic");
  M3_model.push_back("impute_m_faminc79_m");
  M3_model.push_back("impute_m_broken");
  M3_model.push_back("impute_m_afqt_rawcorr");
  M3_model.push_back("PIAT_recog_age10");
  //  M3_model.push_back("f");
  M3_model.push_back("constant");
  mymodel.AddModel("m3","M3 model","linear",M3_model,norm);
  mymodel.LastModel_Detailsim();
  
  
  /////////////////////////////////////////////
  // FACTOR 2
  /////////////////////////////////////////////

  // 2 factor
  norm[0]=0.0;
  norm[1]=-1.0;

  // 3 factor
  norm[2]=0.0;
  
  vector<TString> M4_model;
  M4_model.push_back("ANTI");
  M4_model.push_back("ind_ANTI"); 
  M4_model.push_back("mom_smoke_preg"); 
  M4_model.push_back("m_faminc79_mi");  
  M4_model.push_back("m_broken");
  M4_model.push_back("m_afqt_rawcorr");
  M4_model.push_back("magebir1");
  M4_model.push_back("magebir3");
  M4_model.push_back("magebir4"); 
  M4_model.push_back("boy");
  M4_model.push_back("black");
  M4_model.push_back("hispanic");
  M4_model.push_back("impute_m_faminc79_m");  
  M4_model.push_back("impute_m_broken");
  M4_model.push_back("impute_m_afqt_rawcorr");
  M4_model.push_back("ANTI_age10");
  //  M4_model.push_back("f");
  M4_model.push_back("constant");
 

  //  mymodel.AddModel("m4","M4 model","probit",M4_model,norm);
  mymodel.AddModel("m4","M4 model","linear",M4_model,norm);
  mymodel.LastModel_Detailsim();

 // 2 factor
 //norm[1]=-9999;

  // 3 factor
  norm[1]=-9999;
  norm[2]=-1.0;

  vector<TString> M5_model;
  M5_model.push_back("ANX");
  M5_model.push_back("ind_ANX");
  M5_model.push_back("mom_smoke_preg");
  M5_model.push_back("m_faminc79_mi");
  M5_model.push_back("m_broken");
  M5_model.push_back("m_afqt_rawcorr");
  M5_model.push_back("magebir1");
  M5_model.push_back("magebir3");
  M5_model.push_back("magebir4"); 
  M5_model.push_back("boy");
  M5_model.push_back("black");
  M5_model.push_back("hispanic"); 
  M5_model.push_back("impute_m_faminc79_m");  
  M5_model.push_back("impute_m_broken");
  M5_model.push_back("impute_m_afqt_rawcorr");
  M5_model.push_back("ANX_age10");
  //  M5_model.push_back("f");
  M5_model.push_back("constant");
  mymodel.AddModel("m5","M5 model","linear",M5_model,norm);
  mymodel.LastModel_Detailsim();

  // 3 factor
  norm[2]=-9999;

  vector<TString> M6_model;
  M6_model.push_back("DEP");
  M6_model.push_back("ind_DEP");
  M6_model.push_back("mom_smoke_preg");
  M6_model.push_back("m_faminc79_mi");  
  M6_model.push_back("m_broken");
  M6_model.push_back("m_afqt_rawcorr");
  M6_model.push_back("magebir1");
  M6_model.push_back("magebir3");
  M6_model.push_back("magebir4");
  M6_model.push_back("boy");
  M6_model.push_back("black");
  M6_model.push_back("hispanic");
  M6_model.push_back("impute_m_faminc79_m");  
  M6_model.push_back("impute_m_broken");
  M6_model.push_back("impute_m_afqt_rawcorr");
  M6_model.push_back("DEP_age10"); 
  //  M6_model.push_back("f");
  M6_model.push_back("constant");
  mymodel.AddModel("m6","M6 model","linear",M6_model,norm);  
  mymodel.LastModel_Detailsim();

 // 3 factor
  norm[1]=-9999;
  norm[2]=0.0;

  vector<TString> M7_model;
  M7_model.push_back("HEAD");
  M7_model.push_back("ind_HEAD");
  M7_model.push_back("mom_smoke_preg");
  M7_model.push_back("m_faminc79_mi");  
  M7_model.push_back("m_broken");
  M7_model.push_back("m_afqt_rawcorr");
  M7_model.push_back("magebir1");
  M7_model.push_back("magebir3");
  M7_model.push_back("magebir4");
  M7_model.push_back("boy");
  M7_model.push_back("black");
  M7_model.push_back("hispanic");
  M7_model.push_back("impute_m_faminc79_m");  
  M7_model.push_back("impute_m_broken");
  M7_model.push_back("impute_m_afqt_rawcorr");
  M7_model.push_back("HEAD_age10"); 
  //  M7_model.push_back("f");
  M7_model.push_back("constant");
  mymodel.AddModel("m7","M7 model","linear",M7_model,norm); 
  mymodel.LastModel_Detailsim();

  // 3 factor
  norm[1]=-9999;

  vector<TString> M8_model;
  M8_model.push_back("HYPR");
  M8_model.push_back("ind_HYPR");
  M8_model.push_back("mom_smoke_preg");
  M8_model.push_back("m_faminc79_mi");  
  M8_model.push_back("m_broken");
  M8_model.push_back("m_afqt_rawcorr");
  M8_model.push_back("magebir1");
  M8_model.push_back("magebir3");
  M8_model.push_back("magebir4");
  M8_model.push_back("boy");
  M8_model.push_back("black");
  M8_model.push_back("hispanic");
  M8_model.push_back("impute_m_faminc79_m");  
  M8_model.push_back("impute_m_broken");
  M8_model.push_back("impute_m_afqt_rawcorr");
  M8_model.push_back("HYPR_age10"); 
  //  M8_model.push_back("f");
  M8_model.push_back("constant");
  mymodel.AddModel("m8","M8 model","linear",M8_model,norm);
  mymodel.LastModel_Detailsim();

  // 3 factor
  norm[1]=-9999;
  norm[2]=-9999;
  
  vector<TString> M9_model;
  M9_model.push_back("PEER");
  M9_model.push_back("ind_PEER");
  M9_model.push_back("mom_smoke_preg");
  M9_model.push_back("m_faminc79_mi");  
  M9_model.push_back("m_broken");
  M9_model.push_back("m_afqt_rawcorr");
  M9_model.push_back("magebir1");
  M9_model.push_back("magebir3");
  M9_model.push_back("magebir4");  
  M9_model.push_back("boy");
  M9_model.push_back("black");
  M9_model.push_back("hispanic"); 
  M9_model.push_back("impute_m_faminc79_m");  
  M9_model.push_back("impute_m_broken");
  M9_model.push_back("impute_m_afqt_rawcorr");
  M9_model.push_back("PEER_age10");
  //  M9_model.push_back("f");
  M9_model.push_back("constant");
  mymodel.AddModel("m9","M9model","linear",M9_model,norm); 
  mymodel.LastModel_Detailsim();

    mymodel.StartNewEstGroup();
    mymodel.StartNewPrintGroup();

   vector<TString> crime;
   crime.push_back("crime");
   crime.push_back("ind_crime");
   crime.insert(crime.end(), bkgnd.begin(), bkgnd.end());
   mymodel.AddModel("crime","Self-reported crime before age 18, probit","probit",crime);
  mymodel.LastModel_Detailsim();

//   mymodel.StartNewEstGroup();
  
//   vector<TString> crime2;
//   crime2.push_back("crime");
//   crime2.push_back("ind_crime");   
//   crime2.insert(crime2.end(), bkgnd2.begin(), bkgnd2.end());
//   mymodel.AddModel("crime","Self-reported crime before age 18, no covar","probit",crime2);

//   mymodel.StartNewEstGroup();
  
//   vector<TString> crime3;
//   crime3.push_back("crime");
//   crime3.push_back("ind_crime");   
//   crime3.insert(crime3.end(), bkgnd3.begin(), bkgnd3.end());
//   mymodel.AddModel("crime","Self-reported crime before age 18, boy","probit",crime3);

//   mymodel.StartNewEstGroup();
  
//   vector<TString> crime4;
//   crime4.push_back("crime");
//   crime4.push_back("ind_crime");
//   crime4.insert(crime4.end(), bkgnd4.begin(), bkgnd4.end());
//   mymodel.AddModel("crime","Self-reported crime before age 18, boy+black","probit",crime4);

//   mymodel.StartNewEstGroup();

   vector<TString> violence;
   violence.push_back("violence");
   violence.push_back("ind_violence");
   violence.insert(violence.end(), bkgnd.begin(), bkgnd.end());
//   mymodel.AddModel("violence","Self-reported violence before age 18","probit",violence);
//   mymodel.LastModel_Detailsim();

//   mymodel.StartNewEstGroup();

   vector<TString> property;
   property.push_back("property");
   property.push_back("ind_property");
   property.insert(property.end(), bkgnd.begin(), bkgnd.end());
//   mymodel.AddModel("property","Self-reported property before age 18","probit",property);
//   mymodel.LastModel_Detailsim();

//   mymodel.StartNewEstGroup();

   vector<TString> drugs;
   drugs.push_back("drugs");
   drugs.push_back("ind_drugs");
   drugs.insert(drugs.end(), bkgnd.begin(), bkgnd.end());
//   mymodel.AddModel("drugs","Self-reported drug related crime  before age 18","probit",drugs);
//   mymodel.LastModel_Detailsim();

      // mymodel.AllModel_Detailsim();
      // mymodel.SimIncludeData();

  if ((mpirank==0)&&(mymodel.GetPrintLevel()==1)) mymodel.PrintModels();

}
