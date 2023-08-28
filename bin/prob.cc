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

TMinLkhd mymodel("Major Choice analysis","MajorChoice",1,0,1,8,0);
void init(int mpirank) {
  // do not touch this line:
  mymodel.SetMPRank(mpirank);

// ============================================================
// ***********************************************************
// Editing starts from here:
// ************************************************************
// ============================================================

    mymodel.SetWorkingDir("/home/Projects/MajorChoice/estimation");
    mymodel.SetSimNobs(10000000);

    mymodel.SetPrintLevel(2);

    //    mymodel.SetCPULimit(3600);
    mymodel.SetStochasticDeriv(0.1);

    
// ********************************************************
// Some Controls to Help Convergece
// *******************************************************  
 // Defines how initial parameters are set for measurement system estimation
 // When running bootstrap samples, fullmodel_par.txt will be used instead of init_par.txt
    //    mymodel.SetInitLoadingMultiple(0.10);
    //    mymodel.SetInitFixedLoading();
    mymodel.SetInitBetaOLS();
    mymodel.SetInitLoadingMultiple(0.1);


//Define dataset*****************************************************
    mymodel.SetData((char *)"/home/Projects/MajorChoice/intdata/MajorChoice/majorchoice_wfac_v8.raw",(char *) "/home/Projects/MajorChoice/intdata/MajorChoice/varlist_majorchoice_wfac_v8.txt");
  mymodel.SetObsIndex("lopnr");

// How Measurement system is estimated *******************************
   // mymodel.EstSequenttialMixtures();

//For Bootstrapping***************************************************
//  mymodel.BootStrap(2,0);
//    mymodel.SetBootStrapStartSample(124); // first sample is 0

// Defining Vector of Factors *****************************************
  //Double_t norm[1];
  Double_t norm[3];
  
// ********************************************************************  
// Defining Lists of Parameters ***************************************  
// ********************************************************************  
  if (mpirank==0) cout << "Hello from code point 1" << endl; 
  // BACKGROUND ***
  vector<TString> bkgnd1990;
  bkgnd1990.push_back("college_mom");
  // bkgnd.push_back("college_mom_miss");
   bkgnd1990.push_back("hs_mom");
   bkgnd1990.push_back("hs_mom_miss");
   bkgnd1990.push_back("college_dad");
   bkgnd1990.push_back("hs_dad");
   bkgnd1990.push_back("hs_dad_miss");
   //   bkgnd1990.push_back("stockholm"); ***
   //   bkgnd1990.push_back("metro"); ***
  // bkgnd1990.push_back("college_mom1990");
  // // bkgnd1990.push_back("college_mom_miss1990");
  //  bkgnd1990.push_back("hs_mom1990");
  //  bkgnd1990.push_back("hs_mom_miss1990");
  //  bkgnd1990.push_back("stockholm1990");
  //  bkgnd1990.push_back("metro1990");
   // bkgnd1990.push_back("bvikt_lbv");
   // bkgnd1990.push_back("bvikt_lbv_miss");
   // bkgnd1990.push_back("malder");
   // bkgnd1990.push_back("malder_miss");
//    bkgnd1990.push_back("fam_inc73_mom");
//    bkgnd1990.push_back("fam_inc73_mom_squared");
//    //   bkgnd1990.push_back("fam_inc73_mom_miss");
//    // bkgnd1990.push_back("mfod_swe");
//    // bkgnd1990.push_back("mfod_swe_miss");
//    bkgnd1990.push_back("mmarried");
//    bkgnd1990.push_back("mmarried_miss");
//    bkgnd1990.push_back("childhealthy");
//    bkgnd1990.push_back("childhealthy_miss");
//    // bkgnd1990.push_back("APGAR1");
//    // bkgnd1990.push_back("APGAR1_miss");
//    //   bkgnd1990.push_back("b1973");
//    bkgnd1990.push_back("malder_lt20");
//    bkgnd1990.push_back("malder_20to24");
//    bkgnd1990.push_back("malder_30to34");
//    bkgnd1990.push_back("malder_gt34");
// // bkgnd1990.push_back("schoolave_advmath");
// // bkgnd1990.push_back("schoolave_adveng");
//    //   bkgnd1990.push_back("malder_miss");
//    bkgnd1990.push_back("b1974");
//    bkgnd1990.push_back("b1975");
//    //   bkgnd1990.push_back("pf");
//    //   bkgnd1990.push_back("pf_miss");
//    // bkgnd1990.push_back("lead");
//    // bkgnd1990.push_back("lead_miss");
//    //   bkgnd1990.push_back("constant");

//  vector<TString> bkgnd1990_nocons = bkgnd1990;

  bkgnd1990.push_back("constant");

  // vector<TString> bkgnd2010;
  // bkgnd2010.push_back("college_mom");
  // // bkgnd2010.push_back("college_mom_miss");
  // bkgnd2010.push_back("hs_mom");
  // bkgnd2010.push_back("hs_mom_miss");
  // bkgnd2010.push_back("stockholm");
  // bkgnd2010.push_back("metro");
  // bkgnd2010.push_back("constant");
 
  // vector<TString> enrollmaj_indic;
  // //  enrollmaj_indic.push_back("ind_enroll_maj1");
  // enrollmaj_indic.push_back("Denroll2");
  // enrollmaj_indic.push_back("Denroll3");
  // enrollmaj_indic.push_back("Denroll4");
  // enrollmaj_indic.push_back("Denroll5");
  // enrollmaj_indic.push_back("Denroll6");
  // enrollmaj_indic.push_back("Denroll7");
  // enrollmaj_indic.push_back("Denroll8");
  // enrollmaj_indic.push_back("Denroll9");
  // enrollmaj_indic.push_back("Denroll10");
  // enrollmaj_indic.push_back("Denroll11");
  // enrollmaj_indic.push_back("Denroll12");
  // enrollmaj_indic.push_back("Denroll13");
  

  if (mpirank==0) cout << "Hello from code point 2" << endl;

  // *****************************************************
  // 4 education choices
  // ****************************************************
  
if (mpirank==0) cout << "Hello from code point 3" << endl;

  // ************************************************************************
  // GPA 
  // ************************************************************************

   if (mpirank==0) cout << "Hello from code point 5.2" << endl;
   
// *****************************************************************
// Measurement system for cognitive ability
// ******************************************************************

// if (mpirank==0) cout << "Hello from code point 5a" << endl;

//  mymodel.StartNewPrintGroup();

   norm[0] = -9999.0; 
   norm[1] = 0.0;
   norm[2] =0.788518272521786;

   /*
   vector<TString> advmath;
   advmath.push_back("advmath");
   advmath.push_back("ind_advmath");
   advmath.insert(advmath.end(), bkgnd1990.begin(), bkgnd1990.end());
   //   advmath.push_back("advmathIV");
   mymodel.AddModel("advmath","9th Grade Advanced Math","probit",advmath,norm);
   mymodel.LastModel_Detailsim();
   uint advmathModelN = mymodel.GetNmodels()-1;
   norm[0] = -9999.0; 
   norm[1] = 0.0;
   norm[2] = -9999.0; 
   
   vector<TString> adveng;
   adveng.push_back("adveng");
   adveng.push_back("ind_adveng");
   adveng.insert(adveng.end(), bkgnd1990.begin(), bkgnd1990.end());
   //   adveng.push_back("advengIV");
   mymodel.AddModel("adveng","9th Grade Advanced English","probit",adveng,norm);
   mymodel.LastModel_Detailsim();
   uint advengModelN = mymodel.GetNmodels()-1;
   
   vector<TString> seng9;
   seng9.push_back("seng9");
   seng9.push_back("ind_seng9");
   seng9.insert(seng9.end(), bkgnd1990.begin(), bkgnd1990.end());
   seng9.push_back("adveng");
   //   seng9.push_back("adveng_miss");
   mymodel.AddModel("seng9","9th Grade English Grade","linear",seng9,norm);
   mymodel.LastModel_SetEndogenousReg(advengModelN);
   uint sengModelN = mymodel.GetNmodels()-1;
   
   vector<TString> smath9;
   smath9.push_back("smath9");
   smath9.push_back("ind_smath9");
   smath9.insert(smath9.end(), bkgnd1990.begin(), bkgnd1990.end());
   smath9.push_back("advmath");
   //   smath9.push_back("advmath_miss");
   mymodel.AddModel("smath9","9th Grade Math Grade","linear",smath9,norm);
   mymodel.LastModel_SetEndogenousReg(advmathModelN);
   uint smathModelN = mymodel.GetNmodels()-1;
   
   vector<TString> ssports9;
   ssports9.push_back("ssports9");
   ssports9.push_back("ind_ssports9");
   ssports9.insert(ssports9.end(), bkgnd1990.begin(), bkgnd1990.end());
   mymodel.AddModel("ssports9","9th Grade Sports Grade","linear",ssports9);
   uint ssportsModelN = mymodel.GetNmodels()-1;
   
   vector<TString> sgpa9;
   sgpa9.push_back("sgpa9");
   sgpa9.push_back("ind_sgpa9");
   sgpa9.insert(sgpa9.end(), bkgnd1990.begin(), bkgnd1990.end());
   sgpa9.push_back("adveng");
   //   sgpa9.push_back("adveng_miss");
   sgpa9.push_back("advmath");
   //   sgpa9.push_back("advmath_miss");
   sgpa9.push_back("smath9");
   //   sgpa9.push_back("smath9_miss");
   sgpa9.push_back("seng9");
   //   sgpa9.push_back("seng9_miss");
   sgpa9.push_back("ssports9");
   sgpa9.push_back("ssports9_miss");
   mymodel.AddModel("sgpa9","9th Grade GPA","linear",sgpa9);
   mymodel.LastModel_SetEndogenousReg(advengModelN);
   mymodel.LastModel_SetEndogenousReg(advmathModelN);
   mymodel.LastModel_SetEndogenousReg(smathModelN);
   mymodel.LastModel_SetEndogenousReg(sengModelN);
   mymodel.LastModel_SetEndogenousReg(ssportsModelN);
   mymodel.LastModel_Detailsim();

   mymodel.StartNewPrintGroup();
   
   vector<TString> HS_track;
   HS_track.push_back("HS_track");
   HS_track.push_back("none");
   HS_track.insert(HS_track.end(), bkgnd1990.begin(), bkgnd1990.end());
   HS_track.push_back("adveng");
   //   HS_track.push_back("adveng_miss");
   HS_track.push_back("advmath");
   //   HS_track.push_back("advmath_miss");
   // HS_track.push_back("schoolave_hstrk2");
   // HS_track.push_back("hstrk2IV");
   // HS_track.push_back("schoolave_hstrk3");
   //   HS_track.push_back("hstrk3IV");
   mymodel.AddModel("HS_track","High School Track","logit",HS_track,NULL,3);
   mymodel.LastModel_SetEndogenousReg(advengModelN);
   mymodel.LastModel_SetEndogenousReg(advmathModelN);
//   mymodel.LastModel_FixParamValue("schoolave_hstrk3", 0.0,2);
//   mymodel.LastModel_FixParamValue("hstrk3IV", 0.0,2);
//   mymodel.LastModel_FixParamValue("schoolave_hstrk2", 0.0,3);
//   mymodel.LastModel_FixParamValue("hstrk2IV", 0.0,3);
   mymodel.LastModel_Detailsim();
   uint HSTrackModelN = mymodel.GetNmodels()-1;
   
   vector<TString> ssports10;
   ssports10.push_back("ssports10");
   ssports10.push_back("ind_ssports10");
   ssports10.insert(ssports10.end(), bkgnd1990.begin(), bkgnd1990.end());
   ssports10.push_back("HS_track2");
   ssports10.push_back("HS_track3");
   mymodel.AddModel("ssports10","10th Grade Sports Grade","linear",ssports10);
   mymodel.LastModel_SetEndogenousReg(HSTrackModelN);
   mymodel.LastModel_Detailsim();
   
   vector<TString> sgpahs;
   sgpahs.push_back("sgpahs");
   sgpahs.push_back("ind_sgpahs");
   sgpahs.insert(sgpahs.end(), bkgnd1990.begin(), bkgnd1990.end());
   sgpahs.push_back("adveng");
   //   sgpahs.push_back("adveng_miss");
   sgpahs.push_back("advmath");
   //   sgpahs.push_back("advmath_miss");
   sgpahs.push_back("HS_track2");
   sgpahs.push_back("HS_track3");
   mymodel.AddModel("sgpahs","HS GPA","linear",sgpahs);
   mymodel.LastModel_SetEndogenousReg(advengModelN);
   mymodel.LastModel_SetEndogenousReg(advmathModelN);
   mymodel.LastModel_SetEndogenousReg(HSTrackModelN);
   mymodel.LastModel_Detailsim();   
   uint HSgpaModelN = mymodel.GetNmodels()-1;

   mymodel.StartNewPrintGroup();

   */
   
   norm[0] = 0.674820084556316;
   norm[1] = 0.0;
   norm[2] = 0.0;

  // //JE's model:
  // norm[1] = -9999.0;

   if (mpirank==0) cout << "Hello from code point 5.2" << endl;
 
   vector<TString> cog1;
   cog1.push_back("cog1");
   cog1.push_back("ind_cog1");
   cog1.insert(cog1.end(), bkgnd1990.begin(), bkgnd1990.end());
   // cog1.push_back("adveng");
   // //   cog1.push_back("adveng_miss");
   // cog1.push_back("advmath");
   // //   cog1.push_back("advmath_miss");
   // cog1.push_back("HS_track2");
   // cog1.push_back("HS_track3");
   mymodel.AddModel("cog1","Cognitive measure 1","linear",cog1,norm);
   // mymodel.LastModel_SetEndogenousReg(advengModelN);
   // mymodel.LastModel_SetEndogenousReg(advmathModelN);
   // mymodel.LastModel_SetEndogenousReg(HSTrackModelN);
   
   norm[0] = -9999.0;

   vector<TString> cog2;
   cog2.push_back("cog2");
   cog2.push_back("ind_cog2");
   cog2.insert(cog2.end(), bkgnd1990.begin(), bkgnd1990.end());
   // cog2.push_back("adveng");
   // //   cog2.push_back("adveng_miss");
   // cog2.push_back("advmath");
   // //   cog2.push_back("advmath_miss");
   // cog2.push_back("HS_track2");
   // cog2.push_back("HS_track3");
   mymodel.AddModel("cog2","Cognitive measure 2","linear",cog2,norm);
   // mymodel.LastModel_SetEndogenousReg(advengModelN);
   // mymodel.LastModel_SetEndogenousReg(advmathModelN);
   // mymodel.LastModel_SetEndogenousReg(HSTrackModelN);
   
   // //JE's model:
   // norm[1] = 0.0;

   vector<TString> cog3;
   cog3.push_back("cog3");
   cog3.push_back("ind_cog3");
   cog3.insert(cog3.end(), bkgnd1990.begin(), bkgnd1990.end());
   // cog3.push_back("adveng");
   // //   cog3.push_back("adveng_miss");
   // cog3.push_back("advmath");
   // //   cog3.push_back("advmath_miss");
   // cog3.push_back("HS_track2");
   // cog3.push_back("HS_track3");
   mymodel.AddModel("cog3","Cognitive measure 3","linear",cog3,norm);
   // mymodel.LastModel_SetEndogenousReg(advengModelN);
   // mymodel.LastModel_SetEndogenousReg(advmathModelN);
   // mymodel.LastModel_SetEndogenousReg(HSTrackModelN);
   
   // //JE's model:
   // norm[1] = -9999.0; 

   vector<TString> cog4;
   cog4.push_back("cog4");
   cog4.push_back("ind_cog4");
   cog4.insert(cog4.end(), bkgnd1990.begin(), bkgnd1990.end());
   // cog4.push_back("adveng");
   // //   cog4.push_back("adveng_miss");
   // cog4.push_back("advmath");
   // //   cog4.push_back("advmath_miss");
   // cog4.push_back("HS_track2");
   // cog4.push_back("HS_track3");
   mymodel.AddModel("cog4","Cognitive measure 4","linear",cog4,norm);
   // mymodel.LastModel_SetEndogenousReg(advengModelN);
   // mymodel.LastModel_SetEndogenousReg(advmathModelN);
   // mymodel.LastModel_SetEndogenousReg(HSTrackModelN);
   
   /*
   mymodel.StartNewPrintGroup();

   norm[0] = -9999.0;
   norm[1] = 0.0;
   norm[2] = 0.0;
   vector<TString> leadmiss;
   leadmiss.push_back("miss_lead");
   leadmiss.push_back("none");
   leadmiss.insert(leadmiss.end(), bkgnd1990.begin(), bkgnd1990.end());
   leadmiss.push_back("adveng");
   //   leadmiss.push_back("adveng_miss");
   leadmiss.push_back("advmath");
   //   leadmiss.push_back("advmath_miss");
   leadmiss.push_back("HS_track2");
   leadmiss.push_back("HS_track3");
   mymodel.AddModel("misslead","lead missing","probit",leadmiss,norm);
   mymodel.LastModel_SetEndogenousReg(advengModelN);
   mymodel.LastModel_SetEndogenousReg(advmathModelN);
   mymodel.LastModel_SetEndogenousReg(HSTrackModelN);
   
  // // Baseline model 
     norm[0] = -9999.0; 
     norm[1] = 0.801113447673862;
     norm[2] = 0.0;

   vector<TString> stresstol;
   stresstol.push_back("pf");
   stresstol.push_back("ind_pf");
   stresstol.insert(stresstol.end(), bkgnd1990.begin(), bkgnd1990.end());
   stresstol.push_back("adveng");
   //   stresstol.push_back("adveng_miss");
   stresstol.push_back("advmath");
   //   stresstol.push_back("advmath_miss");
   stresstol.push_back("HS_track2");
   stresstol.push_back("HS_track3");
   mymodel.AddModel("stresstol","Emotional Stability","linear",stresstol,norm);
   mymodel.LastModel_SetEndogenousReg(advengModelN);
   mymodel.LastModel_SetEndogenousReg(advmathModelN);
   mymodel.LastModel_SetEndogenousReg(HSTrackModelN);
   
     norm[0] = -9999.0; 
     norm[1] = -9999.0;
     norm[2] = 0.0;

   vector<TString> lead;
   lead.push_back("lead");
   lead.push_back("ind_lead");
   lead.insert(lead.end(), bkgnd1990.begin(), bkgnd1990.end());
   lead.push_back("adveng");
   //   lead.push_back("adveng_miss");
   lead.push_back("advmath");
   //   lead.push_back("advmath_miss");
   lead.push_back("HS_track2");
   lead.push_back("HS_track3");
   mymodel.AddModel("lead","Leadership","linear",lead,norm);
   mymodel.LastModel_SetEndogenousReg(advengModelN);
   mymodel.LastModel_SetEndogenousReg(advmathModelN);
   mymodel.LastModel_SetEndogenousReg(HSTrackModelN);


//   // ***********************************************************************
//   // *****************OUTCOMES**********************************************
//   // ***********************************************************************

   mymodel.StartNewPrintGroup();

   // Comment this out for estimating income
   mymodel.StartNewEstGroup();

   
  // vector<TString> educD1;
  // educD1.push_back("DEnroll");
  // educD1.push_back("ind_DEnroll");
  // educD1.insert(educD1.end(), bkgnd1990.begin(), bkgnd1990.end());
  // mymodel.AddModel("educDEnroll","Enroll College","probit",educD1,NULL);
  // mymodel.LastModel_Detailsim();

  //  mymodel.StartNewPrintGroup();

  vector<TString> educD2;
  educD2.push_back("Denroll");
  educD2.push_back("ind_Denroll");
  educD2.insert(educD2.end(), bkgnd1990.begin(), bkgnd1990.end());
   educD2.push_back("adveng");
   //   educD2.push_back("adveng_miss");
   educD2.push_back("advmath");
   //   educD2.push_back("advmath_miss");
  educD2.push_back("HS_track2");
  educD2.push_back("HS_track3");
  educD2.push_back("sgpahs");
  mymodel.AddModel("educDenroll","Major","logit",educD2,NULL,13);  
  mymodel.LastModel_Detailsim();
  mymodel.LastModel_SetEndogenousReg(advengModelN);
  mymodel.LastModel_SetEndogenousReg(advmathModelN);
  mymodel.LastModel_SetEndogenousReg(HSTrackModelN);
  mymodel.LastModel_SetEndogenousReg(HSgpaModelN);

  mymodel.StartNewPrintGroup();
  
  uint MajorChoiceModelN = mymodel.GetNmodels()-1;

  //  mymodel.StartNewPrintGroup();


  string outcomevar, indvar, name, title;
  string major[12] = {"1","2","3","4","5","6","7","8","9","10","11","12"};

  // for (int imajor = 0; imajor < 12 ; imajor++) {
    
  //   if (imajor==6) mymodel.StartNewPrintGroup();
  //   // *** Estimate each model conditional on major:
    
  //   outcomevar = "Dswitch_" + major[imajor];
  //   indvar = "ind_Dswitch_" + major[imajor];
  //   name = "educDSwitchfromMaj" + major[imajor];
  //   title = "Switch from Major " + major[imajor];
    
  //   // cout << "Outcome: " << outcomevar.c_str() << endl;
  //   // cout << "Indicator: " << indvar.c_str() << endl;
  //   // cout << "Name: " << name.c_str() << endl;
  //   // cout << "Title: " << title.c_str() << endl << endl;;
    
  //   vector<TString> modelvars;
  //   modelvars.push_back(outcomevar.c_str());
  //   modelvars.push_back(indvar.c_str());
  //   modelvars.insert(modelvars.end(), bkgnd1990.begin(), bkgnd1990.end());
  //   mymodel.AddModel(name.c_str(),title.c_str(),"probit", modelvars);
  //   mymodel.LastModel_Detailsim();
  // }

  //  mymodel.StartNewPrintGroup();

  vector<TString> educD3;
  educD3.push_back("DSwitchMajor");
  educD3.push_back("ind_DSwitchMajor");
  educD3.insert(educD3.end(), bkgnd1990.begin(), bkgnd1990.end());
  educD3.push_back("adveng");
  //  educD3.push_back("adveng_miss");
  educD3.push_back("advmath");
  //  educD3.push_back("advmath_miss");
  educD3.push_back("HS_track2");
  educD3.push_back("HS_track3");
  educD3.push_back("sgpahs");
  educD3.insert(educD3.end(), enrollmaj_indic.begin(), enrollmaj_indic.end()); 
  mymodel.AddModel("educDSwitchtoMajor","Switch to Major","logit",educD3,NULL,12);
  //  mymodel.AddModel("educD2","Education Model D2 (VOC HSG)","orderedprobit",educD2,NULL,3);
  mymodel.LastModel_Detailsim();
  mymodel.LastModel_SetEndogenousReg(advengModelN);
  mymodel.LastModel_SetEndogenousReg(advmathModelN);
  mymodel.LastModel_SetEndogenousReg(HSTrackModelN);
  mymodel.LastModel_SetEndogenousReg(HSgpaModelN);
  mymodel.LastModel_SetEndogenousReg(MajorChoiceModelN);


  // Set omitted categories to be transitions with less than 10 observations
  // Choice 2:
  mymodel.LastModel_FixParamValue("Denroll7", 0.0,2);
  mymodel.LastModel_FixParamValue("Denroll11", 0.0,2);
  mymodel.LastModel_FixParamValue("Denroll13", 0.0,2);

  // Choice 3:
  mymodel.LastModel_FixParamValue("Denroll5", 0.0,3);
  mymodel.LastModel_FixParamValue("Denroll7", 0.0,3);
  mymodel.LastModel_FixParamValue("Denroll11", 0.0,3);
  mymodel.LastModel_FixParamValue("Denroll13", 0.0,3);
  
  // // Choice 4:
  mymodel.LastModel_FixParamValue("Denroll2", 0.0,4);
  mymodel.LastModel_FixParamValue("Denroll4", 0.0,4);
  mymodel.LastModel_FixParamValue("Denroll7", 0.0,4);
  mymodel.LastModel_FixParamValue("Denroll10", 0.0,4);
  mymodel.LastModel_FixParamValue("Denroll11", 0.0,4);
  mymodel.LastModel_FixParamValue("Denroll13", 0.0,4);


  // // Choice 5:
  mymodel.LastModel_FixParamValue("Denroll11", 0.0,5);
  mymodel.LastModel_FixParamValue("Denroll13", 0.0,5);

  // // Choice 6:
  mymodel.LastModel_FixParamValue("Denroll4", 0.0,6);
  mymodel.LastModel_FixParamValue("Denroll5", 0.0,6);
  mymodel.LastModel_FixParamValue("Denroll11", 0.0,6);
  mymodel.LastModel_FixParamValue("Denroll12", 0.0,6);
  mymodel.LastModel_FixParamValue("Denroll13", 0.0,6);
  
  // // Choice 7:
  mymodel.LastModel_FixParamValue("Denroll7", 0.0,7);

  // // Choice 8:
  mymodel.LastModel_FixParamValue("Denroll4", 0.0,8);
  mymodel.LastModel_FixParamValue("Denroll7", 0.0,8);
  mymodel.LastModel_FixParamValue("Denroll11", 0.0,8);

  // // Choice 9: 3, 6, 10 (4, 12)
  mymodel.LastModel_FixParamValue("Denroll4", 0.0,9);
  mymodel.LastModel_FixParamValue("Denroll5", 0.0,9);
  mymodel.LastModel_FixParamValue("Denroll7", 0.0,9);
  mymodel.LastModel_FixParamValue("Denroll11", 0.0,9);

  // // Choice 10: 4, 6, 11 (2, 3, 5, 12)
  mymodel.LastModel_FixParamValue("Denroll2", 0.0,10);
  mymodel.LastModel_FixParamValue("Denroll5", 0.0,10);
  mymodel.LastModel_FixParamValue("Denroll7", 0.0,10);
  mymodel.LastModel_FixParamValue("Denroll12", 0.0,10);

  // // Choice 11: 6 (3, 5, 10)
  mymodel.LastModel_FixParamValue("Denroll4", 0.0,11);
  mymodel.LastModel_FixParamValue("Denroll7", 0.0,11);
  mymodel.LastModel_FixParamValue("Denroll11", 0.0,11);

  // // Choice 12: 4, 6, 10 (2,3, 5)
  mymodel.LastModel_FixParamValue("Denroll5", 0.0,12);
  mymodel.LastModel_FixParamValue("Denroll7", 0.0,12);
  mymodel.LastModel_FixParamValue("Denroll11", 0.0,12);
  
  mymodel.StartNewPrintGroup();

  
  for (int imajor = 0; imajor < 12 ; imajor++) {
    
    if (imajor==6) mymodel.StartNewPrintGroup();

    // *** Estimate each model conditional on major:
    
    outcomevar = "DGradColl" + major[imajor];
    indvar = "ind_DGradColl" + major[imajor];
    name = "educDGradCollMaj" + major[imajor];
    title = "Graduate Major " + major[imajor];
    
    // cout << "Outcome: " << outcomevar.c_str() << endl;
    // cout << "Indicator: " << indvar.c_str() << endl;
    // cout << "Name: " << name.c_str() << endl;
    // cout << "Title: " << title.c_str() << endl << endl;;
    
    vector<TString> modelvars;
    modelvars.push_back(outcomevar.c_str());
    modelvars.push_back(indvar.c_str());
    modelvars.insert(modelvars.end(), bkgnd1990.begin(), bkgnd1990.end());
    modelvars.push_back("adveng");
    //    modelvars.push_back("adveng_miss");
    modelvars.push_back("advmath");
    //    modelvars.push_back("advmath_miss");
    modelvars.push_back("HS_track2");
    modelvars.push_back("HS_track3");
    mymodel.AddModel(name.c_str(),title.c_str(),"probit", modelvars);
    mymodel.LastModel_Detailsim();
    mymodel.LastModel_SetEndogenousReg(advengModelN);
    mymodel.LastModel_SetEndogenousReg(advmathModelN);
    mymodel.LastModel_SetEndogenousReg(HSTrackModelN);
  }

   */


  /*
  string income_meas[5] = {"logpvdispinc", "logwage", "logincome", "logpvincome", "logdispinc"};
  string education[17] = {"all", "HS", "CollDO_low", "CollDO_high", "CollGrad", "Maj1","Maj2", "Maj3","Maj4","Maj5","Maj6","Maj7","Maj8", "Maj9","Maj10","Maj11","Maj12"};


  for (int imeas = 0; imeas < 5 ; imeas++) {

    mymodel.StartNewPrintGroup();
    
    for (int ieduc = 0; ieduc < 17 ; ieduc++) {

    // *** Estimate each model conditional on major:
      mymodel.StartNewEstGroup();
      if (ieduc==5 || ieduc==11) mymodel.StartNewPrintGroup();
      
      outcomevar = income_meas[imeas];
      indvar = "ind_" + income_meas[imeas] + "_" + education[ieduc];
      name = income_meas[imeas] + education[ieduc];
      title = income_meas[imeas] + " (" + education[ieduc] +")";
      
      // cout << "Outcome: " << outcomevar.c_str() << endl;
      // cout << "Indicator: " << indvar.c_str() << endl;
      // cout << "Name: " << name.c_str() << endl;
      // cout << "Title: " << title.c_str() << endl << endl;;
      
      vector<TString> modelvars;
      modelvars.push_back(outcomevar.c_str());
      modelvars.push_back(indvar.c_str());
      modelvars.insert(modelvars.end(), bkgnd1990.begin(), bkgnd1990.end());
      modelvars.push_back("adveng");
      //      modelvars.push_back("adveng_miss");
      modelvars.push_back("advmath");
      //      modelvars.push_back("advmath_miss");
      modelvars.push_back("HS_track2");
      modelvars.push_back("HS_track3");
      mymodel.AddModel(name.c_str(),title.c_str(),"linear", modelvars);
      mymodel.LastModel_SetEndogenousReg(advengModelN);
      mymodel.LastModel_SetEndogenousReg(advmathModelN);
      mymodel.LastModel_SetEndogenousReg(HSTrackModelN);
      mymodel.LastModel_Detailsim();
    }
  }

  */
  //    mymodel.AllModel_Detailsim();
 //    mymodel.SimIncludeData();

  if ((mpirank==0)&&(mymodel.GetPrintLevel()>=1)) mymodel.PrintModels();
}
