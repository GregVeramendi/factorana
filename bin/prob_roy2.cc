#ifndef __CINT__

#include "Riostream.h"
#include "Rtypes.h"
#include <iostream>
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

#endif

//TMinLkhd mymodel("Econ440 HW3","hw3",1,0,10,1);
TMinLkhd mymodel("2-factor Roy Model","2-factor Roy Model",2,1,1,4,1);

void init() {
  // TMinLkhd mymodel("2-factor Roy Model","2-factor Roy Model",2,0,6);
  mymodel.SetData(11,(char *)"roy_2fac_data.txt",(char *)"roy_2fac_var.txt");
  
  vector<TString> testC;
  testC.push_back("test");
  testC.push_back("none");
  testC.push_back("cons");
  testC.push_back("Q");
  Double_t normC[2] = {1.0,0.0};
  mymodel.AddModel(0,"testcog","Test Cog","linear",testC,normC);
  
  vector<TString> testSE1;
  testSE1.push_back("B1");
  testSE1.push_back("none");
  testSE1.push_back("cons");
  testSE1.push_back("Q");
  Double_t normSE1[2] = {0.0,1.0};
  mymodel.AddModel(0,"testse1","Test SE1","probit",testSE1,normSE1);
  
  vector<TString> testSE2;
  testSE2.push_back("B2");
  testSE2.push_back("none");
  testSE2.push_back("cons");
  testSE2.push_back("Q");
  Double_t normSE2[2] = {0.0,-9999};
  mymodel.AddModel(0,"testse2","Test SE2","probit",testSE2,normSE2);
  
  vector<TString> testSE3;
  testSE3.push_back("B3");
  testSE3.push_back("none");
  testSE3.push_back("cons");
  testSE3.push_back("Q");
  Double_t normSE3[2] = {0.0,-9999};
  mymodel.AddModel(0,"testse3","Test SE3","probit",testSE3,normSE3);
  
  vector<TString> sector;
  sector.push_back("sector1");
  sector.push_back("none");
  sector.push_back("cons");
  sector.push_back("Z");
  mymodel.AddModel(0,"sector","Sector","probit",sector);
  
  vector<TString> wages0;
  wages0.push_back("wage");
  wages0.push_back("sector0");
  wages0.push_back("cons");
  wages0.push_back("X");
  mymodel.AddModel(0,"wages0","Wage S0","linear",wages0);
  
  vector<TString> wages1;
  wages1.push_back("wage");
  wages1.push_back("sector1");
  wages1.push_back("cons");
  wages1.push_back("X");
  mymodel.AddModel(0,"wages1","Wage S1","linear",wages1);
  
   mymodel.PrintModels();
   mymodel.ResetFitInfo();
}
