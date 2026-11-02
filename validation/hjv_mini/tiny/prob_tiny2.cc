#ifndef __CINT__
#include "Riostream.h"
#include "Rtypes.h"
#include <vector>
#include "TString.h"
#include "TMinLkhd.hh"
#endif
using namespace std;

// name, title, nfactors, ntypes, fcorr, fnmix, nquad, stage (1 = measurement system)
TMinLkhd mymodel("tiny","tiny",1,1,0,1,8,1);

void init(int mpirank) {
  mymodel.SetMPRank(mpirank);
  mymodel.SetPrintLevel(1);
  mymodel.SetWorkingDir("results2/");
  mymodel.SetData((char *)"tiny.raw",(char *)"tiny_varlist.txt");

  Double_t norm[1]; norm[0] = 1.0;           // anchor: M1 loading fixed at 1, variance free
  vector<TString> m1 = {"M1","none","const","x1","x2"};
  vector<TString> m2 = {"M2","none","const","x1","x2"};
  vector<TString> m3 = {"M3","none","const","x1","x2"};
  mymodel.AddModel("m1","M1 model","linear",m1);
  mymodel.AddModel("m2","M2 model","linear",m2);
  mymodel.AddModel("m3","M3 model","linear",m3);
  if (mpirank==0) mymodel.PrintModels();
}
