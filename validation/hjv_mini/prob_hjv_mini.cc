// Mini-HJV problem for the legacy code. Compile with -DSTAGE=<1|2> -DNOBS=<N> [-DNTYPES=<T>].
// Mirrors model_r.R component by component. Groups: 0 = measurement, 1 = education, 2 = labor market.
// With NTYPES > 1 (T9): data files carry the suffix _t<T>, the measurement models fix their type intercepts at 0,
// and education and labor market form ONE final group that starts with the type probability model, as in bin/prob_hjv.cc.
#ifndef __CINT__
#include "Riostream.h"
#include "Rtypes.h"
#include <vector>
#include "TString.h"
#include "TMinLkhd.hh"
#endif
using namespace std;
#ifndef STAGE
#define STAGE 1
#endif
#ifndef NOBS
#define NOBS 2000
#endif
#define STR_(x) #x
#define STR(x) STR_(x)
#ifndef NTYPES
#define NTYPES 1
#endif
#if NTYPES > 1
#define DTAG "_t" STR(NTYPES)   // data file suffix
#define WTAG "t" STR(NTYPES)    // working directory tag
#elif defined(DSUF)
#define DTAG DSUF               // untyped model on typed data (legacy cannot run stage 1 with ntypes > 1 and a skipped type model)
#define WTAG DSUF
#else
#define DTAG ""
#define WTAG ""
#endif

// name, title, nfactors, ntypes, fcorr, fnmix, nquad, stage
#ifdef ADAPT
TMinLkhd mymodel("mini","mini",3,NTYPES,0,1,-8,STAGE);   // nquad < 0: adaptive integration, max 8 points
#else
TMinLkhd mymodel("mini","mini",3,NTYPES,0,1,8,STAGE);
#endif

static vector<TString> V(std::initializer_list<const char*> l) { vector<TString> v; for (auto s : l) v.push_back(s); return v; }

void init(int mpirank) {
  mymodel.SetMPRank(mpirank);
  mymodel.SetPrintLevel(1);
  #if defined(MEAS_ONLY)
  mymodel.SetWorkingDir("results/legacy_s" STR(STAGE) "m" WTAG "_" STR(NOBS) "/");
#elif defined(ADAPT)
  mymodel.SetWorkingDir("results/legacy_s" STR(STAGE) "a" ADAPT WTAG "_" STR(NOBS) "/");   // ADAPT = score source tag
#elif defined(PRIOR)
  mymodel.SetWorkingDir("results/legacy_s" STR(STAGE) "p" WTAG "_" STR(NOBS) "/");
#else
  mymodel.SetWorkingDir("results/legacy_s" STR(STAGE) WTAG "_" STR(NOBS) "/");
#endif
  mymodel.SetData((char *)"data/mini_" STR(NOBS) DTAG ".raw",(char *)"data/mini_varlist" DTAG ".txt");
  mymodel.SetObsIndex("id");
#ifdef ADAPT
  mymodel.SetFactorScores();            // reads factor_predictions.txt from the working directory
  mymodel.SetAdaptIntThresh(0.5);
#endif
#ifdef PRIOR
  mymodel.SetIncludePriorFactorScore();   // posterior-mode scores (R default); legacy default is ML scores
#endif

  // Loading normalizations: nfac entries, then NTYPES-1 zeros so that measurement models carry no type intercepts.
  auto NM = [](std::initializer_list<double> l) { vector<Double_t> v(l); for (int t = 1; t < NTYPES; t++) v.push_back(0.0); return v; };
  vector<Double_t> free3v = NM({-9999.0,-9999.0,-9999.0}), f1v = NM({-9999.0, 0.0, 0.0}), f12v = NM({-9999.0,-9999.0, 0.0});
  Double_t *free3 = free3v.data(), *f1 = f1v.data(), *f12 = f12v.data();

  // ---- group 0: measurement system ----
  auto X  = [](const char* y, const char* ind) { return V({y, ind, "intercept","x1","x2"}); };
  auto XH = [](const char* y, const char* ind) { return V({y, ind, "intercept","x1","x2","HS_track2","HS_track3","HS_track4"}); };
  vector<TString> m;
  m = X("g_eng9","none");    mymodel.AddModel("g_eng9","g_eng9","linear",m,free3);
  m = X("g_math9","none");   mymodel.AddModel("g_math9","g_math9","linear",m,free3);
  m = X("g_swe9","none");    mymodel.AddModel("g_swe9","g_swe9","linear",m,free3);
  m = X("g_sports9","none"); mymodel.AddModel("g_sports9","g_sports9","linear",m,free3);
  m = X("g_gpa9","none");    mymodel.AddModel("g_gpa9","g_gpa9","linear",m,free3);
  m = XH("g_swe10","ind_hs");    mymodel.AddModel("g_swe10","g_swe10","linear",m,free3);
  m = XH("g_sports10","ind_hs"); mymodel.AddModel("g_sports10","g_sports10","linear",m,free3);
  m = XH("g_gpahs","ind_hs");    mymodel.AddModel("g_gpahs","g_gpahs","linear",m,free3);
  m = X("cog1","ind_cog"); mymodel.AddModel("cog1","cog1","linear",m,f1);
  m = X("cog2","ind_cog"); mymodel.AddModel("cog2","cog2","linear",m,f1);
  m = X("cog3","ind_cog"); mymodel.AddModel("cog3","cog3","linear",m,f1);
  m = X("misslead","none"); mymodel.AddModel("misslead","misslead","probit",m,f1);
  m = X("lead","ind_lead");      mymodel.AddModel("lead","lead","linear",m,f12);
  m = X("stresstol","ind_lead"); mymodel.AddModel("stresstol","stresstol","linear",m,f12);

#ifndef MEAS_ONLY
  // ---- group 1: education (with NTYPES > 1: the single final group, opened by the type probability model) ----
  mymodel.StartNewEstGroup();
#if NTYPES > 1
  m = V({"none","intercept"}); mymodel.AddTypesModel(m);   // P(type | f): logit with an intercept and free factor loadings
#endif
  auto XZ = [](const char* y, const char* ind) { return V({y, ind, "intercept","x1","x2","advmath","adveng","HS_track3","HS_track4"}); };
  m = X("advmath","none"); mymodel.AddModel("advmath","advmath","probit",m);
  m = X("adveng","none");  mymodel.AddModel("adveng","adveng","probit",m);
  m = V({"HS_track","none","intercept","x1","x2","advmath","adveng"});
  mymodel.AddModel("HS_track","HS_track","logit",m,NULL,4);
  if (STAGE != 1) mymodel.LastModel_FixParamValue("x2", 0.0, 4);  // choice 4 = third non-reference alternative.
  // In stage 1 AddModel skips outcome models, so an unguarded call would bind to the last measurement model.
  m = XZ("apply","none");     mymodel.AddModel("apply","apply","probit",m);
  m = XZ("app1","apply");     mymodel.AddModel("application","application","logit",m,NULL,5,5);  // app1..app5 adjacent
  m = XZ("Dswesat","none");   mymodel.AddModel("Dswesat","Dswesat","probit",m);
  m = XZ("swesat","Dswesat"); mymodel.AddModel("swesat","swesat","linear",m);
  m = XZ("enroll","apply");   mymodel.AddModel("enroll","enroll","probit",m);
  m = XZ("switch_major","ind_enroll"); mymodel.AddModel("switch_major","switch_major","logit",m,NULL,4);
  m = XZ("grad1","ind_enroll"); mymodel.AddModel("grad1","grad1","probit",m);
  m = XZ("grad2","ind_enroll"); mymodel.AddModel("grad2","grad2","probit",m);
  m = XZ("grad3","ind_enroll"); mymodel.AddModel("grad3","grad3","probit",m);

  // ---- group 2: labor market (joint with education when NTYPES > 1) ----
#if NTYPES == 1
  mymodel.StartNewEstGroup();
#endif
  m = V({"wage1","none","intercept","x1","x2","HS_track3","HS_track4"}); mymodel.AddModel("wage1","wage1","linear",m);
  m = V({"wage2","none","intercept","x1","x2","HS_track3","HS_track4"}); mymodel.AddModel("wage2","wage2","linear",m);

#endif  // MEAS_ONLY
  if (mpirank==0) mymodel.PrintModels();
}
