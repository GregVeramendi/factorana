#ifndef TMODEL_HH
#define TMODEL_HH

#include "Rtypes.h" // for Double32_t
#include "TNamed.h"
#include "TString.h"
#include <vector>
#include <iostream>

class TModel : public TNamed {
   
private:
 
  Int_t modtype;  // 1 = linear, 2 = probit, 3=logit, 4=ordered probit
  Int_t modgroup;  // 0 = measurement system, 1+ = outcomes to be estimated together
  Int_t printgroup;  // so that similar models end up in the same table
  Int_t outcome;
  Int_t missing;
  Int_t nregressors;
  Int_t splitsim;
  UInt_t detailsim;
  std::vector<Int_t> regressors;
  std::vector<Double_t> facnorm;
  Int_t numfac;
  Int_t numchoice; // for multivariate models
  Bool_t ignore;

  public:
  TModel() {};
  TModel(const char *name, const char *title, Int_t modeltype, Int_t modelgroup, Int_t prntgroup, std::vector<Int_t> & moddata, Int_t nfac, Double_t * thisnormfac = NULL, UInt_t nchoice = 2);
  virtual ~TModel();

  void SetIgnore(Bool_t ig) {ignore=ig;}
  Bool_t GetIgnore() {return ignore;}
  Int_t GetNreg() {return nregressors;}
  Int_t GetNchoice() {return numchoice;}
  Int_t GetOutcome() {return outcome;}
  Int_t GetMissing() {return missing;}
  Int_t GetType() {return modtype;}
  Int_t GetGroup() {return modgroup;}
  Int_t GetPrintGroup() {return printgroup;}
  //  Int_t * GetReg() {return regressors;}
  std::vector<Int_t> GetReg() {return regressors;}
  std::vector<Double_t> GetNorm() {return facnorm;}
  Int_t GetReg(UInt_t i) {return regressors[i];}
  Int_t GetSplitSim() {return (splitsim>-1);}

  void DetailSim(UInt_t ivar) {detailsim = ivar;}
  UInt_t GetDetailSim() {return detailsim;}

  void SplitSim(UInt_t ivar);
  void PrintModel(std::vector<TString> vartab);
  void Eval(UInt_t iobs_offset, const std::vector<Double_t> & data, const std::vector<Double_t> & param, UInt_t firstpar, const std::vector <Double_t> & fac, std::vector<Double_t> & modeval, std::vector<Double_t> & hess, Int_t gradflag);

  void Sim(UInt_t iobs_offset, const std::vector<Double_t> & data, const std::vector<Double_t> & param, UInt_t firstpar, const std::vector <Double_t> & fac, FILE * pfile);

//   std::vector<Double_t> Eval(UInt_t iobs_offset, Double_t * data, std::vector<Double_t> param, UInt_t firstpar, std::vector <Double_t> fac);
//  void Print(TString * vartab);

//  ClassDef(TModel,1) //TModel
};

#endif
