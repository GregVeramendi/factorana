#ifndef TMODEL_HH
#define TMODEL_HH

#include "Rtypes.h" // for Double32_t
#include "TNamed.h"
#include "TString.h"
#include <vector>
#include <iostream>
#include <random>
#include <list>
#include <utility>

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
  Int_t numtyp;
  Int_t numchoice; // for multivariate models
  Int_t numrank; // for exploded logit model
  Int_t ranksharevar;
  std::list<std::pair<int,int>  > excludelist; //list of regressors to exclude from certain ranked choices

  Bool_t ignore;

  Double_t simresult;
  Double_t simresidual;
  std::vector<Int_t> endogVarList;
  std::vector<Int_t> endogRegList;
  std::vector<Int_t> endogModelList;
  std::vector<Int_t> endogChoiceList;
  std::vector<Double_t> endogXtiles;
  
  //static std::mt19937 mt{std::random_device{}()};

  public:
  TModel() {};
  TModel(const char *name, const char *title, Int_t modeltype, Int_t modelgroup, Int_t prntgroup, std::vector<Int_t> & moddata, Int_t nfac, Int_t ntyp, Double_t * thisnormfac = NULL, UInt_t nchoice = 2, UInt_t nrank = 1);
  virtual ~TModel();

  void SetIgnore(Bool_t ig) {ignore=ig;}
  Bool_t GetIgnore() {return ignore;}
  Int_t GetNreg() {return nregressors;}
  Int_t GetNchoice() {return numchoice;}
  Int_t GetOutcome() {return outcome;}
  Int_t GetMissing() {return missing;}
  Int_t GetType() {return modtype;}
  Int_t GetNrank() {return numrank;}
  Int_t GetGroup() {return modgroup;}
  Int_t GetPrintGroup() {return printgroup;}
  //  Int_t * GetReg() {return regressors;}
  std::vector<Int_t> GetReg() {return regressors;}
  std::vector<Double_t> GetNorm() {return facnorm;}
  Int_t GetReg(UInt_t i) {return regressors[i];}
  Int_t GetSplitSim() {return (splitsim>-1);}
  Bool_t ModelHasXtiles() {return (endogXtiles.size()>0);}
  
  void DetailSim(UInt_t ivar) {detailsim = ivar;}
  UInt_t GetDetailSim() {return detailsim;}
  Double_t GetSimResult() {return simresult;}
  void ClearSimResult() {simresult = -9999.0;}

  //  void SetRankedChoice(UInt_t Nrank; UInt_t rankshare = -9999.0) { if (modtype==3) { numrank = Nrank; ranksharevar = rankshare;} else assert(0) }

  void SetRankShareVar(UInt_t sharevarnum);
  void ExcludeRegressorRank(Int_t regnum, Int_t rank) {excludelist.push_back(std::make_pair(regnum,rank));}
  void SetXtileThresholds(std::vector<Double_t> & thresh) {endogXtiles = thresh;}
  void AddEndogenousReg(UInt_t endogvarnum) {endogVarList.push_back(endogvarnum);}
  void SetEndogenousRegs(std::vector<TModel> & models,std::vector<TString> & vartab);
  std::vector<Int_t> GetEndogenousReg() {return endogRegList;}
  void SplitSim(UInt_t ivar);
  void PrintModel(std::vector<TString> vartab);
  void Eval(UInt_t iobs_offset, const std::vector<Double_t> & data, const std::vector<Double_t> & param, UInt_t firstpar, const std::vector <Double_t> & fac, std::vector<Double_t> & modeval, std::vector<Double_t> & hess, Int_t gradflag);

  void Sim(UInt_t iobs_offset, const std::vector<Double_t> & data, std::vector<TModel> & models, const std::vector<Double_t> & param, UInt_t firstpar, const std::vector <Double_t> & fac, FILE * pfile, UInt_t gof = 0);
  Double_t GetPdf(UInt_t iobs_offset, const std::vector<Double_t> & data, const std::vector<Double_t> & param, UInt_t firstpar, const std::vector <Double_t> & fac);

//   std::vector<Double_t> Eval(UInt_t iobs_offset, Double_t * data, std::vector<Double_t> param, UInt_t firstpar, std::vector <Double_t> fac);
//  void Print(TString * vartab);

//  ClassDef(TModel,1) //TModel
};

#endif
