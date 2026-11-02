#ifndef TMINLKHD_HH
#define TMINLKHD_HH

#include "Rtypes.h" // for Double32_t
#include "TNamed.h"
#include "TModel.hh"
#include "TString.h"

#include <vector>

void SetBitOn(UInt_t & var, UInt_t shift);
void SetBitOff(UInt_t & var, UInt_t shift);
UInt_t GetBit(UInt_t var, UInt_t shift);

//class TModel;

class TMinLkhd : public TNamed {
   
private:
 
  Int_t printlvl;
  Int_t printmargeffect;
  Int_t cpurank;
  Int_t counter;
  Bool_t initializing;
  Bool_t predicting;
  UInt_t seed;
  TString workingdir;

  // variables for subsampling and bootstrap:
  Int_t current_sample;
  Int_t current_printgroup;
  Int_t nprintgroup;
  Int_t current_estimationgroup;

  std::vector <Double_t> EstGroup_lkhd;

  UInt_t newflag;  // bitmap for new settings: subsample(0), bootstrap(1), fixed param(2), ignore model(3), est_nmix(4)

  UInt_t nsubsamples;
  Double_t subsample_percent;
  //  UInt_t newskipobs;
  std::vector <UInt_t> skipobs;

  UInt_t nbootsamples;
  UInt_t bootstrapstart;
  //  UInt_t newbootstrap;
  std::vector <UInt_t> bootstrapobs;
  UInt_t nchsamples;
  
  UInt_t predictobs;
  std::vector<Double_t> facprediction;
  Bool_t includePriorFactorScore;
  
  //Settings for initial parameters
  Bool_t initBetaOLS;
  Bool_t initEstOutcomeLoadings;
  Bool_t initFixedLoadings;
  UInt_t HessStdErr;
  Double_t loadingMultiplier;
  Double_t initVariance; 
  //  Bool_t AutoRestartMeas;

  UInt_t EstSequentMix;

  Double_t max_cpu_limit;
  Double_t stoch_deriv_frac;
  UInt_t CalcHess;
  UInt_t HessFactorScores;
  
  Bool_t adapt_int;
  Double_t adapt_int_thresh;
  UInt_t nquad_points;
  std::vector<UInt_t> fac_npoints;
  UInt_t stage;
  std::vector<Double_t> x, xadapt, w, wadapt; //Quadrature constants         

  std::vector<std::vector<Double_t> >  fscore;
  std::vector<std::vector<Double_t> >  fstderr;
  
  UInt_t nfac;
  Int_t fac_corr;
  UInt_t fac_nmix;
  UInt_t est_nmix;
  UInt_t f_nvariance;  // number of variance parameters
  UInt_t nfac_param;
  std::vector<Int_t> norm_models; 

  UInt_t ntyp;
  UInt_t type_model;
  
  UInt_t simIncData;
  UInt_t simWithData;
  UInt_t sim_nobs;
  UInt_t sampleposterior;
  std::vector<UInt_t> sim_modelorder;
  
  UInt_t nobs;
  UInt_t nvar;
  std::vector<Double_t> data;
  std::vector<TString> var_table;
//   Double_t * data;
//   TString * var_table;
  Int_t indexvar;
  Int_t weightvar;

  UInt_t nparam;
  UInt_t nfreeparam;
  std::vector<TModel> models;
  std::vector<UInt_t> fparam_models;
  std::vector<UInt_t> nparam_models;
  std::vector<UInt_t> nobs_models;
  std::vector <Double_t> param;
  std::vector <Double_t> param_err;
  std::vector <Double_t> param_fixval;
  std::vector <UInt_t> parfixed;
  std::vector <Int_t> parconstrained;
  //  UInt_t newfixpar;
  //  UInt_t newmodignore;

  // Get parameters for mixture:
  Double_t Getfvar(Int_t imix,Int_t ivar) {return param[imix*(f_nvariance+nfac+1)+ivar];}
  Double_t Getfmean(Int_t imix,Int_t ifac) {return param[imix*(f_nvariance+nfac+1)+f_nvariance+ifac];}
  Double_t Getfw(Int_t imix) {return param[imix*(f_nvariance+nfac+1)+f_nvariance+nfac];}
  
//   Double_t Getfvar(Int_t imix,Int_t ivar) {return param.at(imix*(f_nvariance+nfac+1)+ivar);}
//   Double_t Getfmean(Int_t imix,Int_t ifac) {return param.at(imix*(f_nvariance+nfac+1)+f_nvariance+ifac);}
//   Double_t Getfw(Int_t imix) {return param.at(imix*(f_nvariance+nfac+1)+f_nvariance+nfac);}

  // Get indices for parameters for mixture:
  Int_t Getifvar(UInt_t imix,UInt_t ivar) {return imix*(f_nvariance+nfac+1)+ivar;}
  Int_t Getifmean(UInt_t imix,UInt_t ifac) {return imix*(f_nvariance+nfac+1)+f_nvariance+ifac;}
  Int_t Getifw(UInt_t imix) {return imix*(f_nvariance+nfac+1)+f_nvariance+nfac;}

  //  UInt_t Getx1i(UInt_t imix, UInt_t xi) {return (imix*nquad_points + xi);}
  UInt_t Getx2i(UInt_t imix, UInt_t x1i, UInt_t x2i) {
    if (fac_corr==0) return (imix*nquad_points + x2i);
    else return (imix*nquad_points*nquad_points + x1i*nquad_points + x2i);
  }

  //Set value of param
  void Setparam(UInt_t ipar, double value) {
    if (param_fixval.at(ipar)<-9998.0) param.at(ipar) = value;
  }

  void Setparam_err(UInt_t ipar, double value) {
    if (param_fixval.at(ipar)<-9998.0) param_err.at(ipar) = value;
  }

  //Function to get quadrature points
  Double_t GetHGQx(UInt_t nquad, UInt_t ipoint) { return x[(nquad-1)*nquad_points + ipoint];}
  Double_t GetHGQw(UInt_t nquad, UInt_t ipoint) { return w[(nquad-1)*nquad_points + ipoint];}
  
public:
  TMinLkhd();
  TMinLkhd(const char *name, const char *title, Int_t nfactors, Int_t ntypes = 1, Int_t fcorr = 1, Int_t fnmix=0, Int_t nquad = 8, Int_t thisstage = 0);
  virtual ~TMinLkhd();
  

//   UInt_t GetNewBootStrap() {return newbootstrap;}
//   void ResetBootStrap() {newbootstrap = 0;}
//  void SetBootStrap(Int_t * thisbootstrapobs) {newbootstrap = 1; for (UInt_t iobs = 0; iobs < nobs ; iobs++) bootstrapobs[iobs] = thisbootstrapobs[iobs];}

  void SetMPRank(int mpirank) {cpurank = mpirank;}
  Int_t GetMPRank() {return cpurank;}

  void SetPrintLevel(int printlevel) {printlvl = printlevel;}
  Int_t GetPrintLevel() {return printlvl;}

  void SetWorkingDir(TString thisdir);
  TString GetWorkingDir() {return workingdir;}

  void SetInitLoadingMultiple(Double_t newmult) {loadingMultiplier = newmult;}
  void SetInitVariance(Double_t initvar) {initVariance = initvar;}
  void SetInitFixedLoading() {initFixedLoadings=1;}
  void SetInitBetaOLS() {initBetaOLS=1;}
  void SetIncludePriorFactorScore() {includePriorFactorScore = 1;}
  void SetInitEstOutcomeLoadings() {initEstOutcomeLoadings=1;}
  
  void CalcHessStdErr(UInt_t useHess = 1) {HessStdErr = useHess;}
  //  void SetAutoRestartMeasEstimation(); {AutoRestartMeas=1;}


  void StartNewPrintGroup() {nprintgroup = 0; current_printgroup++;}
  void StartNewEstGroup() {current_estimationgroup++;}

  UInt_t GetNewBootStrap() {return GetBit(newflag,1);}
  void ResetBootStrap() {SetBitOff(newflag,1);}
  void SetBootStrap(Int_t * thisbootstrapobs) {SetBitOn(newflag,1); for (UInt_t iobs = 0; iobs < nobs ; iobs++) bootstrapobs[iobs] = thisbootstrapobs[iobs];}
  std::vector<UInt_t> GetBootStrapObs() {return bootstrapobs;}
  UInt_t GetNBootSample() {return nbootsamples;}
  void SetBootStrapStartSample(int BSstart) {bootstrapstart = BSstart;}
  UInt_t GetBootStrapStartSample() {return bootstrapstart;}

  void SetNCHSamples(int nsamples) {nchsamples = nsamples;}
  
  UInt_t GetNewSkipObs() {return GetBit(newflag,0);}
  void ResetSkipObs() {SetBitOff(newflag,0);}
  void SetSkipObs(Int_t * thisobsskip) {SetBitOn(newflag,0); for (UInt_t iobs = 0; iobs < nobs ; iobs++) skipobs[iobs] = thisobsskip[iobs];}
  std::vector<UInt_t> GetSkipObs() {return skipobs;}
  Int_t GetNSubSample() {return nsubsamples;}

  Int_t GetCurrentSample() {return current_sample;}
  void SetCurrentSample(Int_t thissample) {current_sample = thissample;}

  Int_t GetCounter() {return counter;}
  Int_t GetInitializing() {return initializing;}
  Int_t GetPredicting() {return predicting;}
  Double_t Getparam(Int_t ipar) {return param.at(ipar);}
  Double_t Getparam_err(Int_t ipar) {return param_err.at(ipar);}
  Int_t GetNparam() {return nparam;}
  Int_t GetNobs() {return nobs;}
  Int_t GetNfreeparam() { return nfreeparam;}
  Int_t GetNfac() {return nfac;}
  Int_t GetNfacparam() {return nfac_param;}
  Int_t GetNmodels() {return models.size();}
  Int_t GetModel_Nparam(Int_t i) {return nparam_models[i];}
  Int_t GetModel_type(Int_t i) {return models[i].GetType();}
  Int_t GetModel_nreg(Int_t i) {return models[i].GetNreg();}
  Int_t GetModel_outcome(Int_t i) {return models[i].GetOutcome();}
  Int_t GetModel_missing(Int_t i) {return models[i].GetMissing();}

  std::vector<Int_t> GetModel_reg(Int_t i) {return models[i].GetReg();}
  void  PrintModel(Int_t i) {models[i].PrintModel(var_table);}
  void  PrintModels() {for (UInt_t i = 0 ; i < models.size() ; i++) models[i].PrintModel(var_table);}
  
  void ClearDataMembers();
  void ResetFitInfo();

  void ResetNewFixPar() {SetBitOff(newflag,2);}
  UInt_t GetNewFixPar() {return GetBit(newflag,2);}
  void FixPar(UInt_t iparam) {if (parfixed.at(iparam)==0) { nfreeparam--; parfixed.at(iparam) = 1; SetBitOn(newflag,2);} }
  void FixParPerm(UInt_t iparam) {if (parfixed.at(iparam)==0) { nfreeparam--; parfixed.at(iparam) = 2; SetBitOn(newflag,2);}}
  void ReleasePar(UInt_t iparam) { if (parfixed.at(iparam)==1) {nfreeparam++; parfixed.at(iparam) = 0; SetBitOn(newflag,2);}}
  UInt_t GetFixPar(UInt_t iparam) {return parfixed.at(iparam);}
  std::vector<UInt_t> GetFixPar() {return parfixed;}
  void SetFixPar(Int_t * thisfixpar) {
    for (UInt_t ipar = 0; ipar < nparam ; ipar++) {
      if (thisfixpar[ipar]==1) FixPar(ipar); 
      else if (thisfixpar[ipar]==0) ReleasePar(ipar);
    }
  }
  void FixParamValue(int parnum, double value) {param_fixval.at(parnum) = value;}
  void LastModel_FixParamValue(TString fixvar, double value, int choice = -1);
  
  
  void SetIgnoreMod(UInt_t imod) {SetBitOn(newflag,3); models[imod].SetIgnore(1);}
  void SetAllIgnoreMod(Int_t * thisignoremod) { SetBitOn(newflag,3); for (UInt_t imod = 0; imod < models.size() ; imod++) models[imod].SetIgnore(thisignoremod[imod]);}
  void RemoveIgnoreMod(UInt_t imod) {SetBitOn(newflag,3); models[imod].SetIgnore(0);}
  void ResetNewIgnoreMod() {SetBitOff(newflag,3);}
  UInt_t GetNewIgnoreMod() {return GetBit(newflag,3);}
  std::vector<UInt_t> GetIgnoreMod() { std::vector<UInt_t> modlist; for (UInt_t imod = 0 ; imod < models.size() ; imod++) modlist.push_back(models[imod].GetIgnore()); return modlist; }
  Int_t GetNignoremod() {Int_t Nmod = 0; for (UInt_t imod = 0 ; imod < models.size() ; imod++) if (models[imod].GetIgnore()==1) Nmod++; return Nmod;}

  void ResetNewEstNMix() {SetBitOff(newflag,4);}
  UInt_t GetNewEstNMix() {return GetBit(newflag,4);}
  UInt_t GetEstNMix() {return est_nmix;}
  void SetEstNMix(UInt_t nmix) {
    SetBitOn(newflag,4); 
    est_nmix = nmix;
    for (UInt_t ipar = 0 ; ipar < nfac_param ; ipar++) 
      if (ipar < UInt_t(Getifmean(est_nmix-1,0))) ReleasePar(ipar); 
      else FixPar(ipar); 
  }

  void PrintParam(Int_t imod = -1);
  void PrintParamTab(int pmod = 0);
  void PrintParam_Varlist();
  void SetData(char * filename, char * var_table);
  void AddModel(const char *name, const char *title, TString modeltype, std::vector<TString> & moddata,Double_t * normfac = NULL, UInt_t nchoice = 2, UInt_t nrank = 1);
  void AddTypesModel(std::vector<TString> & typedata, Double_t * typenorm = NULL);
  void LastModel_Splitsim(TString splitvar);

  void ConstrainFactorCorrelations();
  //  void ConstrainLastBetaToModel(const Int_t imod, const TString cov);
  void ConstrainLastFactorLoadingToModel(const UInt_t targetmod, UInt_t ifac);

  void SubSample(UInt_t nsamples, Double_t keep_perc = 0.7);
  void BootStrap(UInt_t nsamples, UInt_t gensample = 0);
  void EstSequenttialMixtures() {EstSequentMix = 1;}
  UInt_t GetNewFlag() {return newflag;}
  void SetNewFlag(UInt_t nflag) {newflag = nflag;}

  void LastModel_Detailsim() {models.back().DetailSim(1);}
  void LastModel_SetXtileThresh(std::vector<Double_t> & threshlist) {models.back().SetXtileThresholds(threshlist);}
  void LastModel_SetEndogenousReg(TString endogvar);
  void LastModel_SetEndogenousMajor(UInt_t modelN);
  void LastModel_SetRankShareVar(TString sharevar);
  
  void AllModel_Detailsim() {for (UInt_t imod = 0 ; imod < models.size() ; imod++) models[imod].DetailSim(1);}
  void SimIncludeData() {simIncData = 1;}
  void SimWithData() {simWithData = 1;}
  void SimSamplePosterior() {sampleposterior = 1;}
  void SetSimNobs(UInt_t nobs) {sim_nobs = nobs;}
  void SetSimModelOrder (std::vector<TString> & modorder);
  
  void SetFactorScores();
  void SetFactorSpecificQuadPoints(std::vector<UInt_t> & quadlist);
  void SetAsymmetricQuadPoints(const int asym);
  void SetObsIndex(const TString index);
  void UseWeights(const TString weight);
  void PrintMargEffect() {printmargeffect=1;}
  void SetCPULimit(const double sec) { if (sec>0) max_cpu_limit = sec;}
  void SetStochasticDeriv(const double frac) { if (frac>0.0 && frac<1.0) stoch_deriv_frac = frac;}
  void SetAdaptIntThresh(const double thresh) { if (thresh>0.0) adapt_int_thresh = thresh;}
  void UseHessatFactorScores() {HessFactorScores = 1;}
  void DoNotUseHessian() {CalcHess = 0;}
  
  Int_t Minimize(Int_t printlevel = 1);
  Int_t Est_measurementsys(Int_t printlevel = 1);
  Int_t Est_outcomes(Int_t printlevel = 1);
  Int_t Simulate(Int_t printlevel = 1);
  Int_t PredictFactors(Int_t printlevel = 1);
  Int_t GenerateCHsamples(Int_t printlevel = 1);
  Int_t TestCode(Int_t printlevel = 1);
  void EvalLkhd(Double_t &f, Double_t *gradL, Double_t *hessL, Double_t *minpar, Int_t iflag, Int_t rank, Int_t np);
  void CalcLkhd(Double_t &f, Double_t *gradL, Double_t *hessL, Double_t *minpar, Int_t iflag, Int_t rank, Int_t np);
  Double_t CalcStdErrLkhdRatio(UInt_t iparam);
  Int_t Min_Minuit(Int_t printlevel);
  Int_t Min_Knitro(Int_t printlevel);
  Int_t Min_Ipopt(Int_t printlevel);
  Int_t TestGradient();
  Int_t TestHessian(); 
//  ClassDef(TMinLkhd,1) //TMinLkhd
};

#endif
