#include "minlkhd_nlp.hh"

#include "Rtypes.h" // for Double32_t
#include <cassert>
#include <iostream>
#include <cmath>


void LkhdFcn(Int_t &npar, Double_t *gin, Double_t &f,
		   Double_t *par, Int_t iflag, Double_t *hess = NULL);


using namespace Ipopt;

// constructor
MINLKHD_NLP::MINLKHD_NLP()
{}


MINLKHD_NLP::MINLKHD_NLP(const int npar, double * daXLo, double * daXUp, 
			double *daXInit, 
			int Nhess, int * naHessRows, int * naHessCols)
{
  nfreeparam = npar;

  Xinit.resize(nfreeparam,-9999.0);	
  XLo.resize(nfreeparam,-9999.0);
  XUp.resize(nfreeparam,-9999.0);

  for (int ipar = 0 ; ipar < npar ; ipar++ ) {
      Xinit[ipar] = daXInit[ipar];	
      XLo[ipar] = daXLo[ipar];
      XUp[ipar] = daXUp[ipar];
  }

  //  printf("Ipopt bounds are %8.5f and %8.5f and I'm using %8.5f\n",nlp_upper_bound_inf,nlp_lower_bound_inf,2e19);


  nHessElem = Nhess;
  HessCols.resize(Nhess,-9999);
  HessRows.resize(Nhess,-9999);

  if (nHessElem!=(nfreeparam*(nfreeparam+1))/2) {
      printf("ERROR in MINLKHD_NLP Constructor! NHessElem is wrong!\n");
      assert(0);  
  }

  for (int ihess = 0; ihess < nHessElem ; ihess++) {
     HessCols[ihess] = naHessCols[ihess];
     HessRows[ihess] = naHessRows[ihess];
  }
}

//destructor
MINLKHD_NLP::~MINLKHD_NLP()
{}

// returns the size of the problem
bool MINLKHD_NLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                             Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  //  printf("MINLKHD_NLP::get_nlp_info\n");
  // The problem described in MINLKHD_NLP.hpp has 4 variables, x[0] through x[3]
  n = nfreeparam;

  // one equality constraint and one inequality constraint
  m = 0;

  // in this example the jacobian is dense and contains 8 nonzeros
  nnz_jac_g = 0;

  // the hessian is also dense and has 16 total nonzeros, but we
  // only need the lower left corner (since it is symmetric)
  nnz_h_lag =  nHessElem;

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  return true;
}

// returns the variable bounds
bool MINLKHD_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  //  printf("MINLKHD_NLP::get_bounds_info\n");
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == nfreeparam);
  assert(m == 0);

  //  printf("Getting lower bounds\n");
  // the variables have lower bounds of 0.0
  for (Index i=0; i<n; i++) {
    x_l[i] = XLo[i];
    //    printf("%5d %8.3f %8.3f\n",i,XLo[i],x_l[i]);
  }

  //  printf("Getting bounds\n");
  // the variables have no upper bound
  for (Index i=0; i<n; i++) {
    x_u[i] = XUp[i];
    //    printf("%5d %8.3f %8.3f %8.3f %8.3f\n",i,XLo[i],x_l[i],XUp[i],x_u[i]);
  }

//   // the first constraint g1 has a lower bound of 25
//   g_l[0] = 0.0;
//   // the first constraint g1 has NO upper bound, here we set it to 2e19.
//   // Ipopt interprets any number greater than nlp_upper_bound_inf as
//   // infinity. The default value of nlp_upper_bound_inf and nlp_lower_bound_inf
//   // is 1e19 and can be changed through ipopt options.
//   g_u[0] = 0.0;

  // the second constraint g2 is an equality constraint, so we set the
  // upper and lower bound to the same value
  //  g_l[1] = g_u[1] = 40.0;
//   for (int iconstr = 0 ; iconstr < nconstr; iconstr++) {	
//       g_l[iconstr] =  g_u[iconstr] = 0.0;	
//   }
  return true;
}

// returns the initial point for the problem
bool MINLKHD_NLP::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);
  assert(n == nfreeparam);

  //  printf("MINLKHD_NLP::get_starting_point\n");
  // initialize to the given starting point
  for (Index i=0; i<n; i++) {
    x[i] = Xinit[i];
    //    printf("Starting par%d=%8.5f %8.5f\n",i,x[i],Xinit[i]);
  }

  return true;
}

// returns the value of the objective function
bool MINLKHD_NLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == nfreeparam);
  
  //  printf("MINLKHD_NLP::eval_f\n");

	
  Double_t obj;
  Double_t * grad = new Double_t[nfreeparam];
  Double_t * hess = new Double_t[nHessElem];

  Double_t * par = new Double_t[nfreeparam]; 
  for (int ipar = 0 ; ipar < nfreeparam ; ipar++) {
	par[ipar] = x[ipar];
  }

  LkhdFcn(n,grad,obj,par,1,hess);

  obj_value = obj;
  //  printf("Eval_F: f=%8.5f\n",obj_value);

  delete [] par;
  delete [] grad;
  delete [] hess;

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool MINLKHD_NLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == nfreeparam);

  //  printf("MINLKHD_NLP::eval_grad_f\n");

  Double_t obj;
  Double_t * grad = new Double_t[nfreeparam];
  Double_t * hess = new Double_t[nHessElem];

  Double_t * par = new Double_t[nfreeparam]; 
  for (int ipar = 0 ; ipar < nfreeparam ; ipar++) {
	par[ipar] = x[ipar];
	//	printf("Param#%5d=%8.4f %8.4f\n",ipar, par[ipar], x[ipar]);
  }

  LkhdFcn(n,grad,obj,par,2,hess);

  for (int ipar = 0 ; ipar < nfreeparam ; ipar++) {
        grad_f[ipar] = grad[ipar];
	//	printf("Grad Par# %5d:%8.4f %8.4f\n",ipar, grad_f[ipar], grad[ipar]);
  }


  delete [] par;
  delete [] grad;
  delete [] hess;

  return true;
}

// return the value of the constraints: g(x)
bool MINLKHD_NLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  assert(n == nfreeparam);
  assert(m == 0);
  printf("MINLKHD_NLP::eval_g\n");
  return false;
}

// return the structure or values of the jacobian
bool MINLKHD_NLP::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
  assert(m==0);

  //  printf("MINLKHD_NLP::eval_jac_g");
  // if (values==NULL) printf(": values=NULL\n");
  // else printf("\n");
  //  printf("Eval_jac_g with m=%d\n",m);
  // if (values == NULL) {
  //   // return the structure of the jacobian
  // 	assert(0);    
  // }
  // else {
  //   // return the values of the jacobian of the constraints
  // 	assert(0);
  // }

  return false;
}

//return the structure or values of the hessian
bool MINLKHD_NLP::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{

  assert(n == nfreeparam);
  assert(nele_hess == nHessElem);

  //  printf("MINLKHD_NLP::eval_h\n");
  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    // the hessian for this problem is actually dense

   for (int ihess = 0; ihess < nHessElem ; ihess++) {
     // Use upper right triangle
     // jCol[ihess] = HessCols[ihess];
     // iRow[ihess] = HessRows[ihess];

     // Use lower left triangle
     jCol[ihess] = HessRows[ihess];
     iRow[ihess] = HessCols[ihess];
   }
  }
  else {
    // return the values. This is a symmetric matrix, fill the lower left
    // triangle only

    Double_t obj;
    Double_t * par = new Double_t[nfreeparam]; 
    Double_t * grad = new Double_t[nfreeparam];
    Double_t * hess = new Double_t[nHessElem];

    for (int ipar = 0 ; ipar < nfreeparam ; ipar++) {
	par[ipar] = x[ipar];
    }

    LkhdFcn(n,grad,obj,par,3,hess);
 
    for (int ipar = 0 ; ipar < nHessElem ; ipar++) {
	values[ipar] = hess[ipar];
	//	printf("Hess element %4d: %8.5f\n",ipar,values[ipar]); 
    }

    delete [] par;
    delete [] grad;
    delete [] hess;
  }

  return true;
}

void MINLKHD_NLP::finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L, const Number* z_U,
                                  Index m, const Number* g, const Number* lambda,
                                  Number obj_value,
				  const IpoptData* ip_data,
				  IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

  printf("MINLKHD_NLP::finalize_solution\n");
  finalparam.resize(n);
  // For this example, we write the solution to the console
  //  std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
  for (Index i=0; i<n; i++) {
    //     std::cout << "x[" << i << "] = " << x[i] << std::endl;
     finalparam[i] = x[i];
  }

//   std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
//   for (Index i=0; i<n; i++) {
//     std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
//   }
//   for (Index i=0; i<n; i++) {
//     std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
//   }

   std::cout << std::endl << std::endl << "Final Objective value" << std::endl;
   std::cout << "f(x*) = " << obj_value << std::endl;
//  cost = obj_value;

//   std::cout << std::endl << "Final value of the constraints:" << std::endl;
//   for (Index i=0; i<m ;i++) {
//     std::cout << "g(" << i << ") = " << g[i] << std::endl;
//   }
}

void MINLKHD_NLP::getparam(double * thisparam) { 
  for (UInt_t ipar = 0; ipar < finalparam.size() ; ipar++) thisparam[ipar] = finalparam[ipar];
}
