#include <sys/types.h>  // For stat().
#include <sys/stat.h> 
#include <cmath>
#include <iostream>
#include <sstream>
#include <assert.h>

void calcgausshermitequadrature(int n, double * x, double * w);

int main() {


    double * getx = new double[20];
    double * getw = new double[20];

    for (int iquad = 2 ; iquad < 10; iquad++ ) {
      printf("*****Getting numbers for %2d GH quadrature points.\n",iquad);
      calcgausshermitequadrature(iquad,getx,getw);
      for (int ipoint = 0; ipoint < iquad/2 + 1; ipoint++)
	printf("points [%2d]= %22.18g %22.18g\n",iquad,getx[ipoint],getw[ipoint]);      
      printf("\n\n");
    }    
  return 0;
}

