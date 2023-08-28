
/*************************************************************************
Computation of nodes and weights for a Gauss-Hermite quadrature formula

The  algorithm  calculates  the  nodes  and  weights  of the Gauss-Hermite
quadrature  formula on domain  (-infinity, +infinity) with weight function
W(x)=Exp(-x*x).

Input parameters:
    n   –   a required number of nodes.
            1 <= n <= 190.

Output parameters:
    x   -   array of nodes.
            Array whose index ranges from 0 to N-1.
    w   -   array of weighting coefficients.
            Array whose index ranges from 0 to N-1.

The algorithm was designed by using information from the QUADRULE library.
*************************************************************************/
//#include <iostream.h>
#include <cmath>

using namespace std;


void calcgausshermitequadrature(int n,
     double * x,
     double * w)
{
    int i;
    int j;
    double r = 0;
    double r1;
    double p1;
    double p2;
    double p3;
    double dp3;
    double pipm4;
    double pi = 3.14159265358979323846;
    
    pipm4 = pow(pi, -0.25);
    for(i = 0; i <= (n+1)/2-1; i++)
    {
      //      printf("Calculating element %d \n",i);
        if( i==0 )
        {
            r = sqrt(double(2*n+1))-1.85575*pow(double(2*n+1), -double(1)/double(6));
        }
        else
        {
            if( i==1 )
            {
                r = r-1.14*pow(double(n), 0.426)/r;
            }
            else
            {
                if( i==2 )
                {
                    r = 1.86*r-0.86*x[0];
                }
                else
                {
                    if( i==3 )
                    {
                        r = 1.91*r-0.91*x[1];
                    }
                    else
                    {
                        r = 2*r-x[i-2];
                    }
                }
            }
        }
        do
        {
            p2 = 0;
            p3 = pipm4;
            for(j = 0; j <= n-1; j++)
            {
                p1 = p2;
                p2 = p3;
                p3 = p2*r*sqrt(double(2)/double(j+1))-p1*sqrt(double(j)/double(j+1));
            }
            dp3 = sqrt(double(2*j))*p2;
            r1 = r;
            r = r-p3/dp3;
	    //	    printf("calculating x%d: %e, diff=%e\n",i,r,fabs(r-r1));
        }
       while((fabs(r-r1)>=fabs(r)*1e-15)&&(fabs(r-r1)>=1e-100));
        x[i] = r;
        w[i] = 2/(dp3*dp3);
        x[n-1-i] = -x[i];
        w[n-1-i] = w[i];
    }
}



