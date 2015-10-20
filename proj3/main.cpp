#include <iostream>
#include <math.h>

//#include "armadillo"
#include "lib.h"
#include "gauss-laguerre.cpp"

// MC
#include <fstream>
#include <iomanip>



/*
//for laguerre:
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#define EPS 3.0e-14
#define MAXIT 10

*/

using namespace std;
//using namespace arma;

double Legendre(int n, double x); // this is a function copied from lecture notes page 121:
double int_function(double x1, double y1, double z1, double x2, double y2, double z2);
double int_function_B(double phi_1, double phi_2, double th_1, double th_2, double r1, double r2);


int main()
{
    //                                           PROBLEM A: SOLVED USING LEGENDRE
    /*
    int N=30;

    double a=-3;double b=3;

    //reserve memory for mesh points x and weights w
    double *x=new double[N];
    double *w=new double[N];
    //set up mesh points and weights:
    gauleg(a,b,x,w,N);

    //evaluating integral with gauss_legendre method:
    double int_gauss = 0.0;
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            for(int k = 0; k<N; k++){
                for(int l = 0; l<N; l++){
                    for(int m = 0; m<N; m++){
                        for(int n = 0; n<N; n++){
                               int_gauss+=w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*int_function(x[i],x[j],x[k],x[l],x[m],x[n]);
                        }
                    }
                }
        }
    }

    }
    //final output:
    cout <<"Problem A, the integral using Legendre"<< int_gauss<<endl;
     */
    /*                                       PROBLEM B !!!!!!! Laguerre !!!!!!!!!!!!!!                   */
    /*
    int N= 5;
    double *xgl = new double [N+1];
    double *wgl = new double [N+1];

//   set up the mesh points and weights and the power of x^alf for both r1 and r2
    double alf = 2.0;
    gauss_laguerre(xgl,wgl, N, alf);

//   mesh points for legendre -1 to 1: = anlge tetha
    double a_th=-1;double b_th=1;

    //reserve memory for mesh points x and weights w
    double *x_th=new double[N];
    double *w_th=new double[N];
    //set up mesh points and weights:
    gauleg(a_th,b_th,x_th,w_th,N);

//  mesh points for legendre -1 to 1: = angle phi
    double a_phi=0;double b_phi=2*M_PI;

    //reserve memory for mesh points x and weights w
    double *x_phi=new double[N];
    double *w_phi=new double[N];
    //set up mesh points and weights:
    gauleg(a_phi,b_phi,x_phi,w_phi,N);


    //evaluating integral with gauss_legendre and laugerre methods:
    double int_gauss = 0.0;
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            for(int k = 0; k<N; k++){
                for(int l = 0; l<N; l++){
                    for(int m = 0; m<N; m++){
                        for(int n = 0; n<N; n++){
                               int_gauss+=w_phi[i]*w_phi[j]*w_th[k]*w_th[l]*wgl[m]*wgl[n]*int_function_B(x_phi[i],x_phi[j],x_th[k],x_th[l],xgl[m],xgl[n]);

                        }
                    }
                }
        }
    }
}
    cout << int_gauss/(pow(4.,5))<<endl;//(pow(2.,6))<<
    */

    //                                           PROBLEM C
    /*
    int n= 10000000;
    long idum = -1;
    double int_mc=0; double x[6]; double sum_sigma=0;double fx; double variance=0.;
    double length=3.;double jacobidet=pow((2*length),6);
    //evaluating integral with the crude Montecarlo method:
    for (int i= 1;i<=n;i++){
        for(int j=0;j<6;j++){
            x[j]=-length+2*length*ran0(&idum);

        }
        fx=int_function(x[0],x[1],x[2],x[3],x[4],x[5]);
        int_mc+=fx;
        sum_sigma +=fx*fx;

    }
    int_mc = int_mc/((double) n);
    sum_sigma = sum_sigma/((double) n);
    variance=sum_sigma-int_mc*int_mc;

    //final output:
    cout<< setiosflags(ios::showpoint | ios::uppercase);
    cout<< " Monte carlo result " << setw(10) << setprecision(8)<< jacobidet*int_mc <<endl;
    cout<< " Exact = " << (5./(16.*16.))*M_PI*M_PI<<endl;
    cout<< " Sigma = " << setw(10)<< setprecision(8)<< jacobidet*sqrt(variance/((double) n))<<endl;

    */







    //                                           PROBLEM D


    return 0;
}

double int_function(double x1,double y1, double z1, double x2, double y2, double z2)
{
    double alpha=2; //corresponds to charge
    // radius lengths of the electron positions
    double r1=sqrt(x1*x1+y1*y1+z1*z1);

    double r2=sqrt(x2*x2+y2*y2+z2*z2);

    double numerator=exp(-2*alpha*(r1+r2));
    double denominator=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
    double value;
    if(denominator < pow(10.,-6.))
    {
       // cout<<"B";
        value=0;

    }
    else
    {
        //cout<< "A";
        value=numerator/denominator;
    }
    return value;
}

double int_function_B(double phi_1, double phi_2, double th_1, double th_2, double r1, double r2)
{
    double cos_Bet;double r12;
    cos_Bet=th_1*th_2+abs(sqrt(1-th_1*th_1))*abs(sqrt(1-th_2*th_2))*cos(phi_1-phi_2);

    r12=sqrt(r1*r1+r2*r2-2*r1*r2*cos_Bet);
    //cout<<r12<<endl;
    double value;
    if(r12< pow(10.,-6))
    {
        value =0;//cout<< r12<< " " << (r12< pow(10.,-6))<< " "<< pow(10.,-6)<<endl;


    }
    else{value=1./r12;}
    return value;
}



















