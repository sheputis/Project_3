#include <iostream>
#include "armadillo"
#include "lib.h"
using namespace std;
using namespace arma;

double Legendre(int n, double x); // this is a function copied from lecture notes page 121:
double int_function(double x1, double y1, double z1, double x2, double y2, double z2);

int main()
{

    int N=30;
    double a=-3;double b=3;

    //reserve memory for mesh points x and weights w
    double *x=new double[N];
    double *w=new double[N];
    //set up mesh points and weights:
    gauleg(a,b,x,w,N);
    /*
    for(int p=0;p<N;p++){
        cout<< x[p]<<"  "<<w[p]<<endl;
    }
    */
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
    cout << int_gauss<<endl<< pow(2,1);


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





















double Legendre(int n, double x)// computes Legendre polynomial of degree N
{
    double s,r,t; //s is Lj+1(x); r is Lj(x); t is Lj-1(x)
    int m;
    // Use recursion relation to generate p1 and p2
    for (m=0; m<n; m++)
    {
        t=r;r=s;
        s=(2*m+1)*x*r - m*t;
        s /=(m+1);
    }//end of do loop
}//end of function Legendre


double L_matrix(double **A)
{

}

