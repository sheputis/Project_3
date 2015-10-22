
// Arnoldas Seputis , project 3. The different sections of code, a, b, c and d can be run independently by uncommenting them
#include <iostream>
#include <math.h>

//#include "armadillo"
#include "lib.h"
#include "gauss-laguerre.cpp"

// MC
#include <fstream>
#include <iomanip>


// keeping track of time
#include "time.h"

using namespace std;
//using namespace arma;

double int_function(double x1, double y1, double z1, double x2, double y2, double z2);
double int_function_B(double phi_1, double phi_2, double th_1, double th_2, double r1, double r2);
double int_function_D(double phi_1, double phi_2, double th_1, double th_2, double r1, double r2);

int main()
{

    //TIME:
    clock_t start;
    clock_t finish;
    start = clock();
    //                                           PROBLEM A: SOLVED USING LEGENDRE
/*
    int N=30;

    double a=-2.5;double b=2.5;

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
    double exact= (5./(16.*16.))*M_PI*M_PI;
    cout <<"Problem A, the integral using Legendre : "<< int_gauss<<endl;
    cout <<" Number of points per dimension : " << N<<endl;
    cout<< " Exact = " << exact <<endl;
    cout<< " Relative error = " <<(int_gauss-exact)/exact;
    */

    /*                                       PROBLEM B !!!!!!! Laguerre !!!!!!!!!!!!!!                   */
/*
    int N= 30;
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

//  mesh points for legendre -0 to 2pi : = angle phi
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
    int_gauss=int_gauss/(pow(4.,5));
   // cout << int_gauss/(pow(4.,5))<<endl;//(pow(2.,6))<<
    //final output:
    double exact= (5./(16.*16.))*M_PI*M_PI;
    cout <<"Problem B, the integral using Legendre and Laguerre : "<< int_gauss<<endl;
    cout <<" Number of points per dimension : " << N <<endl;
    cout<< " Exact = " << exact <<endl;
    cout<< " Relative error = " <<(int_gauss-exact)/exact;
*/

    //                                           PROBLEM C
/*
    int n= 7.3*pow(10.,3);
    long idum = -1;
    double int_mc=0; double x[6]; double sum_sigma=0;double fx; double variance=0.;
    double length=2.5;double jacobidet=pow((2*length),6);
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
    variance=abs(sum_sigma-int_mc*int_mc);
    int_mc=jacobidet*int_mc;
    //final output:
    double exact=(5./(16.*16.))*M_PI*M_PI;
    cout<< " n :" << n <<endl;
    cout<< setiosflags(ios::showpoint | ios::uppercase);
    cout<< " Monte carlo result " << setw(10) << setprecision(8)<< int_mc <<endl;
    cout<< " Exact = " << exact <<endl;
    cout<< " Sigma = " << setw(10)<< setprecision(8)<< jacobidet*sqrt(variance/((double) n))<<endl;
    cout<< " Relative error =" << (int_mc-exact)/exact<<endl;


*/






    //                                           PROBLEM D
/*
    int n= 7.3*pow(10.,8);
    long idum = -1;
    double int_mc=0; double x[6]; double sum_sigma=0;double fx; double variance=0.;
    double jacobidet=(M_PI*M_PI)/(pow(4.,3));
    //evaluating integral with the crude Montecarlo method:
    for (int i= 1;i<=n;i++){

        x[0]=ran0(&idum);  // this is x1 [0,1] converted from r1 [0,inf)
        x[1]=ran0(&idum);  // x2
        //converting to y:
        x[0]=-log(1-x[0]);
        x[1]=-log(1-x[1]);

        x[2]=-1+2*ran0(&idum); //dcos_theta_1
        x[3]=-1+2*ran0(&idum); //dcos_theta_2
        x[4]=2*M_PI*ran0(&idum); //d_phi_1
        x[5]=2*M_PI*ran0(&idum); //d_phi_2

      //  cout<<x[0]<<endl;


        fx=int_function_D(x[4],x[5],x[2],x[3],x[0],x[1]);
      //  fx=int_function_B(x[0],x[1],x[2],x[3],x[4],x[5]);
        int_mc+=fx;
        sum_sigma +=fx*fx;

    }
    int_mc = int_mc/((double) n);
    sum_sigma = sum_sigma/((double) n);
    variance=sum_sigma-int_mc*int_mc;

    double exact=(5./(16.*16.))*M_PI*M_PI;
    //final output:
    cout<< " n :" << n<< endl;
    cout<< setiosflags(ios::showpoint | ios::uppercase);
    cout<< " Monte carlo result " << setw(10) << setprecision(8)<< jacobidet*int_mc <<endl;
    cout<< " Exact = " << exact <<endl;
    cout<< " Sigma = " << setw(10)<< setprecision(8)<< jacobidet*sqrt(variance/((double) n))<<endl;
    cout<< " Relative error =" << (jacobidet*int_mc-exact)/exact<<endl;


 */
    //.................... end of problem 4

    finish = clock();
    double time = ((finish-start)/CLOCKS_PER_SEC);
    cout<< "Time used: " << time << endl;
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

        value=0;

    }
    else
    {

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




double int_function_D(double phi_1, double phi_2, double th_1, double th_2, double r1, double r2)
{

    double numerator=(r1*r1)*(r2*r2);

    double cos_Bet;
    cos_Bet=th_1*th_2+abs(sqrt(1-th_1*th_1))*abs(sqrt(1-th_2*th_2))*cos(phi_1-phi_2);
    double denominator=sqrt(r1*r1+r2*r2-2*r1*r2*cos_Bet);
    double value;
    if(denominator< pow(10.,-6))
    {
        value =0;//cout<< r12<< " " << (r12< pow(10.,-6))<< " "<< pow(10.,-6)<<endl;


    }
    else{value=numerator/denominator;}
    return value;
}









