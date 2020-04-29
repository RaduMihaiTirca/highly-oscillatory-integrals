

# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <complex>
# include "Faddeeva.hh"
# include "matplotlibcpp.h"



using namespace std;

namespace plt = matplotlibcpp;


const complex<double> i(0,1);
const double pi      = acos(-1);
const double C       = 137.036;
const double k       = 0.05/C;
const double tau     = 5*C*(2.*pi/0.05);
const double A0      = C;
const double phiMax  = 2*tau;
const double phiMin  = -2*tau;
      double p       = 0.001;

double *AphiTab(int n, double x[]){

  double *fx;
  
  fx = new double[n+1];
  
  for (int j = 0; j <= n; j++){

    fx[j] = A0 * exp(-((x[j] * x[j]) / (tau * tau))) * sin(k * x[j]);

  }
  return fx;
}

complex<double> *complexAphiTab(int n, double x[]){

  complex<double> *fx;
  
  fx = new complex<double>[n+1];
  
  for (int j = 0; j <= n; j++){

    fx[j] = A0 * exp(-((x[j] * x[j]) / (tau * tau))) * sin(k * x[j]);

  }
  return fx;
}


double derS(double x){

  double tmp = A0 * exp( -1 * ( (x * x) / (tau * tau) ) ) * sin(k * x);

    return (1/(2.*C))*tmp*tmp + p;

}

complex<double> S(double phi){
  
  return p * phi + (1/(2*C)) * (1./8.) * exp(-0.5*k*k*tau*tau) * sqrt(pi/2.) *
          tau * ( -2 + 2 * exp( (k * k * tau * tau ) / 2.) * 
          (1 + Faddeeva::erf( (sqrt(2)*phi) / tau )) 
         - i * Faddeeva::erfi( (k * tau * tau - 2.* i * phi) / (sqrt(2)*tau)) 
         + i * Faddeeva::erfi( (k * tau * tau + 2.* i * phi) / (sqrt(2)*tau)) );
}

complex<double> *hasPhiTab (int n, double phi[]){

  complex<double> *fx;

  fx = new complex<double>[n];

  for(int j = 0; j < n; j++)
    fx[j] =  S(phi[j]) ;

  return fx;

}



complex<double> filonQuadrature(int n, double a, double b, double *f(int n , double x[]), double omega){


    double h = (b-a) / (2.*n);
    double *x;
    double *v;
    double *ftab;

    double theta = omega*h;
    double sint = sin(theta), cost = cos(theta);
    double alpha, beta, gamma;

    if ( 6.0 * fabs ( theta ) <= 1.0 )
  {
    alpha = 2.0 * pow ( theta, 3 ) /   45.0 
          - 2.0 * pow ( theta, 5 ) /  315.0 
          + 2.0 * pow ( theta, 7 ) / 4725.0;
  
    beta =  2.0                    /     3.0 
          + 2.0 * pow ( theta, 2 ) /    15.0 
          - 4.0 * pow ( theta, 4 ) /   105.0 
          + 2.0 * pow ( theta, 6 ) /   567.0 
          - 4.0 * pow ( theta, 8 ) / 22275.0;

    gamma = 4.0                    /      3.0 
          - 2.0 * pow ( theta, 2 ) /     15.0 
          +       pow ( theta, 4 ) /    210.0 
          -       pow ( theta, 6 ) /  11340.0;
  }
  else
  {
    alpha = ( pow ( theta, 2 ) + theta * sint * cost 
      - 2.0 * sint * sint ) / pow ( theta, 3 );

    beta = ( 2.0 * theta + 2.0 * theta * cost * cost
      - 4.0 * sint * cost ) / pow ( theta, 3 );

    gamma = 4.0 * ( sint - theta * cost ) / pow ( theta, 3 );
   }

    v = new double [2*n+1];

    for (int j = 0; j <= 2*n; j++)
        v[j] = (double) a + j*h; 

    x = new double [2*n+1];

    x[0] = phiMin;

    for (int j = 0; j < 2*n; j++){
        double tmp = x[j];
        x[j+1] = tmp + (omega/derS(tmp))*h;
    }
    
    cout << x[0] << " " << x[1] << "  " << h << " " ;

    ftab = f(2*n,x); 

    for (int j = 0; j <= 2*n; j++){
      ftab[j] *= ftab[j];
      ftab[j] /= derS(x[j]);
    }
    complex<double> sigma1(0,0);

    for (int j = 0; j<=n; j++)
        sigma1 += ftab[2*j] * exp(i*omega*v[2*j]);

    complex<double> sigma2(0,0);

    for (int j = 1; j<=n; j++)
        sigma2 += ftab[2*j-1] * exp(i*omega*v[2*j-1]);



    complex<double> result = h * omega * (i * alpha * (ftab[0] - ftab[2*n] * exp(i* omega)) 
                                   
                                    + beta  * (sigma1 - 0.5 * (ftab[0] + ftab[2*n] * exp(i* omega)))
                                   
                                    + gamma * sigma2

                                  );

    delete [] ftab;
    delete [] x;

    return result;

}

complex<double> trapezoidalMethod(int n, double a, double b, complex<double> *f(int n , double x[])){
  
      
  double h = (double) (b-a) / n;
  double *x;
  complex<double> *hasTab;
  complex<double> *ftab;

  hasTab = new complex<double>[n+1];
  x   = new double[n+1];

    for (int j = 0; j <= n; j++)
        x[j] = (double) a + (double) h*j; 

    ftab = f(n,x); 
    hasTab = hasPhiTab(n+1,x) ;
               
     for (int j = 0 ; j < n+1 ; j ++){                              
        ftab[j] *= ftab[j];
        ftab[j] *= exp(i*hasTab[j]) ;
     }
  complex<double> result(0,0) ;

    for (int j = 1; j < n ; j++)
      result += ftab[j]*h;
  
  result += h/2.*(ftab[0] +ftab[n]);

  delete [] ftab;
  delete [] x;

  return result;

}

int main(){

    int n = 1000000;
    double a = phiMin;
    double b = phiMax;
    ofstream file;
    double *tmp;
    double *g;
    
    double h = (b-a)/(n);
    
    g = new double[n];
    
    for (int j = 0; j < n ; j ++){
        g[j] = a+ j*h;
    }
    tmp = AphiTab(n, g);
    
    file.open("output/filon0.8/sampleOutput.txt");

    vector<double> x(n),y(n); 
    
    for(int j=0; j<n; ++j) {
        x.at(j) = a + h*j;
        y.at(j) = real(S(x.at(j)));//tmp[j]*tmp[j]*cos(real(S(x.at(j))));
       // y.at(j) = real(S(x.at(j)));//+x.at(j);
    }
    plt::plot(x, y);

    // complex<double> filonRes, trapRes;
    // filonRes = filonQuadrature(n, 0, 1, AphiTab, real(S(phiMax)));
    // //n = 1000000;
    // trapRes  = trapezoidalMethod(n, phiMin, phiMax, complexAphiTab);

    // cout << filonRes << "\n" << trapRes << "  " << (S(phiMax)/derS(phiMin));
    
    
    
    
    
    // complex<double> result = filonRes;

    // file << "     I(A^2e^ih(phi))        (I(...)(n))-(I(...)(n-1))      N        P" <<"\n"; 
    
    // n = 111;

    // for (int j = n; j <= 1111 ; j += 20){
    
    //     complex<double> tmp = result;
        
    //     result = filonQuadrature(j, a, b, AphiTab, real(S(phiMax)));
    //     double absDiff = abs(abs(result)-abs(tmp)); 
    //     file << setprecision(10) << setfill (' ') << setw(20) <<  "\n" << result <<"        " << absDiff << "             " << j;
    // }

    // n = 1111;
    // for (int j = 1; j <= 1000 ; j ++){
    
    //     complex<double> tmp = result;
    //     p = j ;
    //     result = filonQuadrature(n, a, b, AphiTab, real(S(phiMax)));
         
    //     file << setprecision(10) << setfill (' ') << setw(20) <<  "\n" << result <<"                                 " << n << "             " << j;
    // }
    plt::show();
    file.close();
    return 0;

}