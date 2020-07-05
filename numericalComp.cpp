# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>
# include <vector>
# include "matplotlibcpp.h"
# include "csvIO.h"

using namespace std;
namespace plt = matplotlibcpp;

const complex<double> iComplex(0,1);
const double pi = acos(-1);


complex<double> filonQuadrature(int n, double a, double b, complex<double> *f(int n , double x[]), double omega){


    double h = (b-a) / (2*n);
    double *x;
    complex<double> *ftab;

    double theta = omega*h;
    double sint = sin(theta), cost = cos(theta);
    double alpha, beta, gamma;

    if ( 6.0 * fabs ( theta ) <= 1.0 ){

        alpha = 2.0 * pow ( theta, 3 )   /   45.0 
              - 2.0 * pow ( theta, 5 )   /  315.0 
              + 2.0 * pow ( theta, 7 )   / 4725.0;
    
        beta  = 2.0                      /     3.0 
              + 2.0 * pow ( theta, 2 )   /    15.0 
              - 4.0 * pow ( theta, 4 )   /   105.0 
              + 2.0 * pow ( theta, 6 )   /   567.0 
              - 4.0 * pow ( theta, 8 )   / 22275.0;

        gamma = 4.0                      /     3.0 
              - 2.0 * pow ( theta, 2 )   /    15.0 
              +       pow ( theta, 4 )   /   210.0 
              -       pow ( theta, 6 )   / 11340.0;
    }
    else{

        alpha = ( pow ( theta, 2 ) + theta * sint * cost 
        - 2.0 * sint * sint ) / pow ( theta, 3 );

        beta = ( 2.0 * theta + 2.0 * theta * cost * cost
        - 4.0 * sint * cost ) / pow ( theta, 3 );

        gamma = 4.0 * ( sint - theta * cost ) / pow ( theta, 3 );
    }

    x = new double [2*n+1];

    for (int j = 0; j <= 2*n; j++)
        x[j] = (double) a + (double) h*j; 

    ftab = f(2*n,x); 
    
    //cout << ftab[1] << endl;

    complex<double> sigma1(0,0);

    for (int j = 0; j<=n; j++)
        sigma1 += ftab[2*j] * exp(iComplex*omega*x[2*j]);

    complex<double> sigma2(0,0);

    for (int j = 1; j<=n; j++)
        sigma2 += ftab[2*j-1] * exp(iComplex*omega*x[2*j-1]);



    complex<double> result = h * (iComplex * alpha * (ftab[0] * exp(iComplex*omega*x[0]) - ftab[2*n] * exp(iComplex*omega*x[2*n])) 
                                   
                                    + beta  * (sigma1 - 0.5 * (ftab[2*n] * exp(iComplex*omega*x[2*n]) + ftab[0] * exp(iComplex*omega*x[0])))
                                   
                                    + gamma  * sigma2

                                  );

    delete [] ftab;
    delete [] x;

    return result;

}

complex<double> trapezoidalMethod(int n, double a, double b, complex<double> *f(int n , double x[]), double omega){
  
      
  double h = (double) (b-a) / n;
  double *x;
  complex<double> *ftab;

 
  x = new double [n+1];

  for (int j = 0; j <= n; j++)
    x[j] = (double) a + (double) h*j; 

  ftab = f(n,x); 

               
  for (int j = 0 ; j < n+1 ; j ++)                              
      ftab[j] *= exp(iComplex*omega*x[j]) ;
    
  complex<double> result = 0;

  for (int j = 1; j < n ; j++)
    result += (double) real(ftab[j])*h;
  
  result += h/2.*(ftab[0]+ftab[n]);

  delete [] ftab;
  delete [] x;

  return result;

}

complex<double> *toIntegratePoly (int n, double x[]){
  
  complex<double> *fx;
  int j;

  fx = new complex<double> [n+1];

  for ( j = 0; j <= n; j++ ){
    fx[j] = x[j] * x[j] + 3*x[j] + pow(x[j],3);
  }

  return fx;
}

complex<double> *toIntegrateCos (int n, double x[]){
  
  complex<double> *fx;
  int j;

  fx = new complex<double> [n+1];

  for ( j = 0; j <= n; j++ ){
    fx[j] = cos(x[j]);
  }

  return fx;
}


int main(){

  int n = 9;
  double a = 0;
  double b = 100;
  double omega = 10;
  double h = (double) (b-a)/n;
  
  vector<double> x(n), y(n), z(n);
  
  complex<double> exactCos = (exp(iComplex*a*omega)*(iComplex*omega*cos(a) 
                         + sin(a)) - iComplex*exp(iComplex*b*omega)*(omega*cos(b) 
                         - iComplex*sin(b))) 
                         / (-1 + pow(omega,2)); 

  
  complex<double> exactPoly = (exp(iComplex*a*omega)*(6.0 - omega*(iComplex*2.0 + 3*omega + a*(iComplex*6.0 + omega*(2 + 3*a - iComplex*(3 + a + pow(a,2))*omega)))) 
                         +  exp(iComplex*b*omega)*(-6.0 + omega*(iComplex*2.0 + 3*omega + b*(iComplex*6.0 + omega*(2 + 3*b - iComplex*(3 + b + pow(b,2))*omega))))) 
                         /  pow(omega,4);


  for (int i = 0; i < n; i++){
    complex<double> trap = trapezoidalMethod(pow(10,i), a, b, toIntegrateCos, omega);
    complex<double> filon = filonQuadrature(pow(10,i)+1, a, b, toIntegrateCos, omega);
    x.at(i) = pow(10,i);
    y.at(i) = fabs(real((trap-exactCos)/exactCos)*100);
    z.at(i) = fabs(real((filon-exactCos)/exactCos)*100);
  }
  
  vector<pair<string, vector<double>>> toWrite = {{"x",x} , {"y",y}, {"z",z}};
  writeToCsv("output/relErrorCosOmega10.csv", toWrite);
  writeToCsv("/home/minastirith/Desktop/Licenta/latex/c4/relErrorCosOmega10.csv", toWrite);

  for (int i = 0; i < n; i++){
    complex<double> trap = trapezoidalMethod(pow(10,i), a, b, toIntegratePoly, omega);
    complex<double> filon = filonQuadrature(pow(10,i)+1, a, b, toIntegratePoly, omega);
    x.at(i) = pow(10,i);
    y.at(i) = fabs(real((trap-exactPoly)/exactPoly)*100);
    z.at(i) = fabs(real((filon-exactPoly)/exactPoly)*100);
  }
  
  toWrite = {{"x",x} , {"y",y}, {"z",z}};
  writeToCsv("output/relErrorPolyOmega10.csv", toWrite);
  writeToCsv("/home/minastirith/Desktop/Licenta/latex/c4/relErrorPolyOmega10.csv", toWrite);
  

  omega = 1000;

  exactCos = (exp(iComplex*a*omega)*(iComplex*omega*cos(a) 
                         + sin(a)) - iComplex*exp(iComplex*b*omega)*(omega*cos(b) 
                         - iComplex*sin(b))) 
                         / (-1 + pow(omega,2));

  exactPoly = (exp(iComplex*a*omega)*(6.0 - omega*(iComplex*2.0 + 3*omega + a*(iComplex*6.0 + omega*(2 + 3*a - iComplex*(3 + a + pow(a,2))*omega)))) 
                         +  exp(iComplex*b*omega)*(-6.0 + omega*(iComplex*2.0 + 3*omega + b*(iComplex*6.0 + omega*(2 + 3*b - iComplex*(3 + b + pow(b,2))*omega))))) 
                         /  pow(omega,4);

  for (int i = 0; i < n; i++){
    complex<double> trap = trapezoidalMethod(pow(10,i), a, b, toIntegrateCos, omega);
    complex<double> filon = filonQuadrature(pow(10,i)+1, a, b, toIntegrateCos, omega);
    x.at(i) = pow(10,i);
    y.at(i) = fabs(real((trap-exactCos)/exactCos)*100);
    z.at(i) = fabs(real((filon-exactCos)/exactCos)*100);
  }

  toWrite = {{"x",x} , {"y",y}, {"z",z}};
  writeToCsv("output/relErrorCosOmega1000.csv", toWrite);
  writeToCsv("/home/minastirith/Desktop/Licenta/latex/c4/relErrorCosOmega1000.csv", toWrite);
  

  for (int i = 0; i < n; i++){
    complex<double> trap = trapezoidalMethod(pow(10,i), a, b, toIntegratePoly, omega);
    complex<double> filon = filonQuadrature(pow(10,i)+1, a, b, toIntegratePoly, omega);
    x.at(i) = pow(10,i);
    y.at(i) = fabs(real((trap-exactPoly)/exactPoly)*100);
    z.at(i) = fabs(real((filon-exactPoly)/exactPoly)*100);
  }
  
  toWrite = {{"x",x} , {"y",y}, {"z",z}};
  writeToCsv("output/relErrorPolyOmega1000.csv", toWrite);
  writeToCsv("/home/minastirith/Desktop/Licenta/latex/c4/relErrorPolyOmega1000.csv", toWrite);
  
 return 0 ;
}