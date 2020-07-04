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

complex<double> *toIntegrate (int n, double x[]){
  
  complex<double> *fx;
  int j;

  fx = new complex<double> [n+1];

  for ( j = 0; j <= n; j++ )
  {
    fx[j] = x[j] * x[j] + 3*x[j] + pow(x[j],3);
  }
  return fx;
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

int main(){

  int n = 10000;
  double a = 0;
  double b = 1000;
  double omega = 100;
  double h = (double) (b-a)/n;
  double *xTab;
  complex<double> *tab;

  xTab = new double [n+1];

  for(int i = 0; i < n ; i++)
    xTab[i] = (double) a + i*h;
  
  tab = toIntegrate(n,xTab);  

  for(int i = 0; i < n; i ++)
    tab[i] *= exp(iComplex*xTab[i]*omega);

  vector<double> x(n), y(n), z(n);
  
  for (int i = 0; i < n; i++){
    x.at(i) = a + i*h;
    y.at(i) = real(tab[i]) ; 
  }
  vector<pair<string, vector<double>>> toWrite = {{"x",x} , {"y",y}};
  writeToCsv("data.csv", toWrite);
  vector<pair<string, vector<double>>> data = readFromCsv("data.csv");
  plt::plot(data.at(0).second,data.at(1).second);
  plt::show();

  complex<double> *toIntegrate(int n, double x[]);

  // complex<double> result = (exp(iComplex*a*omega)*(iComplex*omega*cos(a) 
  //                        + sin(a)) - iComplex*exp(iComplex*b*omega)*(omega*cos(b) 
  //                        - iComplex*sin(b))) 
  //                        / (-1 + pow(omega,2)); //filonQuadrature(1111, a, b, toIntegrate, omega);


  complex<double> result = (exp(iComplex*a*omega)*(6.0 - omega*(iComplex*2.0 + 3*omega + a*(iComplex*6.0 + omega*(2 + 3*a - iComplex*(3 + a + pow(a,2))*omega)))) 
                         +  exp(iComplex*b*omega)*(-6.0 + omega*(iComplex*2.0 + 3*omega + b*(iComplex*6.0 + omega*(2 + 3*b - iComplex*(3 + b + pow(b,2))*omega))))) 
                         /  pow(omega,4);

  complex<double> result2 = trapezoidalMethod(100, a, b, toIntegrate, omega);
  cout << setw(24) << real(result) << endl;
  cout << setw(24) << real(result2);

  // for(int i = 0; i < 2; i ++)
  //   tab[i] = trapezoidalMethod(pow(10,i), a, b, toIntegrate, omega);

  // for (int i = 1; i < 2; i++){
  //   x.at(i) = pow(10,i);
  //   y.at(i) = fabs(real(tab[i] - result)/real(result))*100 ;
  //   z.at(i) = fabs(real(filonQuadrature(pow(10,i)+1, a, b, toIntegrate, omega)-result)/real(result))*100;
  // }
  // plt::loglog(x,y);
  // plt::loglog(x,z);
  // plt::show();




 return 0 ;
}