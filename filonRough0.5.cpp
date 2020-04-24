# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>

using namespace std;

const complex<double> i(0,1);
const double pi = acos(-1);


complex<double> filonQuadrature(int n, double a, double b, complex<double> *f(int n , double x[]), double omega){


    double h = (b-a) / (n-1.);
    double *x;
    complex<double> *ftab;

    double theta = omega*h;
    double sint = sin(theta), cost = cos(theta);

    double alpha = (2. * pow(theta, 2.) - 2. * sint * sint + sin(2.*theta) * theta) /
                                            (2. * pow(theta, 3.));
    
    double beta  = (2. * theta * (1. + cost * cost ) - 2. * sin(2.*theta)) /
                                            (2. * pow(theta, 3.));
    
    double gamma = (4. * (-theta * cost) + sint) / 
                            pow(theta, 3.);

    x = new double [2*n+1];

    for (int j = 0; j <= 2*n; j++)
        x[j] = (double) a + (double) h*j; 

    ftab = f(2*n,x); 
    
    //cout << ftab[1] << endl;

    complex<double> sigma1;

    for (int j = 0; j<=n; j++)
        sigma1 += ftab[2*j] * exp(i*omega*x[2*j]);

    complex<double> sigma2;

    for (int j = 1; j<=n; j++)
        sigma1 += ftab[2*j-1] * exp(i*omega*x[2*j-1]);



    complex<double> result = h * (i * alpha * (ftab[0] * exp(i*omega*x[0]) - ftab[2*n] * exp(i*omega*x[2*n])) 
                                   
                                    + beta  * (sigma1 - 0.5 * (ftab[2*n] * exp(i*omega*x[2*n]) + ftab[0] * exp(i*omega*x[0])))
                                   
                                    + gamma  * sigma2

                                  );

    delete [] ftab;
    delete [] x;

    return result;

}


complex<double> *secondIntegrand (int n, double x[]){
  
  complex<double> *fx;
  int j;

  fx = new complex<double> [n+1];

  for ( j = 0; j <= n; j++ )
  {
    fx[j] = x[j] * x[j];
  }
  return fx;
}

double trapezoidalMethod(int n, double a, double b, complex<double> *f(int n , double x[]), double omega){
  
      
  double h = (double) (b-a) / n;
  double *x;
  complex<double> *ftab;

 
  x = new double [n+1];

    for (int j = 0; j <= n; j++)
        x[j] = (double) a + (double) h*j; 

    ftab = f(n,x); 

               
     for (int j = 0 ; j < n+1 ; j ++)                              
        ftab[j] *= cos(omega*x[j]) ;
    
  double result = 0;

    for (int j = 1; j < n ; j++)
      result += (double) real(ftab[j])*h;
  
  result += h/2.*(real(ftab[0])+real(ftab[n]));

  delete [] ftab;
  delete [] x;

  return result;

}

int main(){

complex<double> *secondIntegrand(int n, double x[]);

complex<double> result = filonQuadrature(11, 0, 2.*pi, secondIntegrand, 1111);

double result2 = trapezoidalMethod(10000000, 0, 2.*pi, secondIntegrand, 1111);     
cout << setw(24) << real(result) << endl;
cout << setw(24) << result2;

 return 0 ;
}