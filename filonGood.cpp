
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <complex>
# include <vector>

using namespace std;


const complex<double> iComplex(0,1);
const double pi      = acos(-1);
const double omega   = 1;
const double omega1  = 2;
const double A0      = 0.1;
const double phiMax  = 10;
const double phiMin  = 0;



double Aphi(double x){

  return pow(x,2)+2*pow(x,4);

}

complex<double> *complexAphiTab(int n, double x[]){ 

  complex<double> *fx;
  
  fx = new complex<double>[n+1];
  
  for (int i = 0; i <= n; i++){

    fx[i] = (0,0) + Aphi(x[i]);

  }
  return fx;
}




complex<double> hasPhi(double phi){ 
  
  return (0,0) + omega*phi+A0*sin(omega1*phi);
}

complex<double> *hasPhiTab (int n, double phi[]){ 

  complex<double> *fx;

  fx = new complex<double>[n];

  for(int i = 0; i < n; i++)
    fx[i] =  hasPhi(phi[i]) ;

  return fx;

}

double derHasPhi(double x){

  

    return (0,0) + omega + A0*omega1*cos(omega1*x);

}

double interpolate( vector<double> &xData, vector<double> &yData, double x, bool extrapolate )
{
   int size = xData.size();

   int i = 0;                                                                  // find left end of interval for interpolation
   if ( x >= xData[size - 2] )                                                 // special case: beyond right end
   {
      i = size - 2;
   }
   else
   {
      while ( x > xData[i+1] ) i++;
   }
   double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];      // points on either side (unless beyond ends)
   if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
   {
      if ( x < xL ) yR = yL;
      if ( x > xR ) yL = yR;
   }

   double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

   return yL + dydx * ( x - xL );                                              // linear interpolation
}

complex<double> *transVal(int n, double x[]){
  complex<double> *fx;
  
  fx = new complex<double>[n+1];

    double a = phiMin;
    double b = phiMax;
    double h = (double) (b-a)/n; 

   

    double ya=real(hasPhi(a));
    double yb=real(hasPhi(b));
    double yh= (yb-ya)/n;

   vector<double> Yt(n+1), Y(n+1), Xy(n+1), X(n+1);

    
  for (int i = 0 ; i <= n; i ++){
      X.at(i)  = a + i*h; 
      Y.at(i)  = ya + i*yh;
      Yt.at(i) = real(hasPhi(X.at(i)));
      
  }
     
    fx[0]=Aphi(Y.at(0))*exp(iComplex*Y.at(0))/derHasPhi(Y.at(0));
    fx[n]=Aphi(Y.at(n))*exp(iComplex*Y.at(n))/derHasPhi(Y.at(n)); 

  for(int i = 1 ; i < n; i ++){
    double val = Y.at(i); 
    double y = interpolate(X,Yt,val,true); 
    x[i] = y;
    fx[i]=Aphi(y)*exp(iComplex*y)/derHasPhi(y); 
  }

  return fx;
}


complex<double> simpleFilonQuadrature(int n, double a, double b, complex<double> *f(int n , double x[]), double omega){


    double h = (b-a) / (2*n);
    double *x;
    complex<double> *ftab;

    double theta = h;
    double sint = sin(theta), cost = cos(theta);
    double alpha, beta, gamma;
    
    if ( 6.0 * fabs ( theta ) >= 1.0 ){

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

   
        
    ftab = f(2*n,x); 
    


   

    complex<double> sigma1(0,0);

    for (int i = 0; i<=n; i++)
        sigma1 += ftab[2*i] ;

    complex<double> sigma2(0,0);

    for (int i = 1; i<=n; i++)
        sigma2 += ftab[2*i-1] ;


    complex<double> result = h * (iComplex * alpha * (ftab[0]  - ftab[2*n]) 
                                   
                                    + beta  * (sigma1 - 0.5 * (ftab[2*n]  + ftab[0] ))
                                   
                                    + gamma  * sigma2

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

    for (int i = 0; i <= n; i++)
        x[i] = (double) a + (double) h*i; 

    ftab = f(n,x); 
    hasTab = hasPhiTab(n+1,x) ;
               
     for (int i = 0 ; i < n+1 ; i ++){                              
        ftab[i] *= exp(iComplex*hasTab[i]) ;
     }
  complex<double> result(0,0) ;

    for (int i = 1; i < n ; i++)
      result += ftab[i]*h;
  
  result += h/2.*(ftab[0] +ftab[n]);

  delete [] ftab;
  delete [] x;

  return result;

}

complex<double> trapezoidalMethodDebug(int n, double a, double b, complex<double> *f(int n , double x[])){
  
      
  double h = (double) (b-a) / n;
  double *x;
  complex<double> *ftab;

  x   = new double[n+1]; 

    for (int i = 0; i <= n; i++)
        x[i] = (double) a + (double) h*i; 

    ftab = f(n,x); 


  complex<double> result(0,0) ;

    for (int i = 1; i < n ; i++)
      result += ftab[i]*h;
  
  result += h/2.*(ftab[0] +ftab[n]);

  delete [] ftab;
  delete [] x;

  return result;

}



int main(){

    int n = 101;
    double a = phiMin;
    double b = phiMax;
  

    double ya=real(hasPhi(a));
    double yb=real(hasPhi(b));
  

    
  
    complex<double> filonRes, trapRes, trapResDebug;
    
    
    filonRes = simpleFilonQuadrature(n, ya, yb, transVal, 1);
    
    n = 100000;
    
    trapRes  = trapezoidalMethod(n, phiMin, phiMax, complexAphiTab);
    trapResDebug = trapezoidalMethodDebug(n, ya, yb, transVal);
    
    cout    << trapRes << "\n" << filonRes << "\n" << trapResDebug ;



    return 0;

}