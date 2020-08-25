

# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <complex>
# include <vector>


using namespace std;


const complex<double> iComplex(0,1);
const double pi      = acos(-1);
//const double C       = 137.036;
//const double k       = 0.05/C;
//const double tau     = 5*C*(2.*pi/0.05);
const double omega   = 1;
const double A0      = 0.5;
const double phiMax  = 0;
const double phiMin  = 100;
      //double p       = 1;


double Aphi(double x){

  return pow(x,2)+2*pow(x,4);

}
double *AphiTab(int n, double x[]){ //good!

  double *fx;
  
  fx = new double[n+1];
  
  for (int i = 0; i <= n; i++){

    fx[i] = Aphi(x[i]);

  }
  return fx;
}
complex<double> *complexAphiTab(int n, double x[]){ //good!

  complex<double> *fx;
  
  fx = new complex<double>[n+1];
  
  for (int i = 0; i <= n; i++){

    fx[i] = (0,0) + Aphi(x[i]);

  }
  return fx;
}




complex<double> hasPhi(double phi){ //good!
  
  return omega*phi+A0*sin(omega*phi);
}

complex<double> *hasPhiTab (int n, double phi[]){ //good!

  complex<double> *fx;

  fx = new complex<double>[n];

  for(int i = 0; i < n; i++)
    fx[i] =  hasPhi(phi[i]) ;

  return fx;

}

double derHasPhi(double x){

  //double tmp = A0 * exp( -1 * ( (x * x) / (tau * tau) ) ) * sin(k * x);

    return (0,0) + omega + A0*omega*cos(omega*x);

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
    double yh= (double) (yb-ya)/n;

    //cout << "lim = " << a << " " << ya << " " << b << " "  << yb << "\n"; 
    
    vector<double> Yt(n+1), Y(n+1), Xy(n+1), X(n+1);


  for (int i = 0 ; i <= n; i ++){
      X.at(i)  = a + i*h; 
      Y.at(i)  = ya + i*yh;
      Yt.at(i) = real(hasPhi(X.at(i)));  //cout << Yt.at(i) << " ";
  }
     //cout << Y.at(0) << " " << Y.at(n) <<"\n";
  for(int i = 0 ; i <= n; i ++){
    double x = Y.at(i); cout << i <<" x: " << x << " ";
    double y = interpolate(X, Yt, x, true); cout << "y: " << y << " ";
    fx[i]=Aphi(y)/derHasPhi(y); cout << "fx: " << fx[i] << "\n";
    //Xy.push_back(y);
    
  }
//   cout << "\n\n\n";
//   for (int i = 0; i <= n; i++){
//   //  cout << X.at(i) << " " << Aphi(Xy.at(i))*Aphi(Xy.at(i)) << " " << Xy.at(i) << " " << Yt.at(i) << " "; 
//     fx[i] = (Aphi(Xy.at(i))*Aphi(Xy.at(i)))/derHasPhi(Xy.at(i));
//  //   cout << derHasPhi(Xy.at(i)) << " ";
//    // cout << fx[i] << "\n";
//   }

//   cout << "\n\n\n";
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

    for (int i = 0; i <= 2*n; i++)
        v[i] = (double) a + i*h; 

    x = new double [2*n+1];

    x[0] = real(hasPhi(phiMin));

    for (int i = 0; i < 2*n; i++){
        double tmp = x[i];
        x[i+1] = tmp + (omega/derHasPhi(tmp))*h;
    }
    
    cout << x[0] << " " << x[1] << "  " << h << " " ;

    ftab = f(2*n,x); 

    for (int i = 0; i <= 2*n; i++){
      ftab[i] *= ftab[i];
      ftab[i] /= derHasPhi(x[i]);
    }
    complex<double> sigma1(0,0);

    for (int i = 0; i<=n; i++)
        sigma1 += ftab[2*i] * exp(iComplex*omega*v[2*i]);

    complex<double> sigma2(0,0);

    for (int i = 1; i<=n; i++)
        sigma2 += ftab[2*i-1] * exp(iComplex*omega*v[2*i-1]);



    complex<double> result = h * omega * (iComplex * alpha * (ftab[0] - ftab[2*n] * exp(iComplex* omega)) 
                                   
                                    + beta  * (sigma1 - 0.5 * (ftab[0] + ftab[2*n] * exp(iComplex* omega)))
                                   
                                    + gamma * sigma2

                                  );

    delete [] ftab;
    delete [] x;

    return result;

}

complex<double> simpleFilonQuadrature(int n, double a, double b, complex<double> *f(int n , double x[]), double omega){


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

    for (int i = 0; i <= 2*n; i++)
        x[i] = (double) a + (double) h*i; 
        
    ftab = f(2*n,x); 
    

    //cout << x[0] << " " << x[2*n] << "\n";

    //cout << ftab[1] << endl;

    complex<double> sigma1(0,0);

    for (int i = 0; i<=n; i++)
        sigma1 += ftab[2*i] * exp(iComplex*omega*x[2*i]);

    complex<double> sigma2(0,0);

    for (int i = 1; i<=n; i++)
        sigma2 += ftab[2*i-1] * exp(iComplex*omega*x[2*i-1]);



    complex<double> result = h * (iComplex * alpha * (ftab[0] * exp(iComplex*omega*x[0]) - ftab[2*n] * exp(iComplex*omega*x[2*n])) 
                                   
                                    + beta  * (sigma1 - 0.5 * (ftab[2*n] * exp(iComplex*omega*x[2*n]) + ftab[0] * exp(iComplex*omega*x[0])))
                                   
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
        ftab[i] *= ftab[i];
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




int main(){

    int n = 20;
    double a = phiMin;
    double b = phiMax;
  //   double h = (double) (b-a)/n;

    double ya=real(hasPhi(a));
    double yb=real(hasPhi(b));
  //   double yh= (double) (yb-ya)/n;
    
  //   vector<double> Yt(n), Y(n), Xy(n), X(n);


  // for (int i = 0 ; i < n; i ++){
  //     X.at(i)  = a + i*h; 
  //     Y.at(i)  = real(hasPhi(a)) + i*yh;
  //     Yt.at(i) = real(hasPhi(X.at(i)));
  // }
    
  // for(double x : Y){
  //   double y = interpolate(Y, Yt, x, true);
  //   Xy.push_back(y);
  // }

    
  
    complex<double> filonRes, trapRes;
    filonRes = simpleFilonQuadrature(n, ya, yb, transVal, 1);
    n = 10000000;
    trapRes  = trapezoidalMethod(n, phiMin, phiMax, complexAphiTab);

    cout << filonRes << "\n" << trapRes << "  " ;//<< (S(phiMax)/derS(phiMin));
    // complex<double> result = filonRes;

    // file << "     I(A^2e^ih(phi))        (I(...)(n))-(I(...)(n-1))      N        P" <<"\n"; 
    
    // n = 111;

    // for (int i = n; i <= 1111 ; i += 20){
    
    //     complex<double> tmp = result;
        
    //     result = filonQuadrature(i, a, b, AphiTab, real(S(phiMax)));
    //     double absDiff = abs(abs(result)-abs(tmp)); 
    //     file << setprecision(10) << setfill (' ') << setw(20) <<  "\n" << result <<"        " << absDiff << "             " << i;
    // }

    // n = 1111;
    // for (int i = 1; i <= 1000 ; i ++){
    
    //     complex<double> tmp = result;
    //     p = i ;
    //     result = filonQuadrature(n, a, b, AphiTab, real(S(phiMax)));
         
    //     file << setprecision(10) << setfill (' ') << setw(20) <<  "\n" << result <<"                                 " << n << "             " << i;
    // }
    
    return 0;

}