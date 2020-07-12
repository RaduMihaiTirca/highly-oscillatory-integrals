# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>
# include <vector>
# include "matplotlibcpp.h"
# include "csvIO.h"

using namespace std;
namespace plt = matplotlibcpp;

int main()
{
   
   vector<pair<string, vector<double>>> relErrorPolyOmega10 = readFromCsv("output/relErrorPolyOmega10.csv");
   vector<pair<string, vector<double>>> relErrorPolyOmega1000 = readFromCsv("output/relErrorPolyOmega1000.csv");
   
   vector<pair<string, vector<double>>> relErrorCosOmega10 = readFromCsv("output/relErrorCosOmega10.csv");
   vector<pair<string, vector<double>>> relErrorCosOmega1000 = readFromCsv("output/relErrorCosOmega1000.csv");
   
   vector<double> x = relErrorCosOmega10.at(0).second;
   vector<double> y = relErrorCosOmega10.at(1).second;
   vector<double> z = relErrorCosOmega10.at(2).second;
   // cout << x.size() << "\n" << y.size() << "\n" << z.size();
   //plt::subplots_adjust(sharex='col', sharey='row')
   plt::subplot(2,2,1);
   plt::loglog(x, y, "o-", {{"label", "$Q^T$"}});
   plt::loglog(x, z, "o-", {{"label", "$Q^F$"}});
   plt::xlabel("$N$");
   plt::ylabel("$\\delta$");
   plt::legend("upper left",{1.05, 1}); 
   //plt::loglog(x,z,".-");
   
   // plt::save("output/relErrorCosOmega10.png");
   // plt::save("/home/minastirith/Desktop/Licenta/latex/c4/relErrorCosOmega10.png");
   // plt::show();

   plt::subplot(2,2,2);
   plt::loglog(relErrorCosOmega1000.at(0).second,relErrorCosOmega1000.at(1).second,"o-");
   plt::loglog(relErrorCosOmega1000.at(0).second,relErrorCosOmega1000.at(2).second,"o-");
   plt::xlabel("$N$");
   plt::ylabel("$\\delta$");
   
   // plt::save("output/relErrorCosOmega1000.png");
   // plt::save("/home/minastirith/Desktop/Licenta/latex/c4/relErrorCosOmega1000.png");
   // plt::show();
   plt::subplot(2,2,3);
   plt::loglog(relErrorPolyOmega10.at(0).second,relErrorPolyOmega10.at(1).second,"o-");
   plt::loglog(relErrorPolyOmega10.at(0).second,relErrorPolyOmega10.at(2).second,"o-");
   plt::xlabel("$N$");
   plt::ylabel("$\\delta$");
   
   // plt::save("output/relErrorPolyOmega10.png");
   // plt::save("/home/minastirith/Desktop/Licenta/latex/c4/relErrorPolyOmega10.png");
   // plt::show();
   plt::subplot(2,2,4);
   plt::loglog(relErrorPolyOmega1000.at(0).second,relErrorPolyOmega1000.at(1).second,"o-");
   plt::loglog(relErrorPolyOmega1000.at(0).second,relErrorPolyOmega1000.at(2).second,"o-");
   plt::xlabel("$N$");
   plt::ylabel("$\\delta$");
   
   // plt::save("output/relErrorPolyOmega1000.png");
   // plt::save("/home/minastirith/Desktop/Licenta/latex/c4/relErrorPolyOmega1000.png");
  

   plt::show();
   

   return 0;
}
