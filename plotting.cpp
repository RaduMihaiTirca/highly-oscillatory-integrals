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
   
   plt::loglog(relErrorCosOmega10.at(0).second,relErrorCosOmega10.at(1).second);
   plt::loglog(relErrorCosOmega10.at(0).second,relErrorCosOmega10.at(2).second);
   plt::save("output/relErrorCosOmega10.png");
   plt::save("/home/minastirith/Desktop/Licenta/latex/c4/relErrorCosOmega10.png");
   plt::show();
   
   plt::loglog(relErrorCosOmega1000.at(0).second,relErrorCosOmega1000.at(1).second);
   plt::loglog(relErrorCosOmega1000.at(0).second,relErrorCosOmega1000.at(2).second);
   plt::save("output/relErrorCosOmega1000.png");
   plt::save("/home/minastirith/Desktop/Licenta/latex/c4/relErrorCosOmega1000.png");
   plt::show();

   plt::loglog(relErrorPolyOmega10.at(0).second,relErrorPolyOmega10.at(1).second);
   plt::loglog(relErrorPolyOmega10.at(0).second,relErrorPolyOmega10.at(2).second);
   plt::save("output/relErrorPolyOmega10.png");
   plt::save("/home/minastirith/Desktop/Licenta/latex/c4/relErrorPolyOmega10.png");
   plt::show();

   plt::loglog(relErrorPolyOmega1000.at(0).second,relErrorPolyOmega1000.at(1).second);
   plt::loglog(relErrorPolyOmega1000.at(0).second,relErrorPolyOmega1000.at(2).second);
   plt::save("output/relErrorPolyOmega1000.png");
   plt::save("/home/minastirith/Desktop/Licenta/latex/c4/relErrorPolyOmega1000.png");
   plt::show();
   

   return 0;
}
