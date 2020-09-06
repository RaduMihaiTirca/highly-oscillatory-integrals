#include<vector>
#include<iostream>
#include<cmath>
#include<algorithm>

using namespace std;

const double INF = 1.e100;
vector<pair<double, double> > table;
double interpolate(double x) {
    // Assumes that "table" is sorted by .first
    // Check if x is out of bound
    vector<pair<double, double> >::iterator it, it2;
    if (x > table.back().first) return INF;
    
    if (x < table[0].first) return -INF;
    
    // INFINITY is defined in math.h in the glibc implementation
    it = lower_bound(table.begin(), table.end(), make_pair(x, -INF));
    // Corner case
    if (it == table.begin()) return it->second;
    it2 = it;
    --it2;
    return it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first);
}
int main() {
    table.push_back(make_pair(5., 15.));
    table.push_back(make_pair(7., 18.));
    table.push_back(make_pair(10., 22.));
    // If you are not sure if table is sorted:
    sort(table.begin(), table.end());
    printf("%f\n", interpolate(8.));
    printf("%f\n", interpolate(10.));
    printf("%f\n", interpolate(2.));
}