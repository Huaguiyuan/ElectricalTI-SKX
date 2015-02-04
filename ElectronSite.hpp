/* 
 * File:   ElectronSite.hpp
 * Author: jason
 *
 * Created on October 14, 2014, 10:07 PM
 */

#ifndef ELECTRONSITE_HPP
#define	ELECTRONSITE_HPP
#include <armadillo>
using namespace std;
using namespace arma;
class ElectronSite {
public:
    vec Location;
    vec Current;
    vec SpinXCurrent;
    vec SpinYCurrent;
    vec SpinZCurrent;
    double phi;
    std::vector<int> ListOfTightBindingNeighbors;
    cx_mat OnSiteBlock;
    std::vector<cx_mat> ListOfOutwardsCouplingBlocks;
    vec Spin;
    int SiteIndex;
    int StartingRowNumber;
    int SiteBlockSize;
    bool IsBoundary;
   // int GetNeighbourIndex(ElectronSite A);
    ElectronSite(vec newSpin, vec newLocation, double JH, double phi);
    int IsNeighbour(ElectronSite A);
};


#endif	/* ELECTRONSITE_HPP */

