/* 
 * File:   OpenBoundaries.hpp
 * Author: jason
 *
 * Created on October 16, 2014, 2:50 PM
 */

#ifndef OPENBOUNDARIES_HPP
#define	OPENBOUNDARIES_HPP
#include "ElectronSite.hpp"

using namespace std;
class OpenBoundary 
{
public:
    std::vector<ElectronSite> ListOfBoundarySites;
    cx_mat F00;
    cx_mat F01; // This sub-matrix couples the boundary sites to the lead.
    cx_mat SelfEnergy;
    cx_mat Gamma;
    cx_mat HugeGamma;
    cx_mat SpectralFromThisBoundary;
    int BoundaryIndex;
    int TotalBoundaryMatrixSize;
    double ChemicalPotential;
    vec VirtualBoundaryShift;
    std::vector<int> StartingRowInBoundaryMatrix;
    void ConstructF00F01(double CouplingT, double CouplingCutoff);
    bool BelongsToThisBoundary(int ElectronSiteIndex);
    void GetSelfEnergy(double energy);
    void AddSelfEnergy(cx_mat &OpenHamiltonian);
    void GetGamma();
    void GetHugeGamma(int TotalSize);
    cx_mat GetSpectralFromHere(cx_mat GR);
    cx_mat SurfaceGreen;
    void GetSurfaceGreenFunction(double energy, double Eta);
    
};

#endif	/* OPENBOUNDARIES_HPP */

