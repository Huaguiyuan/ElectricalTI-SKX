/* 
 * File:   ElectronSystem.hpp
 * Author: jason
 *
 * Created on October 14, 2014, 10:06 PM
 */

#ifndef ELECTRONSYSTEM_HPP
#define	ELECTRONSYSTEM_HPP
#include "ElectronSite.hpp"
#include "SpinSystem.hpp"
#include "OpenBoundary.hpp"
#include <iostream>
#include "discpp.h"


using namespace std;

class SpinSystem;
class ElectronSystem {
public:
    std::vector<ElectronSite> ListOfSites;
    std::vector<OpenBoundary> ListOfOpenBoundaries;
    cx_mat Hamiltonian;
    cx_mat GR;
    int TotalMatrixSize;
    int NumSite;
    ElectronSystem(const char* inputFilename, const char* boundaryFilename, 
                   const char* inputBoundaryShiftFilename, double t, double CouplingCutoff);
    void ReadInGeometry(const char* filename);
    void CreateNeighbourList(double CutoffRange, double t);
    void ReadInOpenBoundaries(const char* filename);
    void ReadInOpenBoundaryVirtialShift(const char* filename);
    void GenerateHamiltonian(void);
    void GenerateBoundaryHamiltonians(double CouplingT, double CouplingCutoff);
    void PrintBoundaryList(void);
    void PrintNeighbourList(void);
    void CalculateGR(double energy);
    void RenewGR(double energy);
    //void CalculateGreenNLinearResponse(double energy);
    cx_mat Transmission(double Energy);
    cx_mat ThermalAverageTransmission(double Temperature, double Ef);
    void PlotCurrentMap(int SizeX, int SizeY, int PlotX, int PlotY);
    void PlotSpinXCurrentMap(int SizeX, int SizeY, int PlotX, int PlotY);
    void PlotSpinYCurrentMap(int SizeX, int SizeY, int PlotX, int PlotY);
    void PlotSpinZCurrentMap(int SizeX, int SizeY, int PlotX, int PlotY);
    cx_mat OnSiteCorelationFunctionGn(int SiteI, int SiteJ, double Ef);
    cx_mat OnSiteSpectralFunctionFromBoundary(int I, int J, OpenBoundary &Source);
    //cx_mat SpectralBlockFromBoundary(OpenBoundary &Source);
    cx_mat GnBlockFromBoundary(OpenBoundary &Source, double Energy, double Temperature);
    cx_mat dGnBlockFromBoundary(OpenBoundary &Source, double Ef);
    void InjectionDOS(OpenBoundary &Source, double &Up, double &Down);
    void CalculateTerminalSpinCurrent(OpenBoundary &Source, double Ef, double &I, 
                                       double &Sx, double &Sy, double &Sz);
    void OutputSpinTextureProFit(const char* filename);
    void OutputSpinCurrentMapProFit(const char* filename, int Xstart, int Xend, int Ystart, int Yend);
    vec OnSiteSpin(int SiteIndex, double Ef);
    void UpdateHamiltonian(SpinSystem &SpinTexture, double JH);
    void OutputElectronSpinMapProFit(const char* filename, double Ef);
    void CalculateCurrentDistribution(double Ef);
    
private: 
    cx_mat ObtainGR_AB(OpenBoundary A, OpenBoundary B);
    double dFdE(double energy, double Ef, double temperature); 
    double Fermi_eV(double Energy, double Temperature, double Ef)
    {
        return 1.0/(1.0+exp((Energy-Ef)/8.6173324e-5/Temperature));
    }
};


#endif	/* ELECTRONSYSTEM_HPP */

