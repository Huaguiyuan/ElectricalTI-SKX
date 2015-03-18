/* 
 * File:   SpinSystem.hpp
 * Author: jason
 *
 * Created on July 14, 2014, 12:20 PM
 */

#ifndef SPINSYSTEM_HPP
#define	SPINSYSTEM_HPP
#include "Nodes.hpp"
#include "randomc.h"
#include "ElectronSystem.hpp"
//#include "voro++.hh"

//using namespace voro;
class ElectronSystem;
class SpinSystem 
{
public:
    double CurrentTime;
    double J;
    double D;
    int NumSite;
    //container* VoroContainer;
    //mat** Tensors;
    vector<MagneticNode> NodeList;
    vector<vec> EffectiveField;
    vector<vec> ExternalField;
    vector<vec> RandomField;
    vector<vec> BackgroundField;
    vector<int> ListOfTorqueSiteIndecies;
    CRandomMersenne* RanGen;  //The random number sequence for the stochastic field calculation.
    SpinSystem(const char* filename, double J_initial, double D_initial, double alpha);
    ~SpinSystem();
    void GenerateNeighborList(double ExchangeCutoff, double StrayCutoff, 
        bool PeriodicX, bool PeriodicY, bool PeriodicZ, vec ax, vec ay, vec az);
    void PrintSpinTextureToTextFile(const char* filename);
    void CalculateStrayFieldTensor(void);
    void FormatTheSystem(vec S_input);
    mat Tensor(vec PQ);
    void CalculateRandomField(double TimeStep);
    void Evolve(double TimeStep, double Time, bool AddTorque);
    void Initialize(void);
    void CalculateEffectiveField(bool AddTorque); 
    void UpdateExternalField(double Time);
    void CreateWindow(void);
    void UpdateWindowDisplay(void);
    void SetTemperature(double newTemperature);
    void ReadInElectronSiteIndex(const char* filename);
    void SetBackgroundField(vec BackgroundField);
    void CalculateTorque(ElectronSystem &Electrons, double Ef, double J_Hunds);
    void OutputTextureToTextFile(const char* filename);
    void OutputTorqueFieldToProFitTextFile(const char* filename);
    void OutputEffectiveToProFitTextFile(const char* filename);
    void OutputSpinTextureGIF(double Xmin, double Xmax, double Ymin, double Ymax, const char* title);
    vec SkyrmionLocation(void);
};

#endif	/* SPINSYSTEM_HPP */

