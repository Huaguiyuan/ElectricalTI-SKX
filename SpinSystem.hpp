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
//#include "voro++.hh"

//using namespace voro;

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
    vector<vec> RamdomField;
    CRandomMersenne* RanGen;  //The random number sequence for the stochastic field calculation.
    SpinSystem(const char* filename, double J_initial, double D_initial);
    ~SpinSystem();
    void GenerateNeighborList(double ExchangeCutoff, double StrayCutoff, 
        bool PeriodicX, bool PeriodicY, bool PeriodicZ);
    void PrintSpinTextureToTextFile(const char* filename);
    void CalculateStrayFieldTensor(void);
    void FormatTheSystem(vec S_input);
    mat Tensor(vec PQ);
    void CalculateRandomField(void);
    void Evolve(double deltaT, double* TemperatureList, vec* ExternalField);
    void Initialize(void);
    void CalculateEffectiveField(void);   
};

#endif	/* SPINSYSTEM_HPP */

