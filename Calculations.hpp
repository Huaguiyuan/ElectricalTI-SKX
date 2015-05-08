/* 
 * File:   Calculations.hpp
 * Author: jason
 *
 * Created on February 4, 2015, 12:23 PM
 */

#ifndef CALCULATIONS_HPP
#define	CALCULATIONS_HPP

#include <cstdlib>
#include <armadillo>
#include "InputOutput.hpp"
#include "NeighborList.hpp"
#include "LLG_Equation.hpp"
#include "TopologicalCharge.hpp"
#include "FiniteTemperature.hpp"
#include <stdio.h>
#include <unistd.h>
#include "discpp.h"
#include <complex>
#include "SpinSystem.hpp"
//#include "voro++.hh"
#include "ElectronSystem.hpp"

#define PIPI 3.1415962535897932384626


using namespace std;
using namespace arma;
//haha 
//haha 
//using namespace voro;
typedef complex<double> Complex;
int SkyrmionCalculation(void);
void CalculateCurrentDistribution(double Ef, double Miu1, double Miu2, double Miu4, ElectronSystem &World);
void CalculateHallEffect(double Temperature, double Ef, double Miu1, 
                         double &Miu2, double &Miu4, double &totalI, double &totalSz,
                         double &UpDos, double &DownDos);
void CalculateInjectionBandStructure(double t);
double TopologicalCharge(SpinSystem &input);

#endif	/* CALCULATIONS_HPP */

