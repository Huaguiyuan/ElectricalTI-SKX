/* 
 * File:   LLG_Equation.hpp
 * Author: RAF_Jason
 *
 * Created on 2013年10月4日, 下午1:53
 */

#ifndef LLG_EQUATION_HPP
#define	LLG_EQUATION_HPP

#include <armadillo>
#include "NeighborList.hpp"
#include "InputOutput.hpp"
using namespace arma;
void UpdateExternalFieldAC(vec*** ExtenralFieldToUpdate, int SizeX, int SizeY, int SizeZ, double rate, vec h_AC, vec h0, 
                          double time, double Radius, double CutOff, double CenterX, double CenterY);
void RungeKuttaEvolve_PeriodicBoundary(vec** Plate, vec** ExternalField, NeighborRecord** Neighbors, 
                        double omega, vec h_AC, vec h0, double Radius, double CutOff, double Center1X, double Center1Y, double Center2X, double Center2Y,
                        double J, double D, double time, double TimeStep, double alpha, int SizeX, int SizeY, EnforcedSpin* EnforcedList, int N_Enforced);
void EnforceOpenBoundaries(vec** Plate, int SizeX, int SizeY);
void EnforceThePlate(EnforcedSpin* EnforcedList, int N_Enforced, vec*** Plate);
double TotalEnergy(vec*** Plate, double J, double D, int SizeX, int SizeY, int SizeZ, NeighborRecord*** NeighborList, 
                  vec*** ExternalField, double*** EnergyDistribution, double*** Heisenberg, double*** DM, double*** External);
void CalculateEffectiveField(double J, double D, vec*** PlateInput, NeighborRecord*** NeighborInput, vec*** ExternalFieldInput, vec*** EffectiveFieldOutput, int SizeX, int SizeY, int SizeZ);
void InitializeDotSpinWaveExcitation(vec** Plate, int x, int y, double Sx, double Sy, double Sz);
void EnforceExternalFieldOnSite(int x, int y, double Hy, vec** ExternalField);
double OnSiteEffectiveFieldZ(int i, int j, vec** Plate, NeighborRecord** NeighborInput, vec** ExternalFieldInput, double J, double D,
                           double &Heisenberg, double &DM, double &External, double &TotalXY, double &DMXY);

#endif	/* LLG_EQUATION_HPP */




