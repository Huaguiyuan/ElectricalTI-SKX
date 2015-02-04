/* 
 * File:   FiniteTemperature.hpp
 * Author: jason
 *
 * Created on October 23, 2013, 9:01 PM
 */

#ifndef FINITETEMPERATURE_HPP
#define	FINITETEMPERATURE_HPP
#include "randomc.h"
#include <armadillo>

using namespace arma;

void GetRandomField(vec** RandomField, double temperature, double TimeStep, double alpha, CRandomMersenne &RanGen, 
                     int SizeX, int SizeY);
void Evolve_PeriodicBoundary_FiniteTemperature(CRandomMersenne &RanGen, double temperature, vec*** Plate, vec*** ExternalField, 
                        NeighborRecord*** Neighbors, double rate, vec h_AC, vec h0, double Radius, double CutOff, double CenterX, 
                        double CenterY, double J, double D, double time, 
                        double TimeStep, double alpha, int SizeX, int SizeY, int SizeZ);
void SpinNormalization(vec*** Plate, int SizeX, int SizeY, int SizeZ);


#endif	/* FINITETEMPERATURE_HPP */

