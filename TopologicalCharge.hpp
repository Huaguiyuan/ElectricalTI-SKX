/* 
 * File:   TopologicalCharge.hpp
 * Author: jason
 *
 * Created on October 10, 2013, 10:50 PM
 */

#ifndef TOPOLOGICALCHARGE_HPP
#define	TOPOLOGICALCHARGE_HPP
#include <armadillo>
#include "NeighborList.hpp"

using namespace arma;
double CalculateTopologicalCharge(vec*** PlateInput, double*** TopologicalChargeDensityOutput, NeighborRecord*** Neighbors, 
                              int SizeX, int SizeY, int SizeZ);
double CalculateContinuumTopologicalCharge(vec** PlateInput, double** TopologicalChargeDensityOutput, NeighborRecord** Neighbors, 
                              int SizeX, int SizeY);


#endif	/* TOPOLOGICALCHARGE_HPP */

