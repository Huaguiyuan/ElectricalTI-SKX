/* 
 * File:   Nodes.hpp
 * Author: jason
 *
 * Created on July 14, 2014, 9:25 AM
 */

#ifndef NODES_HPP
#define	NODES_HPP
#include <armadillo>

using namespace std;
using namespace arma;

class MagneticNode
{
public:
    vec Spin;   // dimensionless magnetization, normalized to Ms_Zero;
    int Index;  // The integer that labels the node.
   /* double Hex0_bulk;   // [A/m], zero temperature exchange field.
    double LatticeConstant; // [m]
    double Ku_zero; */        // easy axis anisotropy energy density  at zero temperature. [eV/m^3]
    double Temperature;     // [K]
    bool Pinned;
   /* double Ms_zero; // zero-T magnetization. [A/m]
    double Beta; // The temperature dependency of the magnetization.
    vec EasyAxis;         // a unit vector along the easy axis.
    vec ExchangeStifness; // unit eV/m*/
    std::vector<int> ListOfExchangeNeighbours;
    //std::vector<int> ListOfStrayFieldNeighbours;
    std::vector<vec> ListOfOurwardsVectors;
   // std::vector<mat> ListOfTensors;
    vec Location; // This gives the location of the spin.
    MagneticNode(int new_index, vec new_Location, vec new_S);
    MagneticNode();
};

#endif	/* NODES_HPP */



