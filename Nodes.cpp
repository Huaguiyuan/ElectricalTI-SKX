#include "Nodes.hpp"
#include <armadillo>

using namespace arma;
using namespace std;

MagneticNode::MagneticNode(int new_index, vec new_Location, vec new_S)
{
    Location = new_Location;
    Spin = new_S;
    Temperature = 0.0;
    Index = new_index;
    ElectronSiteIndex = -1; // default value is for the non-electronic site case.
    /*Ku_zero = 1.3e6; // [J/m^3] 
    EasyAxis << 0.0 << 0.0 << 1.0;
    LatticeConstant = 3.9e-10;
    Beta = 2.5;
    Ms_zero = 1.0e7; //[A/m]*/
}

/*MagneticNode::MagneticNode()
{
    Index = -1;
    Spin = zeros<vec>(3);
    Temperature = 0.0;
    Location = zeros<vec>(3);
    Ku_zero = 1.3e6; // [J/m^3] 
    EasyAxis << 0.0 << 0.0 << 1.0;
    LatticeConstant = 3.9e-10;
    Beta = 2.5;
    Ms_zero = 1.0e7; //[A/m]
}*/
