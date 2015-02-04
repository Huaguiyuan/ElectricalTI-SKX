#include "ElectronSite.hpp"
#include <complex>
typedef complex<double> Complex;

ElectronSite::ElectronSite(vec newSpin, vec newLocation, double JH, double newPhi)
{
    cx_mat Pauli_X(2,2);
    cx_mat Pauli_Y(2,2);
    cx_mat Pauli_Z(2,2);
    Pauli_X.zeros();
    Pauli_Y.zeros();
    Pauli_Z.zeros();
    Pauli_X(0,1) = 1.0;
    Pauli_X(1,0) = 1.0;
    Pauli_Y(0,1) = Complex(0.0, -1.0);
    Pauli_Y(1,0) = Complex(0.0, 1.0);
    Pauli_Z(0.0) = 1.0;
    Pauli_Z(1,1) = -1.0;
    
    this->OnSiteBlock =  JH*newSpin(0)*Pauli_X
                        +JH*newSpin(1)*Pauli_Y
                        +JH*newSpin(2)*Pauli_Z;
    this->Spin = newSpin;
    this->phi = newPhi;
    this->Location = newLocation;
    this->IsBoundary = false;
}

/*int ElectronSite::GetNeighbourIndex(ElectronSite A)
{
    for (int i; i<this->ListOfTightBindingNeighbors.size(); i++)
    {
        if (this->ListOfTightBindingNeighbors[i] == A.SiteIndex)
            return i;
    }
    return -1;
}*/

int ElectronSite::IsNeighbour(ElectronSite A)
{
    for (int i=0; i<this->ListOfTightBindingNeighbors.size(); i++)
    {
        if (ListOfTightBindingNeighbors[i] == A.SiteIndex)
            return i;
    }
    return -1;
}
