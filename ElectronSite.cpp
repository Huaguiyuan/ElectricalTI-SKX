#include "ElectronSite.hpp"
#include "ElectronSystem.hpp"
#include <complex>
typedef complex<double> Complex;

/*ElectronSite::ElectronSite(vec newSpin, vec newLocation, double JH, double newPhi)
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
}*/


ElectronSite::ElectronSite(vec newSpin, vec newLocation, double JH, double newPhi)
{
    cx_mat Pauli_X(2,2);
    cx_mat Pauli_Y(2,2);
    cx_mat Pauli_Z(2,2);
    cx_mat Ho(2,2);
    Ho.zeros();
    Pauli_X.zeros();
    Pauli_Y.zeros();
    Pauli_Z.zeros();
    
    //the following variables occur in electronsite, electron system, openboundary
    double deltaX = 15e-10; //interatomic distance along x direction
    double deltaY = deltaX; //intteratomic distance along y direction
    double hv = 3.2955e-10*0.015; //hbar*VF in electron volt, Vf =  0.5*  0.5e6m/s
    double gamma = 5e-10*1;
    
    Pauli_X(0,1) = 1.0;
    Pauli_X(1,0) = 1.0;
    Pauli_Y(0,1) = Complex(0.0, -1.0);
    Pauli_Y(1,0) = Complex(0.0, 1.0);
    Pauli_Z(0.0) = 1.0;
    Pauli_Z(1,1) = -1.0;
    
   
    Ho = (hv*4*gamma/(deltaX*deltaX))*Pauli_Z;
    
    
    
    
    this->OnSiteBlock =  Ho 
                        +JH*newSpin(0)*Pauli_X
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
