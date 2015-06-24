#include "OpenBoundary.hpp"
#include <armadillo>
#define PIPI 3.1415926535897932384626
using namespace arma;
typedef complex<double> Complex;

void OpenBoundary::GetSurfaceGreenFunction(double energy, double Eta)
{
    cx_mat F00Rs_new = F00;
    cx_mat F00R_new = F00;
    cx_mat F01R_new = F01;
    cx_mat F10R_new = trans(F01);
    int Size = F00.n_rows;
    Complex MaxT = 1.0;
    cx_mat F00Rs_old, F00R_old, F01R_old, F10R_old, G_temp;
    MaxT = 1.0;
    while (abs(MaxT) > 1.0e-10)
    {
        F00Rs_old = F00Rs_new;
        F00R_old = F00R_new;
        F01R_old = F01R_new;
        F10R_old = F10R_new;
        G_temp = inv(  eye<cx_mat>(Size, Size)*Complex(energy, Eta) - F00R_old);
        F01R_new = F01R_old*G_temp*F01R_old;
        F10R_new = F10R_old*G_temp*F10R_old;
        F00R_new = F00R_old + F01R_old*G_temp*F10R_old + F10R_old*G_temp*F01R_old;
        F00Rs_new = F00Rs_old + F01R_old*G_temp*F10R_old;
        MaxT = max(max(F01R_new));     
    }
    SurfaceGreen = inv(  eye<cx_mat>(Size, Size)*Complex(energy, Eta) - F00Rs_new);    
}

////////////////

void OpenBoundary::ConstructF00F01(double CouplingT, double CouplingCutoff)
{
    int N_matrix = this->TotalBoundaryMatrixSize;
    this->F00.set_size(N_matrix, N_matrix);
    this->F01.set_size(N_matrix, N_matrix);
    F00.zeros();
    F01.zeros();
    cx_mat T;
    cx_mat I;
    int StartingRowCount = 0;
    for (int i=0; i<ListOfBoundarySites.size(); i++)
    {
        // for each site, first take the on-site block
        this->StartingRowInBoundaryMatrix.push_back(StartingRowCount);
        int BlockSize = ListOfBoundarySites[i].OnSiteBlock.n_rows;
        StartingRowCount += BlockSize;
        I.eye(BlockSize, BlockSize);
        F00.submat(i*BlockSize, i*BlockSize, (i+1)*BlockSize-1, (i+1)*BlockSize-1)
                = ListOfBoundarySites[i].OnSiteBlock;
        // now put the coupling blocks in;
        int NeighbourIndex;
        for (int j=0; j<ListOfBoundarySites.size(); j++)
        {
            if (i == j)
                continue;
            NeighbourIndex = ListOfBoundarySites[i].IsNeighbour(ListOfBoundarySites[j]);
            if (NeighbourIndex >= 0)
            {
                T = ListOfBoundarySites[i].ListOfOutwardsCouplingBlocks[NeighbourIndex];
                F00.submat(i*BlockSize, j*BlockSize, (i+1)*BlockSize-1, (j+1)*BlockSize-1) = T;
            }
        }
       /* I.eye(BlockSize, BlockSize);
        F01.submat(i*BlockSize, i*BlockSize, (i+1)*BlockSize-1, (i+1)*BlockSize-1) = I*CouplingT;*/
        //The above two lines were the cheating code neglecting the detailed atomic structure at the boundary.
        //Now we formally formulate the F01 matrix.
        
        
    }
    
    for (int i=0; i<this->ListOfBoundarySites.size(); i++)
    {
        for (int j=0; j<this->ListOfBoundarySites.size(); j++)
        {
            int I, J;
            I = this->StartingRowInBoundaryMatrix[i];
            J = this->StartingRowInBoundaryMatrix[j];
            int BlockSizeRow = ListOfBoundarySites[i].OnSiteBlock.n_rows;
            int BlockSizeCol = ListOfBoundarySites[j].OnSiteBlock.n_rows;
            T.resize(BlockSizeRow, BlockSizeCol);
            T.eye(BlockSizeRow, BlockSizeCol);
            T = T*CouplingT;
            vec r(3);
            r = (ListOfBoundarySites[i].Location) - (ListOfBoundarySites[j].Location+this->VirtualBoundaryShift);
            double distanceSquared = r(0)*r(0)+r(1)*r(1)+r(2)*r(2);
            if (distanceSquared < CouplingCutoff*CouplingCutoff)
            {
                F01.submat(I, J, I+BlockSizeRow-1, J+BlockSizeCol-1) = T;
            }
        }
    }
}

///////////////
void OpenBoundary::GetSelfEnergy(double energy)
{
    this->GetSurfaceGreenFunction(energy, 1.0e-8);
    this->SelfEnergy = (this->F01)*SurfaceGreen*trans((this->F01));
    //printf("Self energy generated for boundary %d.\n", this->BoundaryIndex);
    //SelfEnergy.print("SelfE =");
}

////////////////
void OpenBoundary::AddSelfEnergy(cx_mat  &OpenHamiltonian)
{
    int NumBoundarySites = this->ListOfBoundarySites.size();
    cx_mat SelfEnergyBlock;
    int SizePerSite = 2;
    //ofstream before("Before.txt");
    //OpenHamiltonian.print(before);
    for (int i=0; i<NumBoundarySites; i++)
    {
        for (int j=0; j<NumBoundarySites; j++)
        {
            int I = this->ListOfBoundarySites[i].SiteIndex;
            int J = this->ListOfBoundarySites[j].SiteIndex;
            SelfEnergyBlock = this->SelfEnergy.submat(i*SizePerSite, j*SizePerSite, (i+1)*SizePerSite-1, (j+1)*SizePerSite-1);
            OpenHamiltonian.submat(I*SizePerSite, J*SizePerSite, (I+1)*SizePerSite-1, (J+1)*SizePerSite-1)
                = 
                OpenHamiltonian.submat(I*SizePerSite, J*SizePerSite, (I+1)*SizePerSite-1, (J+1)*SizePerSite-1)
                +
                    SelfEnergyBlock;
                    
        }
    }
    //ofstream after("After.txt");
    //OpenHamiltonian.print(after);
}
//////
bool OpenBoundary::BelongsToThisBoundary(int ElectronSiteIndex)
{
    for (int i=0; i<this->ListOfBoundarySites.size(); i++)
    {
        if (ElectronSiteIndex == this->ListOfBoundarySites[i].SiteIndex)
            return true;
    }
    return false;
}
//////
void OpenBoundary::GetGamma()
{
    Complex I(0.0, 1.0);
    this->Gamma = I*(this->SelfEnergy - trans(this->SelfEnergy));
}
//////
void OpenBoundary::GetHugeGamma(int TotalSize)
{
    this->GetGamma();
    this->HugeGamma.zeros(TotalSize, TotalSize);
    this->AddSelfEnergy(HugeGamma);
    HugeGamma = trans(HugeGamma);
    HugeGamma = -HugeGamma;
    this->AddSelfEnergy(HugeGamma);
    HugeGamma = HugeGamma*Complex(0.0, 1.0);
}
///////
cx_mat OpenBoundary::GetSpectralFromHere(cx_mat GR)
{
    this->GetHugeGamma(GR.n_rows);
    return GR*HugeGamma*trans(GR);
}
/////
void OpenBoundary::CalculateBandStructure(const char* filename)
{
    FILE *fp;
    fp = fopen(filename, "w");
    for (double ka = 0.0; ka <= PIPI; ka+= PIPI/100.0)
    {
        cx_mat Hk = F00 + F01*exp(Complex(0.0, 1.0)*ka) + trans(F01)*exp(-Complex(0.0, 1.0)*ka);
        vec eigval = eig_sym(Hk);
        fprintf(fp, "% lf\t", ka);
        for (int i=0; i<eigval.n_rows; i++)
        {
            fprintf(fp, "% le\t", eigval(i));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
