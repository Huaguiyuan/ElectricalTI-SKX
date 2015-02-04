


#include "LLG_Equation.hpp"
#include "randomc.h"
#include "InputOutput.hpp"
#include "NeighborList.hpp"
using namespace std;
using namespace arma;

void GetRandomField(vec*** RandomField, double temperature, double TimeStep, double alpha, CRandomMersenne &RanGen, 
                     int SizeX, int SizeY, int SizeZ)
{
    vec L(3);
    double TempRand;
    double Lambda = alpha/(1.+alpha*alpha)*temperature;
    for (int x=0; x<SizeX; x++)
    {
        for (int y=0; y<SizeY; y++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                for (int p=0; p<3; p++)
                {
                    TempRand = RanGen.Random() - 0.5;  // From -0.5 to +0.5;
                    L(p) = sqrt(24.*Lambda*TimeStep)*TempRand;
                }
                RandomField[x][y][k] = L;
            }
        }
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void Evolve_PeriodicBoundary_FiniteTemperature(CRandomMersenne &RanGen, double temperature, vec*** Plate, vec*** ExternalField, 
                        NeighborRecord*** Neighbors, double rate, vec h_AC, vec h0, double Radius, double CutOff, double CenterX, 
                        double CenterY, double J, double D, double time, 
                        double TimeStep, double alpha, int SizeX, int SizeY, int SizeZ)
{
    // This function evolves Plate AND ExternalField from time to time+TimeStep
    // Random force is added to the effective field at each time step.
    vec A(3);
    vec BL(3);
    vec A_tilde(3);
    vec BL_tilde(3);
    vec S(3);
    vec L(3);
    vec Heff(3);
    vec*** EffectiveField;
    vec*** S_tilde;
    vec*** PlateUpdate;
    vec*** RandomField;
    PlateUpdate = Allocate_3D_Vector_Plate(SizeX, SizeY, SizeZ);
    S_tilde = Allocate_3D_Vector_Plate(SizeX, SizeY, SizeZ);
    EffectiveField = Allocate_3D_Vector_Plate(SizeX, SizeY, SizeZ);
    RandomField = Allocate_3D_Vector_Plate(SizeX, SizeY, SizeZ);
    // First, generate the random field;
    GetRandomField(RandomField, temperature, TimeStep, alpha, RanGen, SizeX, SizeY, SizeZ);
    // Now calculate S_tilde and add the contribution at "time"
    UpdateExternalFieldAC(ExternalField, SizeX, SizeY, SizeZ, rate, h_AC, h0, time, Radius, CutOff, CenterX, CenterY);
    //EnforceOpenBoundaries(Plate, SizeX, SizeY);
    CalculateEffectiveField(J, D, Plate, Neighbors, ExternalField, EffectiveField, SizeX, SizeY, SizeZ);
    //printf("haha\n");
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            S_tilde[i][j][0] = zeros(3);
            S_tilde[i][j][SizeZ-1] = zeros(3);
            for (int k=1; k<SizeZ-1; k++)
            {
                S = Plate[i][j][k];
                L = RandomField[i][j][k];
                Heff = EffectiveField[i][j][k];
                A = -1.0/(1+alpha*alpha)*cross(S, Heff)-alpha/(1.0+alpha*alpha)*(S*(dot(S, Heff)) - Heff);
                BL = -1.0/(1+alpha*alpha)*cross(S, L)-alpha/(1.0+alpha*alpha)*(S*(dot(S, L)) - L);
                S_tilde[i][j][k] = S + A*TimeStep + BL;
                PlateUpdate[i][j][k] = 0.5*A*TimeStep + 0.5*BL;
            }
        } 
    }
    //Now calculate the PlateUpdate and ExternalField to time+TimeStep;
    UpdateExternalFieldAC(ExternalField, SizeX, SizeY, SizeZ, rate, h_AC, h0, time+TimeStep, Radius, CutOff, CenterX, CenterY);
//    printf("haha\n");
    //EnforceOpenBoundaries(Plate, SizeX, SizeY);
    CalculateEffectiveField(J, D, S_tilde, Neighbors, ExternalField, EffectiveField, SizeX, SizeY, SizeZ);
 //   printf("haha\n");
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                S = S_tilde[i][j][k];
                L = RandomField[i][j][k];
                Heff = EffectiveField[i][j][k];
                vec PlateUpdatePrevious = PlateUpdate[i][j][k];
                vec PlatePrevious = Plate[i][j][k];
                A = -1.0/(1+alpha*alpha)*cross(S, Heff)-alpha/(1.0+alpha*alpha)*(S*(dot(S, Heff)) - Heff);
                BL = -1.0/(1+alpha*alpha)*cross(S, L)-alpha/(1.0+alpha*alpha)*(S*(dot(S, L)) - L);
                PlateUpdate[i][j][k] = PlateUpdatePrevious + 0.5*A*TimeStep + 0.5*BL;
                Plate[i][j][k] = PlatePrevious + PlateUpdate[i][j][k];
            }
        } 
    }
    Delete_3D_Plate(S_tilde, SizeX, SizeY);
    Delete_3D_Plate(EffectiveField, SizeX, SizeY);
    Delete_3D_Plate(PlateUpdate, SizeX, SizeY);
    Delete_3D_Plate(RandomField, SizeX, SizeY);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SpinNormalization(vec*** Plate, int SizeX, int SizeY, int SizeZ)
{
    double length;
    vec temp(3);
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                temp = Plate[i][j][k];
                length = sqrt(temp(0)*temp(0)+temp(1)*temp(1)+temp(2)*temp(2));
                Plate[i][j][k] = Plate[i][j][k]/length;
            }
        }
    }
}