#include <armadillo>
#include "NeighborList.hpp"
#include "InputOutput.hpp"
#include "LLG_Equation.hpp"
#include <math.h>
#define PIPI 3.1415926535897932384626
using namespace arma;

void CalculateEffectiveField(double J, double D, vec*** PlateInput, 
        NeighborRecord*** NeighborInput, vec*** ExternalFieldInput, 
        vec*** EffectiveFieldOutput, int SizeX, int SizeY, int SizeZ)
{
    vec SpinRight(3);
    vec SpinLeft(3);
    vec SpinUp(3);
    vec SpinDown(3);
    vec SpinTop(3);
    vec SpinBottom(3);
    //vec ExtraTerm(3);  // For eazy-plane;
    vec S(3);
    vec X(3);
    vec Y(3);
    vec Z(3);
    vec Ext(3);
    X.zeros(3);
    Y.zeros(3);
    Z.zeros(3);
    //ExtraTerm.zeros(); // This is for the easy-plane;
    X(0) = 1.0;
    Y(1) = 1.0;
    Z(2) = 1.0;
    //S = PlateInput[i][j];
    //ExtraTerm(2) = S(2);
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                S = PlateInput[i][j][k];
                Ext = ExternalFieldInput[i][j][k];
                //ExtraTerm(2) = S(2);   // This is for easy-plane
                SpinRight = PlateInput[NeighborInput[i][j][k].Right[0]][NeighborInput[i][j][k].Right[1]][NeighborInput[i][j][k].Right[2]];
                SpinLeft = PlateInput[NeighborInput[i][j][k].Left[0]][NeighborInput[i][j][k].Left[1]][NeighborInput[i][j][k].Left[2]];
                SpinUp = PlateInput[NeighborInput[i][j][k].Up[0]][NeighborInput[i][j][k].Up[1]][NeighborInput[i][j][k].Up[2]];
                SpinDown = PlateInput[NeighborInput[i][j][k].Down[0]][NeighborInput[i][j][k].Down[1]][NeighborInput[i][j][k].Down[2]];
                SpinTop = PlateInput[NeighborInput[i][j][k].Top[0]][NeighborInput[i][j][k].Top[1]][NeighborInput[i][j][k].Top[2]];
                SpinBottom = PlateInput[NeighborInput[i][j][k].Bottom[0]][NeighborInput[i][j][k].Bottom[1]][NeighborInput[i][j][k].Bottom[2]];
                EffectiveFieldOutput[i][j][k] = J*(SpinRight+SpinLeft+SpinUp+SpinDown+SpinTop+SpinBottom) - D*cross((SpinRight-SpinLeft), X) - 
                    D*cross((SpinUp-SpinDown),Y) - D*cross((SpinTop-SpinBottom),Z) + Ext; //- 0.7*0.0036*ExtraTerm; 
               // printf("%d %d %d\n",i,j,k);
                
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void UpdateExternalFieldAC(vec*** ExtenralFieldToUpdate, int SizeX, int SizeY, int SizeZ, double rate, vec h_AC, vec h0, 
                          double time, double Radius, double CutOff, double CenterX, double CenterY)
{   
    /* Now start the cylinder electrode method. */
    //Radius is the radius of the linear field. outside the linear field it is 1/r decay. 
    //We use the cut-off distance at 3.0*R
    /*
          { h0*r/R
     H = {
          { h0*R/r
     */
    /* This is for the cylindrical electrode (realistic) */
    
    vec h_AC_temp(3);
    h_AC_temp(2) = 0.0;
    double theta;
    double SineValue;
    double sign = 1.0;
    /*int iPeriod = (int) (time/2.0/PIPI*omega);
    double sign;
    if (iPeriod%2==0)
        sign = 1.0;
    else
        sign = -1.0;
    if (sin(omega*time)>=0.0)
        SineValue = -1.0;//-sin(omega*time);
    else
        SineValue = 0.0;*/
    //if (sin(omega*time)>0.)
    
    /*if (time<3500)//PIPI/omega)//(time < 2.0*PIPI/omega)
        SineValue = 2.0*omega/PIPI*time;//sin(omega*time);
    else
        SineValue = -2.0*omega/PIPI*time + 2.0*2.0*omega/PIPI*3500;//4.0*omega/PIPI*PIPI/omega;//SineValue = 1.0; //1.0;//0.0;
    if (SineValue < 0.0)
        SineValue = 0.0;
    SineValue = -SineValue;*/
    SineValue=1.0;
    double Magnitude=rate*time;
    //SineValue = 1.0;
    // This is for the impulse test;
    /*if (!(fabs(time - 0.0) < 1.0e-5))
    {
        h_AC.zeros();
    }*/
    // impulse test finished here.
    /*if (time > 300.)
        sign = 0.0;
    else
        sign = 1.0;*/
    if (sign>0)
    {
        for (int i=(int)(CenterX-(CutOff)); i<=(int)(CenterX+(CutOff)); i++)
        {
            for (int j=(int)(CenterY-(CutOff)); j<=(int)(CenterY+(CutOff)); j++)
            {
                for (int k=1; k<SizeZ-1; k++)
                {
                //double SinValue = sin(omega*time);
                    theta = atan2((double)((double)j-CenterY), (double)((double)i-CenterX));
                    double DistanceSquare = ((double)i-CenterX)*((double)i-CenterX)+((double)j-CenterY)*((double)j-CenterY);
                    if (DistanceSquare < (CutOff)*(CutOff))
                    {
                        if (DistanceSquare < Radius*Radius)
                        {
                            h_AC_temp(0) = /*sign**/Magnitude/(double)Radius*sqrt((double)DistanceSquare)*-sin(theta);
                            h_AC_temp(1) = /*sign**/Magnitude/(double)Radius*sqrt((double)DistanceSquare)*cos(theta);
                        }
                        else
                        {
                            h_AC_temp(0) = /*sign**/Magnitude*(double)Radius/sqrt((double)DistanceSquare)*-sin(theta);
                            h_AC_temp(1) = /*sign**/Magnitude*(double)Radius/sqrt((double)DistanceSquare)*cos(theta);
                        }
                        ExtenralFieldToUpdate[i][j][k] = h0 - (h_AC_temp*SineValue);/*+1.0*h_AC_temp*SineValue*///*cos(omega*time);
                    }
                    else
                        ExtenralFieldToUpdate[i][j][k] = h0;
                }
            }
        }
    }
    else
    {
        for(int i=0; i<SizeX; i++)
        {
            for( int j=0; j<SizeY; j++)
            {
                for (int k=1; k<SizeZ-1; k++)
                {
                    ExtenralFieldToUpdate[i][j][k] = h0;
                }
            }
        }
    }
     
    
    /* This is for the square sharp unrealistic field.
    int iPeriod = (int) (time/2.0/PIPI*omega);
    double sign;
    if (iPeriod%2==0)
        sign = -1.0;
    else
        sign = 1.0;
    for (int i=Center1X-Radius; i<=Center1X+Radius; i++)
    {
        for (int j=Center1Y-Radius; j<=Center1Y+Radius; j++)
        {
            if (  1.0)//(i-Center1X)*(i-Center1X)+(j-Center1Y)*(j-Center1Y) <= Radius*Radius   )
            {
                double sinValue = sin(omega*time);
                if (sinValue > 0.0)
                    ExtenralFieldToUpdate[i][j] = h0 + h_AC;
                else
                    ExtenralFieldToUpdate[i][j] = h0;
            }
                
        }
    }
    for (int i=Center2X-Radius; i<=Center2X+Radius; i++)
    {
        for (int j=Center2Y-Radius; j<=Center2Y+Radius; j++)
        {
            if ( 1.0)//(i-Center2X)*(i-Center2X)+(j-Center2Y)*(j-Center2Y) <= Radius*Radius   )
            {
                double sinValue = sin(omega*time);
                if (sinValue > 0.0)
                    ExtenralFieldToUpdate[i][j] = h0 +sign*h_AC;
                else
                    ExtenralFieldToUpdate[i][j] = h0;
            }
                
        }
    }
     */
}
////////////////////////////////////////////////////////////////////////////

void EnforceThePlate(EnforcedSpin* EnforcedList, int N_Enforced, vec*** Plate)
{
    int x,y,z;
    vec temp(3);
    for (int i=0; i<N_Enforced; i++)
    {
        x = EnforcedList[i].x;
        y = EnforcedList[i].y;
        z = EnforcedList[i].z;
        temp = (EnforcedList[i].spin);
        Plate[x][y][z] = temp;
    }
}
////////////////////////////////////////////////////////////////////////////
void RungeKuttaEvolve_PeriodicBoundary(vec*** Plate, vec*** ExternalField, NeighborRecord*** Neighbors, 
                        double omega, vec h_AC, vec h0, double Radius, double CutOff, double CenterX, double CenterY,
                        double J, double D, double time, double TimeStep, double alpha, int SizeX, int SizeY, int SizeZ, EnforcedSpin* EnforcedList, 
                        int N_Enforced)
{
    // This function evolves Plate AND ExternalField from time to time+TimeStep
    // The TextureToUpDate and ExternalField are at time.
    // First we need
    vec*** K1;
    vec*** K2;
    vec*** K3;
    vec*** K4;
    vec*** EffectiveField;
    vec*** PlatePerturbed;
    K1 = Allocate_3D_Vector_Plate(SizeX, SizeY, SizeZ);
    K2 = Allocate_3D_Vector_Plate(SizeX, SizeY, SizeZ);
    K3 = Allocate_3D_Vector_Plate(SizeX, SizeY, SizeZ);
    K4 = Allocate_3D_Vector_Plate(SizeX, SizeY, SizeZ);
    PlatePerturbed = Allocate_3D_Vector_Plate(SizeX, SizeY, SizeZ);
    EffectiveField = Allocate_3D_Vector_Plate(SizeX, SizeY, SizeZ);
    UpdateExternalFieldAC(ExternalField, SizeX, SizeY, SizeZ, omega, h_AC, h0, time, Radius, CutOff, CenterX, CenterY);
    EnforceThePlate(EnforcedList, N_Enforced, Plate);
    CalculateEffectiveField(J, D, Plate, Neighbors, ExternalField, EffectiveField, SizeX, SizeY, SizeZ);
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                K1[i][j][k] = TimeStep*(-1.0/(1+alpha*alpha)*cross(Plate[i][j][k], EffectiveField[i][j][k])-
                               alpha/(1.0+alpha*alpha)*(Plate[i][j][k]*(dot(Plate[i][j][k], EffectiveField[i][j][k])) - EffectiveField[i][j][k]));
            }
        }
    }
    UpdateExternalFieldAC(ExternalField, SizeX, SizeY, SizeZ, omega, h_AC, h0, time+TimeStep/2., Radius, CutOff, CenterX, CenterY);
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                PlatePerturbed[i][j][k] = Plate[i][j][k] + 0.5*K1[i][j][k];
            }
        }
    }
    EnforceThePlate(EnforcedList, N_Enforced, Plate);
    CalculateEffectiveField(J, D, PlatePerturbed, Neighbors, ExternalField, EffectiveField, SizeX, SizeY, SizeZ);
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                K2[i][j][k] = TimeStep*(-1.0/(1+alpha*alpha)*cross(PlatePerturbed[i][j][k], EffectiveField[i][j][k])-
                                  alpha/(1.0+alpha*alpha)*(PlatePerturbed[i][j][k]*(dot(PlatePerturbed[i][j][k], EffectiveField[i][j][k])) - EffectiveField[i][j][k]));
            }
        }
    }
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                PlatePerturbed[i][j][k] = Plate[i][j][k] + 0.5*K2[i][j][k];
            }
        }
    }
    EnforceThePlate(EnforcedList, N_Enforced, Plate);
    CalculateEffectiveField(J, D, PlatePerturbed, Neighbors, ExternalField, EffectiveField, SizeX, SizeY, SizeZ);
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                K3[i][j][k] = TimeStep*(-1.0/(1+alpha*alpha)*cross(PlatePerturbed[i][j][k], EffectiveField[i][j][k])-
                       alpha/(1.0+alpha*alpha)*(PlatePerturbed[i][j][k]*(dot(PlatePerturbed[i][j][k], EffectiveField[i][j][k])) - EffectiveField[i][j][k]));
            }
        }
    }
    UpdateExternalFieldAC(ExternalField, SizeX, SizeY, SizeZ, omega, h_AC, h0, time+TimeStep, Radius, CutOff, CenterX, CenterY);
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                PlatePerturbed[i][j][k] = Plate[i][j][k] + K3[i][j][k];
            }
        }
    }
    EnforceThePlate(EnforcedList, N_Enforced, Plate);
    CalculateEffectiveField(J, D, PlatePerturbed, Neighbors, ExternalField, EffectiveField, SizeX, SizeY, SizeZ);
    // By now, ExternalField is updated to time+TimeStep.
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                K4[i][j][k] = TimeStep*(-1.0/(1+alpha*alpha)*cross(PlatePerturbed[i][j][k], EffectiveField[i][j][k])-
                       alpha/(1.0+alpha*alpha)*(PlatePerturbed[i][j][k]*(dot(PlatePerturbed[i][j][k], EffectiveField[i][j][k])) - EffectiveField[i][j][k]));
            }
        }
    }
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                Plate[i][j][k] = Plate[i][j][k] + 1.0/6.0*(K1[i][j][k] + 2.0*K2[i][j][k] + 2.0*K3[i][j][k] + K4[i][j][k]);
            }
        }
    }
    EnforceThePlate(EnforcedList, N_Enforced, Plate);
    // By now, Plate is updated to time + TimeStep;
    Delete_3D_Plate(K1, SizeX, SizeY);
    Delete_3D_Plate(K2, SizeX, SizeY);
    Delete_3D_Plate(K3, SizeX, SizeY);
    Delete_3D_Plate(K4, SizeX, SizeY);
    Delete_3D_Plate(EffectiveField, SizeX, SizeY);
    Delete_3D_Plate(PlatePerturbed, SizeX, SizeY);
}
//////////////////////////////////////////////////////////////////////////////
void EnforceOpenBoundaries(vec*** Plate, int SizeX, int SizeY, int SizeZ)
{
    vec zero(3);
    zero.zeros();
    for (int k=1; k<SizeZ-1; k++)
    {
        for (int i=0; i<SizeX; i++)
        {
            Plate[0][i][k] = zero;
            Plate[SizeY-1][i][k] = zero;
        }
    
    
        for (int i=1; i<SizeY-1; i++)
        {
            Plate[i][0][k] = zero;
            Plate[i][SizeX-1][k] = zero;
        }
    }
}
/////////////////////////////////////////////////////////////////////////
double TotalEnergy(vec*** Plate, double J, double D, int SizeX, int SizeY, 
                   int SizeZ, NeighborRecord*** NeighborList, 
                   vec*** ExternalField, double*** EnergyDistribution, 
                   double*** Heisenberg, double*** DM, double*** External)
{
    double TotalEnergy=0.0;
    vec Left(3);
    vec Right(3);
    vec Up(3);
    vec Down(3);
    vec Top(3);
    vec Bottom(3);
    //vec FM_part = ExternalField[0][0][0];
    //vec AC_part = ExternalField[0][0][0];
    //E_AC = 0.0;
    //FM_part(0) = 0.0;
    //FM_part(1) = 0.0;
    Left.zeros();
    Right.zeros();
    Up.zeros();
    Down.zeros();
    Top.zeros();
    Bottom.zeros();
    Left(0) = -1.0;
    Right(0) = 1.0;
    Up(1) = 1.0;
    Down(1) = -1.0;
    Top(2) = 1.0;
    Bottom(2) = -1.0;
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                int xLeft, yLeft, zLeft;
                int xUp, yUp, zUp;
                int xRight, yRight, zRight;
                int xDown, yDown, zDown;
                int xTop, yTop, zTop;
                int xBottom, yBottom, zBottom;
                
                xLeft = NeighborList[i][j][k].Left[0];
                xRight = NeighborList[i][j][k].Right[0];
                xUp = NeighborList[i][j][k].Up[0];
                xDown = NeighborList[i][j][k].Down[0];
                xTop = NeighborList[i][j][k].Top[0];
                xBottom = NeighborList[i][j][k].Bottom[0];
                
                yLeft = NeighborList[i][j][k].Left[1];
                yRight = NeighborList[i][j][k].Right[1];
                yUp = NeighborList[i][j][k].Up[1];
                yDown = NeighborList[i][j][k].Down[1];
                yTop = NeighborList[i][j][k].Top[1];
                yBottom = NeighborList[i][j][k].Bottom[1];
                
                zLeft = NeighborList[i][j][k].Left[2];
                zRight = NeighborList[i][j][k].Right[2];
                zUp = NeighborList[i][j][k].Up[2];
                zDown = NeighborList[i][j][k].Down[2];
                zTop = NeighborList[i][j][k].Top[2];
                zBottom = NeighborList[i][j][k].Bottom[2];
                
                
                vec ThisSpin = Plate[i][j][k];
                vec LeftSpin = Plate[xLeft][yLeft][zLeft];
                vec RightSpin = Plate[xRight][yRight][zRight];
                vec UpSpin = Plate[xUp][yUp][zUp];
                vec DownSpin = Plate[xDown][yDown][zDown];
                vec TopSpin = Plate[xTop][yTop][zTop];
                vec BottomSpin = Plate[xBottom][yBottom][zBottom];
                
                Heisenberg[i][j][k] =  -J*dot(ThisSpin, LeftSpin)/2.0
                                    -J*dot(ThisSpin, RightSpin)/2.0
                                    -J*dot(ThisSpin, UpSpin)/2.0
                                    -J*dot(ThisSpin, DownSpin)/2.0
                                    -J*dot(ThisSpin, TopSpin)/2.0
                                    -J*dot(ThisSpin, BottomSpin)/2.0;
                vec temp = ExternalField[i][j][k];
                temp(0) = 0.0;
                temp(1) = 0.0;
                External[i][j][k] = -dot(temp, ThisSpin);//-dot(ExternalField[i][j][k], ThisSpin);
                //AC_part = ExternalField[i][j][k];
                //AC_part(2) = 0.0;
                //E_AC += -dot(AC_part, ThisSpin);
                //External[i][j] = -dot(FM_part, ThisSpin);
            
                DM[i][j][k] =         +D*dot(Left, cross(ThisSpin, LeftSpin))/2.0
                                   +D*dot(Right, cross(ThisSpin, RightSpin))/2.0
                                   +D*dot(Up, cross(ThisSpin, UpSpin))/2.0
                                   +D*dot(Down, cross(ThisSpin, DownSpin))/2.0
                                   +D*dot(Top, cross(ThisSpin, TopSpin))/2.0
                                   +D*dot(Bottom, cross(ThisSpin, BottomSpin))/2.0;
                EnergyDistribution[i][j][k] = Heisenberg[i][j][k]  + DM[i][j][k] + External[i][j][k];
                TotalEnergy += EnergyDistribution[i][j][k];
            }
        }
    }
    //TotalMinusAC = TotalEnergy - E_AC;
    return TotalEnergy;
}
/////////////////
/*
void InitializeDotSpinWaveExcitation(vec** Plate, int x, int y, double Sx, double Sy, double Sz)
{
    vec temp(3);
    temp(0) = Sx;
    temp(1) = Sy;
    temp(2) = Sz;
    Plate[x][y] = temp;
}*/
//////////////////////////////////////////////
/*void EnforceExternalFieldOnSite(int x, int y, double Hy, vec** ExternalField)
{
    vec temp(3);
    temp.zeros();
    temp(1) = Hy;
    ExternalField[x][y] = temp;
}*/
////////////////////////////////////////////
/*
double OnSiteEffectiveFieldZ(int i, int j, vec** Plate, NeighborRecord** NeighborInput, vec** ExternalFieldInput, double J, double D,
                           double &HeisenbergZ, double &DMZ, double &ExternalZ, double &TotalXY, double &DMXY)
{
    vec FieldTotal(3);
    vec FieldHeisenberg(3);
    vec FieldDM(3);
    vec FieldExtra(3);
    vec SpinRight(3);
    vec SpinLeft(3);
    vec SpinUp(3);
    vec SpinDown(3);
    vec X(3);
    vec Y(3);
    X.zeros(3);
    Y.zeros(3);
    X(0) = 1.0;
    Y(1) = 1.0;
    SpinRight = Plate[NeighborInput[i][j].Right[0]][NeighborInput[i][j].Right[1]];
    SpinLeft = Plate[NeighborInput[i][j].Left[0]][NeighborInput[i][j].Left[1]];
    SpinUp = Plate[NeighborInput[i][j].Up[0]][NeighborInput[i][j].Up[1]];
    SpinDown = Plate[NeighborInput[i][j].Down[0]][NeighborInput[i][j].Down[1]];
    FieldHeisenberg = J*(SpinRight+SpinLeft+SpinUp+SpinDown);
    FieldDM = - D*cross((SpinRight-SpinLeft), X) - D*cross((SpinUp-SpinDown),Y);
    FieldExtra = ExternalFieldInput[i][j]; 
    FieldTotal = FieldHeisenberg+FieldDM+FieldExtra;
    HeisenbergZ = FieldHeisenberg(2);
    DMZ = FieldDM(2);
    ExternalZ = FieldExtra(2);
    TotalXY = sqrt(FieldTotal(0)*FieldTotal(0)+FieldTotal(1)*FieldTotal(1));
    DMXY = sqrt(FieldDM(0)*FieldDM(0)+FieldDM(1)*FieldDM(1));
    return FieldTotal(2);
}
*/