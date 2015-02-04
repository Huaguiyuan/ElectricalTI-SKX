/* 
 * File:   InputOutput.hpp
 * Author: RAF_Jason
 *
 * Created on 2013年10月1日, 上午10:13
 */

#ifndef INPUTOUTPUT_HPP
#define	INPUTOUTPUT_HPP

#include <armadillo>
#include "discpp.h"
#include "NeighborList.hpp"
using namespace arma;
void ReadInTexture(vec** Plate, int SizeX, int SizeY);
void Delete_Plate(vec **p, int SizeX);
vec** Allocate_Vector_Plate(int SizeX, int SizeY);
vec*** Allocate_3D_Vector_Plate(int SizeX, int SizeY, int SizeZ);
void Delete_3D_Plate(vec ***p, int SizeX, int SizeY);
void InitializeFM_With_Perturbation(vec*** Plate, int SizeX, int SizeY, int SizeZ, double PerturbationFactor);
void InitializeSpinPlot(Dislin &Graph, float* &xp, float* &yp, float* &xv, float* &yv, int SizeX, int SizeY );
void InitializeSpinPlotYZ(Dislin &Graph, float* &xp, float* &yp, float* &xv, float* &yv, int SizeX, int SizeY);
void UpdateSpinPlot(Dislin &Graph, vec*** Plate, double rate, int count, double TimeStep, double TopologicalCharge, 
                   float* xp, float* yp, float* xv, float* yv, int SizeX, int SizeY, int Z, double TotalEnergy);
void UpdateSpinPlotYZ(Dislin &Graph, vec*** Plate, double rate, int count, double TimeStep, double TopologicalCharge, 
                   float* xp, float* yp, float* xv, float* yv, int SizeX, int SizeY, int Z, double TotalEnergy);
void FinalizeSpinPlot(Dislin &Graph, float* xp, float* yp, float* xv, float* yv);
void FinalizeSpinPlotYZ(Dislin &Graph, float* xp, float* yp, float* xv, float* yv);
void InitializeExternalField(vec*** ExternalFieldOutput, int SizeX, int SizeY, int SizeZ, vec h0);
double** AllocateDoublePlate(int SizeX, int SizeY);
void DeleteDoublePlate(double** p, int SizeX);
double*** Allocate_3D_Double_Plate(int SizeX, int SizeY, int NLayers);
void Delete_3D_Double_Plate(double ***p, int SizeX, int NLayers);
void InitializeChargePlot(Dislin &Graph, float* &pointer, int SizeX, int SizeY, int PlotX, int PlotY);
void UpdateChargePlot(Dislin &Graph, float* PlotPointer, double*** TopologicalChargeDistribution, double omega, int count, 
                     double TimeStep, double TopologicalCharge, int SizeX, int SizeY, int Z, int PlotX, int PlotY);
void FinalizeChargePlot(Dislin &Graph, float* PlotPointer);
void WriteToFiles(void);
void OutputOnSiteData(vec** Plate, int CurrentStep, int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4, int x5, int y5, 
                      int StartingStep, int EndingStep, FILE* fp);
void WriteTextureToFile(vec*** Plate, int SizeX, int SizeY, int Z, int PlotX, int PlotY, double Q);
void WriteTextureToFileYZ(vec*** Plate, int SizeY, int SizeZ, int X, int PlotY, int PlotZ, double Q);
struct EnforcedSpin
{
    vec spin;
    int x;
    int y;
    int z;
};
void ReadInEnforcementList(char* InputFileName, EnforcedSpin* &OutputList, int &OutputLength);
void RecordVibration(int x, int y, int count, vec** Plate, cx_vec &Record);
void WriteToTextFiles(vec*** Plate, int SizeX, int SizeY, int Z, char* filename);
vec TotalMagnitization(vec*** Plate, int SizeX, int SizeY, int SizeZ, double Radius, double CenterX, double CenterY);
void WriteTheEnergyAlongX(double*** EnergyDistribution, double*** Heisenberg, double*** DM, double*** External, 
        int SizeX, int Y, int Z, char* filename, vec h0, double J);
void WriteSzToTextFilex(vec**Plate, int SizeX, int SizeY, char* filename);
void WriteEffectiveFieldAlongX(vec** Plate, NeighborRecord** NeighborInput, vec** ExternalFieldInput, 
                             int SizeX, int Y_Coordinate, char* filename, 
                             vec h0, double J, double D);
//void WriteTextureToFile(vec** Plate, int SizeX, int SizeY);
#endif	/* INPUTOUTPUT_HPP */

