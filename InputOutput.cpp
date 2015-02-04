#include <stdio.h>
#include <math.h>
#include <complex>
#include <armadillo>
#include <time.h>
#include "discpp.h"
#include "InputOutput.hpp"
#include "NeighborList.hpp"

#define PIPI 3.1415926535897932384626

using namespace std;
using namespace arma;

typedef complex<double> Complex;



vec** Allocate_Vector_Plate(int SizeX, int SizeY)
{
    vec **A = new vec*[SizeX];
    for (int i=0; i<SizeX; i++)
    {
        A[i] = new vec[SizeY];
    }
    return A;
}
/////////////////////////////////////////////////////////////////////////////////////////
void Delete_Plate(vec **p, int SizeX)
{
    for (int i=0; i<SizeX; i++)
    {
        delete [] p[i];
    }
    delete [] p;
}
////////////////////////////////////////////////////////////////////////////////////////
vec*** Allocate_3D_Vector_Plate(int SizeX, int SizeY, int SizeZ)
{
    vec*** A = new vec**[SizeX];
    for (int i=0; i<SizeX; i++)
    {
        A[i] = new vec*[SizeY];
        for (int j=0; j<SizeY; j++)
        {
            A[i][j] = new vec[SizeZ];
        }    
    }
    return A;
}
/////////////////////////////////////////////////////////////////////////////////
void Delete_3D_Plate(vec*** p, int SizeX, int SizeY)
{
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            delete [] p[i][j];
        } 
        delete [] p[i];
    }
    delete [] p;
}
    
/////////////////////////////////////////////////////////////////////////////
double** AllocateDoublePlate(int SizeX, int SizeY)
{
    double **A = new double*[SizeX];
    for (int i=0; i<SizeX; i++)
    {
        A[i] = new double[SizeY];
    }
    return A;
}
////////////////////////////////////////////////////////////////////////////////////////
void DeleteDoublePlate(double** p, int SizeX)
{
    for (int i=0; i<SizeX; i++)
    {
        delete [] p[i];
    }
    delete [] p;
}
////////////////////////////////////////////////////////////////////////////////////
double*** Allocate_3D_Double_Plate(int SizeX, int SizeY, int SizeZ)
{
    double*** A = new double**[SizeX];
    for (int i=0; i<SizeX; i++)
    {
        A[i] = new double*[SizeY];
        for (int j=0; j<SizeY; j++)
        {
            A[i][j] = new double[SizeZ];
        }    
    }
    return A;
}
/////////////////////////////////////////////////////////////////////////////////
void Delete_3D_Double_Plate(double ***p, int SizeX, int SizeY)
{
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            delete [] p[i][j];
        } 
        delete [] p[i];
    }
    delete [] p;
}
////////////////////////////////////////////////////////////////////////////////////
void ReadInTexture(vec** Plate, int SizeX, int SizeY)
{
    FILE *fp;
    fp = fopen("texture.txt","r");
    double X, Y, Z;
    int i, j;
    double x,y;
    // Here, x and y tracks the coordinate sown up there.
    // i, j numbers the sites for the Hamiltonian Block (i,j);
    for (int p=0; p<SizeX*SizeY; p++)
    {
        fscanf(fp, "%lf%lf%lf%lf%lf\n", &x, &y, &X, &Y, &Z);
        vec temp(3);
        temp(0) = X;
        temp(1) = Y;
        temp(2) = Z;
        Plate[(int)x][(int)y] = temp;
    }
    fclose(fp);
 }
//////////////////////////////////////////////////////////////////////////////////////
void InitializeFM_With_Perturbation(vec*** Plate, int SizeX, int SizeY, int SizeZ, double PerturbationFactor)
{
    // PerturbationFactor is a small number that randomizes the FM state.
    srand (time(NULL));
    double ZeroToOne;
    double Theta;
    double Angle;
    double X, Y, Z;
    vec ZeroVector(3);
    ZeroVector.zeros();
    for (int x=0; x<SizeX; x++)
    {
        for (int y=0; y<SizeY; y++)
        {
            Plate[x][y][0] = ZeroVector;
            Plate[x][y][SizeZ-1] = ZeroVector;
            for (int z=1; z<SizeZ-1; z++)
            {
                ZeroToOne = (double)rand()/(double)RAND_MAX;
                Theta = PIPI*ZeroToOne*PerturbationFactor;
                ZeroToOne = (double)rand()/(double)RAND_MAX;
                Angle = ZeroToOne*2.0*PIPI;
                X = cos(Angle)*sin(Theta);
                Y = sin(Angle)*sin(Theta);
                Z = cos(Theta);
                vec temp(3);
                temp(0) = X;
                temp(1) = Y;
                temp(2) = Z;
                Plate[x][y][z] = temp;
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////////
void InitializeSpinPlot(Dislin &Graph, float* &xp, float* &yp, float* &xv, float* &yv, int SizeX, int SizeY)
{

    // allocate the memory
   xp = new float[SizeX];
   yp = new float[SizeY];
   xv = new float[SizeX*SizeY];
   yv = new float[SizeX*SizeY];
  
  int AxisX;
  int AxisY;
  if (SizeX>SizeY)
  {
      AxisX = 2000;
      AxisY = 2000*SizeY/SizeX+1;
  }
  else
  {
      AxisY = 2000;
      AxisX = 2000*SizeX/SizeY+1;
  }
  
  Graph.winsiz ((AxisX/2), AxisY/2);
  Graph.page ((AxisX+600), AxisY+600);
  
  Graph.sclmod ("full");
  Graph.scrmod ("revers");
  Graph.metafl ("xwin");
  Graph.x11mod("nostore");
  
  Graph.disini ();
  Graph.pagera ();
  Graph.hwfont ();
  
  
  Graph.axspos (350, AxisY+300);
  Graph.axslen (AxisX, AxisY);
  
  Graph.name ("X-axis", "x");
  Graph.name ("Y-axis", "y");
  Graph.vecopt(1.0, "scale");
  Graph.vecopt(0.4, "size");
  Graph.vecopt(1.0, "length");
  Graph.selwin(1);
  Graph.graf (-1, SizeX, 0, 10, -1, SizeY, 0, 10);

  Graph.height (50);
  Graph.vecclr (-2);
}  
//////////////////////////////////////////////////////
void UpdateSpinPlot(Dislin &Graph, vec*** Plate, double rate, int count, double TimeStep, double TopologicalCharge, 
               float* xp, float* yp, float* xv, float* yv, int SizeX, int SizeY, int Z, double TotalEnergy)
{
  
    Graph.erase();
    
    
    vec temp(3);
    char buffer[256];
    float xstep, ystep;
    xstep = 2.0 / (SizeX - 1);; 
    ystep = 2.0 / (SizeY - 1);;


    // memory allocation finished.
    // initialize the values of the starting points. 
   
    for (int i=0; i<SizeX; i++)
        xp[i] = (float)i;//(float) (-1. +  i * xstep);
    for (int i=0; i<SizeY; i++)
        yp[i] = (float)i;// (-1. +  i * ystep);
    
    for (int x=0; x<SizeX; x++)
    {
        for (int y=0; y<SizeY; y++)
        {
            temp = Plate[x][y][Z];
            xv[x*SizeY+y]=(float)temp(0);
            yv[x*SizeY+y]=(float)temp(1);
        }
    }
    
    sprintf(buffer, "Hc=%.4lf (%.4lf) TopologicalCharge=% .4lf  %d E=% .4lf", 
                 (double)count*TimeStep*rate, (double)count*TimeStep, TopologicalCharge, count, TotalEnergy);
    Graph.titlin (buffer, 4);
    Graph.title ();
     
    Graph.vecmat (xv, yv, SizeX, SizeY, xp, yp, 1921);
    //endgrf();
    Graph.sendbf();
    //disfin();

}
/////////////////////////////////////////////////
void FinalizeSpinPlot(Dislin &Graph, float* xp, float* yp, float* xv, float* yv)
{
    delete [] xp;
    delete [] yp;
    delete [] xv;
    delete [] yv;
    Graph.disfin();
}
/////////////////////////////////////////////////////
void InitializeExternalField(vec*** ExternalFieldOutput, int SizeX, int SizeY, int SizeZ, vec h0)
{
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                ExternalFieldOutput[i][j][k] = h0;
            }
        }
    }
}
////////////////////////////////////////////////////////////////

void InitializeChargePlot(Dislin &Graph, float* &pointer, int SizeX, int SizeY, int PlotX, int PlotY)
{
    pointer = new float[PlotX*PlotY];
    int AxisX;
    int AxisY;
    if (PlotX>PlotY)
    {
      AxisX = 2000;
      AxisY = 2000*SizeY/SizeX+1;
    }
    else
    {
      AxisY = 2000;
      AxisX = 2000*SizeX/SizeY+1;
    }
  
    Graph.winsiz ((AxisX/2+200), AxisY/2+200);
    Graph.page ((AxisX+800), AxisY+600);
  
    Graph.sclmod ("full");
    Graph.scrmod ("revers");
    Graph.metafl ("xwin");
    Graph.x11mod("nostore");
  
    Graph.disini ();
    Graph.pagera ();
    Graph.hwfont ();
  
  
    Graph.axspos (350, AxisY+300);
    Graph.axslen (AxisX, AxisY);
  
    Graph.name ("X-axis", "x");
    Graph.name ("Y-axis", "y");
    Graph.vecopt(1.0, "scale");
    Graph.vecopt(0.4, "size");
    Graph.vecopt(1.0, "length");
    Graph.autres(PlotX, PlotY);
    Graph.height (50);
    Graph.title  ();
}
////////////////////////////////////////////////////////////////////
void UpdateChargePlot(Dislin &Graph, float* PlotPointer, double*** TopologicalChargeDistribution, double omega, int count, 
                     double TimeStep, double TopologicalCharge, int SizeX, int SizeY, int Z, int PlotX, int PlotY)
{
    Graph.erase();
    Graph.graf3  (SizeX/2-PlotX/2, SizeX/2+PlotX/2, SizeX/2-PlotX/2, 2, SizeY/2-PlotY/2, SizeY/2+PlotY/2, SizeY/2-PlotY/2, 2,
            -0.5, 0.5, -0.5, 0.1);
    char buffer[256];
    for (int i=0; i<PlotX; i++)
    {
        for(int j=0; j<PlotY; j++)
        {
            int x = SizeX/2-PlotX/2+i;
            int y = SizeY/2-PlotY/2+j;
            PlotPointer[i*PlotY+j] = (float)TopologicalChargeDistribution[x][y][Z];
        }
    }
    sprintf(buffer, "t=%.4lfT (%.4lf) TopologicalCharge=%.4lf Step=%d", (double)count*TimeStep/2.0/PIPI*omega, (double)count*TimeStep, TopologicalCharge, count);
    Graph.titlin (buffer, 4);
    Graph.title ();
    Graph.crvmat (PlotPointer, PlotX, PlotY, 15, 15);
    Graph.endgrf();
}
/////////////////////////////////////////////////////////////////
void FinalizeChargePlot(Dislin &Graph, float* PlotPointer)
{
    delete [] PlotPointer;
    Graph.disfin();
}
///////////////////////////////////////////////////////////////////
void WriteToTextFiles(vec*** Plate, int SizeX, int SizeY, int Z, char* filename)
{
    FILE* fp;
    fp = fopen(filename, "w");
    vec temp(3);
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            temp = Plate[i][j][Z];
            fprintf(fp, "%d\t%d\t% le\t% le\t% le\n", i, j, temp(0), temp(1), temp(2));
        }
    }
    fclose(fp);
}
//////////////////////////////////////////////////////////////////
void WriteSzToTextFilex(vec***Plate, int SizeX, int SizeY, int Z, char* filename)
{
    FILE* fp;
    fp = fopen(filename,"w");
    vec temp(3);
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            temp = Plate[i][j][Z];
            fprintf(fp, "%d\t%d\t%lf\n", i, j, temp(2));
        }
    }
    fclose(fp);
}
//////////////////////////////////////////////////////////////////////////////////
void ReadInEnforcementList(char* InputFileName, EnforcedSpin* &OutputList, int &OutputLength)
{
    FILE* fp;
    fp = fopen(InputFileName,"r");
    int count=0;
    int x, y, z;
    double Sx, Sy, Sz;
    vec temp(3);
    while (!feof(fp))
    {
        fscanf(fp, "%d%d%d%le%le%le", &x, &y, &z, &Sx, &Sy, &Sz);
        count++;
        printf("count = %d\n", count);
    }
    rewind(fp);
    OutputLength = count;
    OutputList = new EnforcedSpin[OutputLength];
    for (int i=0; i<OutputLength; i++)
    {
        fscanf(fp, "%d%d%d%le%le%le", &x, &y, &z, &Sx, &Sy, &Sz);
        OutputList[i].x = x;
        OutputList[i].y = y;
        OutputList[i].z = z;
        temp(0) = Sx;
        temp(1) = Sy; 
        temp(2) = Sz;
        OutputList[i].spin = temp;
    }
    fclose(fp);
    
}
//////////////////////////////////////////////////////////////////////////////////
/*void OutputOnSiteData(vec** Plate, int CurrentStep, int x1, int y1, int x2, int y2,
                      int x3, int y3, int x4, int y4, int x5, int y5, int StartingStep, int EndingStep, FILE* fp)
{
    // This function tracks the spin/effectiveField on site (x,y). 
    // The result will be stored in the FILE pointer fp;
    if (CurrentStep >= StartingStep && CurrentStep <= EndingStep)
    {
        fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
                                 (Plate[x1][y1])(0), (Plate[x1][y1])(1), (Plate[x1][y1])(2),
                                 (Plate[x2][y2])(0), (Plate[x2][y2])(1), (Plate[x2][y2])(2),
                                 (Plate[x3][y3])(0), (Plate[x3][y3])(1), (Plate[x3][y3])(2),
                                 (Plate[x4][y4])(0), (Plate[x4][y4])(1), (Plate[x4][y4])(2),
                                 (Plate[x5][y5])(0), (Plate[x5][y5])(1), (Plate[x5][y5])(2));
    }
}*/
///////////////////
void WriteTextureToFile(vec*** Plate, int SizeX, int SizeY, int Z, int PlotX, int PlotY, double Q)
{
    FILE* fp;
    fp = fopen("SpinTexture.txt", "w");
    vec temp(3);
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            temp = Plate[i][j][Z];
            fprintf(fp, "%d\t%d\t% le\t% le\t% le\n", i, j, temp(0), temp(1), temp(2) );
        }
    }
    fclose(fp);
   float* xp = new float[PlotX];
   float* yp = new float[PlotY];
   float* xv = new float[PlotX*PlotY];
   float* yv = new float[PlotX*PlotY];
   float* pointer = new float[PlotX*PlotY];
  int AxisX;
  int AxisY;
  if (PlotX>PlotY)
  {
      AxisX = 2000;
      AxisY = 2000*SizeY/SizeX+1;
  }
  else
  {
      AxisY = 2000;
      AxisX = 2000*PlotX/PlotY+1;
  }
  
  Dislin Graph;
  Graph.winsiz ((AxisX/2), AxisY/2);
  Graph.page ((AxisX+600), AxisY+600);
  
  Graph.sclmod ("full");
  Graph.scrmod ("black");
  Graph.metafl ("png");
  Graph.x11mod("nostore");
  
  Graph.disini ();
  Graph.pagera ();
  Graph.hwfont ();
  
  
  Graph.axspos (350, AxisY+300);
  Graph.axslen (AxisX, AxisY);
  
  Graph.name ("X-axis", "x");
  Graph.name ("Y-axis", "y");
  Graph.vecopt(1.0, "scale");
  Graph.vecopt(1.0-0.618, "size");
  Graph.vecopt(1.0, "length");
  Graph.vecopt(18.0, "angle");
  Graph.selwin(1);
  Graph.graf (SizeX/2-PlotX/2, SizeX/2+PlotX/2, 0, 10, SizeY/2-PlotY/2, SizeY/2+PlotY/2, 0, 10);

  Graph.height (50);
  Graph.vecclr (-2);
 // vec temp(3);
    char buffer[256];
    float xstep, ystep;
    xstep = 2.0 / (PlotX - 1);; 
    ystep = 2.0 / (PlotY - 1);;


    // memory allocation finished.
    // initialize the values of the starting points. 
   
    for (int i=0; i<PlotX; i++)
        xp[i] = (float)(SizeX/2-PlotX/2+i);//(float) (-1. +  i * xstep);
    for (int i=0; i<PlotY; i++)
        yp[i] = (float)(SizeY/2-PlotY/2+i);// (-1. +  i * ystep);
    
    for (int x=0; x<PlotX; x++)
    {
        for (int y=0; y<PlotY; y++)
        {
            temp = Plate[SizeX/2-PlotX/2+x][SizeY/2-PlotY/2+y][Z];
            xv[x*PlotY+y]=(float)temp(0);
            yv[x*PlotY+y]=(float)temp(1);
            pointer[x*PlotY+y]=(float)temp(2);
        }
    }
    
    sprintf(buffer, "Q = %.1lf", Q);
    Graph.titlin (buffer, 4);
    Graph.title ();
    Graph.vecmat (xv, yv, PlotX, PlotY, xp, yp, 1901);
    //Graph.crvmat (pointer, PlotX, PlotY, 15, 15);
    
    
    //endgrf();
    Graph.sendbf();
  delete [] xp;
  delete [] yp;
  delete [] xv;
  delete [] yv;
  delete [] pointer;
  Graph.disfin();
}
//////////////////////////////////////////
void RecordVibration(int x, int y, int z, int count, vec*** Plate, cx_vec &Record)
{
    if (count < (int)Record.n_elem)
        Record(count) = (Plate[x][y][z])(1);
}
/////////////////////////////////////////
vec TotalMagnitization(vec*** Plate, int SizeX, int SizeY, int SizeZ, double Radius, double CenterX, double CenterY)
{
    vec Total(3);
    Total.zeros();
    int count = 0;
    for (int i=0/*(int)(floor(CenterX))*//*(int)(CenterX-Radius)*/; i<SizeX/*=(int)(CenterX+Radius)*//*(int)(ceil(CenterX))*/; i++)
    {
        for (int j=0/*(int)(floor(CenterY))*//*(int)(CenterY-Radius)*/; j<SizeY/*=(int)(CenterY+Radius)*//*(int)(ceil(CenterY))*/; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
           // double DistanceSquare = ((double)i-CenterX)*((double)i-CenterX)+((double)j-CenterY)*((double)j-CenterY);
           // if (DistanceSquare < Radius*Radius)
           // {
                Total = Total + Plate[i][j][k];
                count++;
            }
            //}
        }
    }
   // printf("%d\t%d\t%d\n", count, (int)(floor(CenterX)), (int)(ceil(CenterX)));
    return Total/(double)count;
}
//////////////////////////////
void WriteTheEnergyAlongX(double*** EnergyDistribution, double*** Heisenberg, double*** DM, double*** External, 
        int SizeX, int Y, int Z, char* filename, vec h0, double J)
{
    FILE* fp;
    fp = fopen(filename, "w");
    for (int i=0; i<SizeX; i++)
    {
         fprintf(fp, "%d\t% le\t% le\t% le\t% le\t% le\n", i, EnergyDistribution[i][Y][Z]-(-h0(2)-J*2.0), Heisenberg[i][Y][Z], 
                                             DM[i][Y][Z], External[i][Y][Z], EnergyDistribution[i][Y][Z]);
    }
    fclose(fp);
}
/////////////////////////////////////
/*
void WriteEffectiveFieldAlongX(vec*** Plate, NeighborRecord*** NeighborInput, vec*** ExternalFieldInput, 
                             int SizeX, int Y_Coordinate, int Z_Coordinate, char* filename, 
                             vec h0, double J, double D)
{
    FILE* fp;
    fp = fopen(filename, "w");
    vec FieldTotal(3);
    vec FieldHeisenberg(3);
    vec FieldDM(3);
    vec FieldExtra(3);
    vec SpinRight(3);
    vec SpinLeft(3);
    vec SpinUp(3);
    vec SpinDown(3);
    vec SpinTop(3);
    vec SpinBottom(3);
    vec X(3);
    vec Y(3);
    vec Z(3);
    X.zeros(3);
    Y.zeros(3);
    Z.zeros(3);
    X(0) = 1.0;
    Y(1) = 1.0;
    Z(2) = 1.0;
    for (int i=0; i<SizeX; i++)
    {
        SpinRight = Plate[NeighborInput[i][Y_Coordinate][Z_Coordinate].Right[0]][NeighborInput[i][Y_Coordinate][Z_Coordinate].Right[1]][Z_Coordinate];
        SpinLeft = Plate[NeighborInput[i][Y_Coordinate][Z_Coordinate].Left[0]][NeighborInput[i][Y_Coordinate][Z_Coordinate].Left[1]][Z_Coordinate];
        SpinUp = Plate[NeighborInput[i][Y_Coordinate][Z_Coordinate].Up[0]][NeighborInput[i][Y_Coordinate][Z_Coordinate].Up[1]][Z_Coordinate];
        SpinDown = Plate[NeighborInput[i][Y_Coordinate][Z_Coordinate].Down[0]][NeighborInput[i][Y_Coordinate][Z_Coordinate].Down[1]][Z_Coordinate];
        SpinTop = Plate[NeighborInput[i][Y_Coordinate][Z_Coordinate].Top[0]][NeighborInput[i][Y_Coordinate][Z_Coordinate].Top[1]][Z_Coordinate];
        SpinBottom = Plate[NeighborInput[i][Y_Coordinate][Z_Coordinate].Bottom[0]][NeighborInput[i][Y_Coordinate][Z_Coordinate].Bottom[1]][Z_Coordinate];
        FieldHeisenberg = J*(SpinRight+SpinLeft+SpinUp+SpinDown);
        FieldDM = - D*cross((SpinRight-SpinLeft), X) - D*cross((SpinUp-SpinDown),Y);
        FieldExtra = ExternalFieldInput[i][Y_Coordinate]; 
        FieldTotal = FieldHeisenberg+FieldDM+FieldExtra;
        fprintf(fp, "% lf\t% lf\t% lf\t% lf\t% lf", (double)i, FieldTotal(2), FieldHeisenberg(2), FieldDM(2), FieldExtra(2));
        fprintf(fp, "% lf\t% lf\t% lf\n", sqrt(FieldDM(0)*FieldDM(0)+FieldDM(1)*FieldDM(1)), 
                                     sqrt(FieldHeisenberg(0)*FieldHeisenberg(0)+FieldHeisenberg(1)*FieldHeisenberg(1)),
                                     sqrt(FieldExtra(0)*FieldExtra(0)+FieldExtra(1)*FieldExtra(1)));
    }
    fclose(fp);
}
//////////////////////////////
*/

void InitializeSpinPlotYZ(Dislin &Graph, float* &xp, float* &yp, float* &xv, float* &yv, int SizeY, int SizeZ)
{

    // allocate the memory
   //SizeZ = SizeZ-2;
   xp = new float[SizeY];
   yp = new float[SizeZ];
   xv = new float[SizeY*(SizeZ)];
   yv = new float[SizeY*(SizeZ)];
  
  int AxisY;
  int AxisZ;
  if (SizeY>SizeZ)
  {
      AxisY = 2000;
      AxisZ = 2000*SizeZ/SizeY+1;
  }
  else
  {
      AxisZ = 2000;
      AxisY = 2000*SizeY/SizeZ+1;
  }

  
  Graph.winsiz ((AxisY/2), AxisZ);
  Graph.page ((AxisY+600), AxisZ+600);
  
  Graph.sclmod ("full");
  Graph.scrmod ("revers");
  Graph.metafl ("xwin");
  Graph.x11mod("nostore");
  
  Graph.disini ();
  Graph.pagera ();
  Graph.hwfont ();
  
  
  Graph.axspos (350, AxisZ+300);
  Graph.axslen (AxisY, AxisZ);
  
  Graph.name ("Y-axis", "y");
  Graph.name ("Z-axis", "z");
  Graph.vecopt(1.0, "scale");
  Graph.vecopt(0.4, "size");
  Graph.vecopt(1.0, "length");
  Graph.selwin(1);
  Graph.graf (-1, SizeY, 0, 10, -1, SizeZ, 0, 10);

  Graph.height (50);
  Graph.vecclr (-2);
}
///////////////
void UpdateSpinPlotYZ(Dislin &Graph, vec*** Plate, double rate, int count, double TimeStep, double TopologicalCharge, 
               float* xp, float* yp, float* xv, float* yv, int SizeY, int SizeZ, int X, double TotalEnergy)
{
  
    Graph.erase();
    //SizeZ = SizeZ;
    
    vec temp(3);
    char buffer[256];
    float ystep, zstep;
    ystep = 2.0 / (SizeY - 1);
    zstep = 2.0 / (SizeZ - 1);
   /* if (ystep < zstep)
    {
        zstep = ystep;
    }
    else
    {
        ystep = zstep;
    }*/

    // memory allocation finished.
    // initialize the values of the starting points. 
   
    for (int i=0; i<SizeY; i++)
        xp[i] = (float)i;//(float) (-1. +  i * xstep);
    for (int i=0; i<SizeZ; i++)
        yp[i] = (float)i;// (-1. +  i * ystep);
    
    for (int y=0; y<SizeY; y++)
    {
        for (int z=0; z<SizeZ; z++)
        {
            temp = Plate[X][y][z];
            xv[y*SizeZ+z]=(float)temp(1);
            yv[y*SizeZ+z]=(float)temp(2);
        }
    }
    
    sprintf(buffer, "Hc=%.4lf (%.4lf) TopologicalCharge=% .4lf  %d E=% .4lf", 
                 (double)count*TimeStep*rate, (double)count*TimeStep, TopologicalCharge, count, TotalEnergy);
    Graph.titlin (buffer, 4);
    Graph.title ();
     
    Graph.vecmat (xv, yv, SizeY, SizeZ, xp, yp, 1921);
    //endgrf();
    Graph.sendbf();
    //disfin();

}
/////////////////////////////////////////////////
void WriteTextureToFileYZ(vec*** Plate, int SizeY, int SizeZ, int X, int PlotY, int PlotZ, double Q)
{
   vec temp(3);
   float* xp = new float[PlotY];
   float* yp = new float[PlotZ];
   float* xv = new float[PlotY*PlotZ];
   float* yv = new float[PlotY*PlotZ];
   float* pointer = new float[PlotY*PlotZ];
  int AxisY;
  int AxisZ;
  if (PlotY>PlotZ)
  {
      AxisY = 2000;
      AxisZ = 2000*PlotZ/PlotY+1;
  }
  else
  {
      AxisZ = 2000;
      AxisY = 2000*PlotY/PlotZ+1;
  }
  
  Dislin Graph;
  Graph.winsiz ((AxisY/2), AxisZ);
  Graph.page ((AxisY+600), AxisZ+600);
  
  Graph.sclmod ("full");
  Graph.scrmod ("black");
  Graph.metafl ("eps");
  Graph.x11mod("nostore");
  
  Graph.disini ();
  Graph.pagera ();
  Graph.hwfont ();
  
  
  Graph.axspos (350, AxisY+300);
  Graph.axslen (AxisY, AxisZ);
  
  Graph.name ("Z-axis", "y");
  Graph.name ("Y-axis", "x");
  Graph.vecopt(1.0, "scale");
  Graph.vecopt(1.0-0.618, "size");
  Graph.vecopt(1.0, "length");
  Graph.vecopt(18.0, "angle");
  Graph.selwin(1);
  Graph.graf (SizeY/2-PlotY/2, SizeY/2+PlotY/2, 0, 10, SizeZ/2-PlotZ/2, SizeZ/2+PlotZ/2, 0, 10);

  Graph.height (50);
  Graph.vecclr (-2);
 // vec temp(3);
    char buffer[256];
    float ystep, zstep;
    ystep = 2.0 / (PlotY - 1);; 
    zstep = 2.0 / (PlotZ - 1);;


    // memory allocation finished.
    // initialize the values of the starting points. 
   
    for (int i=0; i<PlotY; i++)
        xp[i] = (float)(SizeY/2-PlotY/2+i);//(float) (-1. +  i * xstep);
    for (int i=0; i<PlotZ; i++)
        yp[i] = (float)(SizeZ/2-PlotZ/2+i);// (-1. +  i * ystep);
    
    for (int y=0; y<PlotY; y++)
    {
        for (int z=0; z<PlotZ; z++)
        {
            temp = Plate[X][SizeY/2-PlotY/2+y][SizeZ/2-PlotZ/2+z];
            xv[y*PlotZ+z]=(float)temp(1);
            yv[y*PlotZ+z]=(float)temp(2);
            pointer[y*PlotZ+z]=(float)temp(2);
        }
    }
    
    sprintf(buffer, "Q = %.1lf", Q);
    Graph.titlin (buffer, 4);
    Graph.title ();
    Graph.vecmat (xv, yv, PlotY, PlotZ, xp, yp, 1901);
    //Graph.crvmat (pointer, PlotX, PlotY, 15, 15);
    
    
    //endgrf();
    Graph.sendbf();
  delete [] xp;
  delete [] yp;
  delete [] xv;
  delete [] yv;
  delete [] pointer;
  Graph.disfin();
}
