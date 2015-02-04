#include <armadillo>
#include "InputOutput.hpp"
#include "NeighborList.hpp"
#include <unistd.h>

using namespace arma;

/*struct NeighborRecord
{
    int Right[3];
    int Left[3];
    int Up[3];
    int Down[3];
    int Top[3];
    int Bottom[3];
};*/
/////////////////////////////////////////

NeighborRecord*** GenerateNeighbors(int SizeX, int SizeY, int SizeZ)
{
    NeighborRecord*** pointer;
    pointer = new NeighborRecord**[SizeX];
    for (int i=0; i<SizeX; i++)
    {
        pointer[i] = new NeighborRecord*[SizeY];
        for (int j=0; j<SizeY; j++)
        {
            pointer[i][j] = new NeighborRecord[SizeZ];
        }
    }
    for (int i=1; i<SizeX-1; i++)
    {
        for (int j=1; j<SizeY-1; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                pointer[i][j][k].Right[0] = i+1;
                pointer[i][j][k].Right[1] = j;
                pointer[i][j][k].Right[2] = k;
                pointer[i][j][k].Left[0] = i-1;
                pointer[i][j][k].Left[1] = j;
                pointer[i][j][k].Left[2] = k;
                pointer[i][j][k].Up[0] = i;
                pointer[i][j][k].Up[1] = j+1;
                pointer[i][j][k].Up[2] = k;
                pointer[i][j][k].Down[0] = i;
                pointer[i][j][k].Down[1] = j-1;
                pointer[i][j][k].Down[2] = k;
                pointer[i][j][k].Top[0] = i;
                pointer[i][j][k].Top[1] = j;
                pointer[i][j][k].Top[2] = k+1;
                pointer[i][j][k].Bottom[0] = i;
                pointer[i][j][k].Bottom[1] = j;
                pointer[i][j][k].Bottom[2] = k-1;
            }
        }
    }
    // Now deal with the boundaries:
    // First, Left boundary:
    for (int j=1; j<SizeY-1; j++)
    {
        for (int k=1; k<SizeZ-1; k++)
        {
            pointer[0][j][k].Right[0] = 1;
            pointer[0][j][k].Right[1] = j;
            pointer[0][j][k].Right[2] = k;
            pointer[0][j][k].Left[0] = SizeX-1;
            pointer[0][j][k].Left[1] = j;
            pointer[0][j][k].Left[2] = k;
            pointer[0][j][k].Up[0] = 0;
            pointer[0][j][k].Up[1] = j+1;
            pointer[0][j][k].Up[2] = k;
            pointer[0][j][k].Down[0] = 0;
            pointer[0][j][k].Down[1] = j-1;
            pointer[0][j][k].Down[2] = k;
            pointer[0][j][k].Top[0] = 0;
            pointer[0][j][k].Top[1] = j;
            pointer[0][j][k].Top[2] = k+1;
            pointer[0][j][k].Bottom[0] = 0;
            pointer[0][j][k].Bottom[1] = j;
            pointer[0][j][k].Bottom[2] = k-1;
        }
    }
    // Now, right boundary:
    for (int j=1; j<SizeY-1; j++)
    {
        for (int k=1; k<SizeZ-1; k++)
        {
            pointer[SizeX-1][j][k].Right[0] = 0;
            pointer[SizeX-1][j][k].Right[1] = j;
            pointer[SizeX-1][j][k].Right[2] = k;
            pointer[SizeX-1][j][k].Left[0] = SizeX-2;
            pointer[SizeX-1][j][k].Left[1] = j;
            pointer[SizeX-1][j][k].Left[2] = k;
            pointer[SizeX-1][j][k].Up[0] = SizeX-1;
            pointer[SizeX-1][j][k].Up[1] = j+1;
            pointer[SizeX-1][j][k].Up[2] = k;
            pointer[SizeX-1][j][k].Down[0] = SizeX-1;
            pointer[SizeX-1][j][k].Down[1] = j-1;
            pointer[SizeX-1][j][k].Down[2] = k;
            pointer[SizeX-1][j][k].Top[0] = SizeX-1;
            pointer[SizeX-1][j][k].Top[1] = j;
            pointer[SizeX-1][j][k].Top[2] = k+1;
            pointer[SizeX-1][j][k].Bottom[0] = SizeX-1;
            pointer[SizeX-1][j][k].Bottom[1] = j;
            pointer[SizeX-1][j][k].Bottom[2] = k-1;
        }
    }
    // Now, bottom boundary:
    for (int i=1; i<SizeX-1; i++)
    {
        for (int k=1; k<SizeZ-1; k++)
        {
            pointer[i][0][k].Right[0] = i+1;
            pointer[i][0][k].Right[1] = 0;
            pointer[i][0][k].Right[2] = k;
            pointer[i][0][k].Left[0] = i-1;
            pointer[i][0][k].Left[1] = 0;
            pointer[i][0][k].Left[2] = k;
            pointer[i][0][k].Up[0] = i;
            pointer[i][0][k].Up[1] = 1;
            pointer[i][0][k].Up[2] = k;
            pointer[i][0][k].Down[0] = i;
            pointer[i][0][k].Down[1] = SizeY-1;
            pointer[i][0][k].Down[2] = k;
            pointer[i][0][k].Top[0] = i;
            pointer[i][0][k].Top[1] = 0;
            pointer[i][0][k].Top[2] = k+1;
            pointer[i][0][k].Bottom[0] = i;
            pointer[i][0][k].Bottom[1] = 0;
            pointer[i][0][k].Bottom[2] = k-1;
        }
    }
    // Now, top boundary:
    for (int i=1; i<SizeX-1; i++)
    {
        for (int k=1; k<SizeZ-1; k++)
        {
            pointer[i][SizeY-1][k].Right[0] = i+1;
            pointer[i][SizeY-1][k].Right[1] = SizeY-1;
            pointer[i][SizeY-1][k].Right[2] = k;
            pointer[i][SizeY-1][k].Left[0] = i-1;
            pointer[i][SizeY-1][k].Left[1] = SizeY-1;
            pointer[i][SizeY-1][k].Left[2] = k;
            pointer[i][SizeY-1][k].Up[0] = i;
            pointer[i][SizeY-1][k].Up[1] = 0;
            pointer[i][SizeY-1][k].Up[2] = k;
            pointer[i][SizeY-1][k].Down[0] = i;
            pointer[i][SizeY-1][k].Down[1] = SizeY-2;
            pointer[i][SizeY-1][k].Down[2] = k;
            pointer[i][SizeY-1][k].Top[0] = i;
            pointer[i][SizeY-1][k].Top[1] = SizeY-1;
            pointer[i][SizeY-1][k].Top[2] = k+1;
            pointer[i][SizeY-1][k].Bottom[0] = i;
            pointer[i][SizeY-1][k].Bottom[1] = SizeY-1;
            pointer[i][SizeY-1][k].Bottom[2] = k-1;
        }
    }
    // Now, the left-bottom corner:
    for (int k=1; k<SizeZ; k++)
    {
        pointer[0][0][k].Right[0] = 1;
        pointer[0][0][k].Right[1] = 0;
        pointer[0][0][k].Right[2] = k;
        pointer[0][0][k].Left[0] = SizeX-1;
        pointer[0][0][k].Left[1] = 0;
        pointer[0][0][k].Left[2] = k;
        pointer[0][0][k].Up[0] = 0;
        pointer[0][0][k].Up[1] = 1;
        pointer[0][0][k].Up[2] = k;
        pointer[0][0][k].Down[0] = 0;
        pointer[0][0][k].Down[1] = SizeY-1;
        pointer[0][0][k].Down[2] = k;
        pointer[0][0][k].Top[0] = 0;
        pointer[0][0][k].Top[1] = 0;
        pointer[0][0][k].Top[2] = k+1;
        pointer[0][0][k].Bottom[0] = 0;
        pointer[0][0][k].Bottom[1] = 0;
        pointer[0][0][k].Bottom[2] = k-1;
        
    // Now, the right-bottom corner:
        pointer[SizeX-1][0][k].Right[0] = 0;
        pointer[SizeX-1][0][k].Right[1] = 0;
        pointer[SizeX-1][0][k].Right[2] = k;
        pointer[SizeX-1][0][k].Left[0] = SizeX-2;
        pointer[SizeX-1][0][k].Left[1] = 0;
        pointer[SizeX-1][0][k].Left[2] = k;
        pointer[SizeX-1][0][k].Up[0] = SizeX-1;
        pointer[SizeX-1][0][k].Up[1] = 1;
        pointer[SizeX-1][0][k].Up[2] = k;
        pointer[SizeX-1][0][k].Down[0] = SizeX-1;
        pointer[SizeX-1][0][k].Down[1] = SizeY-1;
        pointer[SizeX-1][0][k].Down[2] = k;
        pointer[SizeX-1][0][k].Top[0] = SizeX-1;
        pointer[SizeX-1][0][k].Top[1] = 0;
        pointer[SizeX-1][0][k].Top[2] = k+1;
        pointer[SizeX-1][0][k].Bottom[0] = SizeX-1;
        pointer[SizeX-1][0][k].Bottom[1] = 0;
        pointer[SizeX-1][0][k].Bottom[2] = k-1;
    // Now, the left-top corner:
        pointer[0][SizeY-1][k].Right[0] = 1;
        pointer[0][SizeY-1][k].Right[1] = SizeY-1;
        pointer[0][SizeY-1][k].Right[2] = k;
        pointer[0][SizeY-1][k].Left[0] = SizeX-1;
        pointer[0][SizeY-1][k].Left[1] = SizeY-1;
        pointer[0][SizeY-1][k].Left[2] = k;
        pointer[0][SizeY-1][k].Up[0] = 0;
        pointer[0][SizeY-1][k].Up[1] = 0;
        pointer[0][SizeY-1][k].Up[2] = k;
        pointer[0][SizeY-1][k].Down[0] = 0;
        pointer[0][SizeY-1][k].Down[1] = SizeY-2;
        pointer[0][SizeY-1][k].Down[2] = k;
        pointer[0][SizeY-1][k].Top[0] = 0;
        pointer[0][SizeY-1][k].Top[1] = SizeY-1;
        pointer[0][SizeY-1][k].Top[2] = k+1;
        pointer[0][SizeY-1][k].Bottom[0] = 0;
        pointer[0][SizeY-1][k].Bottom[1] = SizeY-1;
        pointer[0][SizeY-1][k].Bottom[2] = k-1;
    // Now, the right-top corner:
        pointer[SizeX-1][SizeY-1][k].Right[0] = 0;
        pointer[SizeX-1][SizeY-1][k].Right[1] = SizeY-1;
        pointer[SizeX-1][SizeY-1][k].Right[2] = k;
        pointer[SizeX-1][SizeY-1][k].Left[0] = SizeX-2;
        pointer[SizeX-1][SizeY-1][k].Left[1] = SizeY-1;
        pointer[SizeX-1][SizeY-1][k].Left[2] = k;
        pointer[SizeX-1][SizeY-1][k].Up[0] = SizeX-1;
        pointer[SizeX-1][SizeY-1][k].Up[1] = 0;
        pointer[SizeX-1][SizeY-1][k].Up[2] = k;
        pointer[SizeX-1][SizeY-1][k].Down[0] = SizeX-1;
        pointer[SizeX-1][SizeY-1][k].Down[1] = SizeY-2;
        pointer[SizeX-1][SizeY-1][k].Down[2] = k;
        pointer[SizeX-1][SizeY-1][k].Top[0] = SizeX-1;
        pointer[SizeX-1][SizeY-1][k].Top[1] = SizeY-1;
        pointer[SizeX-1][SizeY-1][k].Top[2] = k+1;
        pointer[SizeX-1][SizeY-1][k].Bottom[0] = SizeX-1;
        pointer[SizeX-1][SizeY-1][k].Bottom[1] = SizeY-1;
        pointer[SizeX-1][SizeY-1][k].Bottom[2] = k-1;
    }
    return pointer;
}

void DeleteNeighbors(NeighborRecord*** p, int SizeX, int SizeY)
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

void NeighborTest(void)
{
    int SizeX = 8;
    int SizeY = 8;
    int SizeZ = 5;
    float *xp, *yp, *xv, *yv;
    Dislin TestGraph;
    vec*** Plate = Allocate_3D_Vector_Plate(SizeX, SizeY, SizeZ);
    InitializeSpinPlot(TestGraph, xp, yp, xv, yv, SizeX, SizeY);
    NeighborRecord*** Neighbors;
    Neighbors = GenerateNeighbors(SizeX, SizeY, SizeZ);
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                Plate[i][j][k].zeros(3);
            }
        }
    }
    vec X(3), Y(3);
    X.zeros(3);
    Y.zeros(3);
    X(0) = 1.0;
    Y(1) = 1.0;
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                for (int l=0; l<SizeX; l++)
                {
                    for (int m=0; m<SizeY; m++)
                    {
                        for (int n=1; n<SizeZ-1; n++)
                        {
                            Plate[l][m][n].zeros(3);
                        }
                    }
                }
                Plate[Neighbors[i][j][k].Right[0]][Neighbors[i][j][k].Right[1]][Neighbors[i][j][k].Right[2]] = X;
                Plate[Neighbors[i][j][k].Left[0]][Neighbors[i][j][k].Left[1]][Neighbors[i][j][k].Left[2]] = -X;
                Plate[Neighbors[i][j][k].Up[0]][Neighbors[i][j][k].Up[1]][Neighbors[i][j][k].Up[2]] = Y;
                Plate[Neighbors[i][j][k].Down[0]][Neighbors[i][j][k].Down[1]][Neighbors[i][j][k].Down[2]] = -Y;
                Plate[Neighbors[i][j][k].Top[0]][Neighbors[i][j][k].Top[1]][Neighbors[i][j][k].Top[2]] = Y;
                Plate[Neighbors[i][j][k].Bottom[0]][Neighbors[i][j][k].Bottom[1]][Neighbors[i][j][k].Bottom[2]] = -Y;
                UpdateSpinPlot(TestGraph, Plate, 1.0, k, 0.1, 1.0, xp, yp, xv, yv, SizeX, SizeY, 2, 1000);
                sleep(1);
                
            }
        }
    }
    Delete_3D_Plate(Plate, SizeX, SizeY);
    DeleteNeighbors(Neighbors, SizeX, SizeY);
    FinalizeSpinPlot(TestGraph, xp, yp, xv, yv);
}