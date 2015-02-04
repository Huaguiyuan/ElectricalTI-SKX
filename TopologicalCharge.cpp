#include <armadillo>
#include <complex>
#include "NeighborList.hpp"
#define PIPI 3.1415962535897932384626
using namespace std;
using namespace arma;

typedef complex<double> Complex;

double CalculateTopologicalCharge(vec*** PlateInput, double*** TopologicalChargeDensityOutput, NeighborRecord*** Neighbors, 
                              int SizeX, int SizeY, int SizeZ)
{
    // For each square, we devide the lattice in the following way:
    /*       c
             * 
             * *       
             *   *     
             *     *   
   d * * * * a * * * * b
       *     *
         *   *
           * *
             *
             e 
     The three numbers should be a->b->c and a->d->e
     */
    double TopologicalCharge = 0.0;
    int ax, ay, bx, by, cx, cy, dx, dy, ex, ey;
    Complex rho1, SigmaA1, EXP1;
    Complex rho2, SigmaA2, EXP2;
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            for (int k=1; k<SizeZ-1; k++)
            {
                ax = i;
                ay = j;
                bx = Neighbors[i][j][k].Right[0];
                by = Neighbors[i][j][k].Right[1];
                cx = Neighbors[i][j][k].Up[0];
                cy = Neighbors[i][j][k].Up[1];
                dx = Neighbors[i][j][k].Left[0];
                dy = Neighbors[i][j][k].Left[1];
                ex = Neighbors[i][j][k].Down[0];
                ey = Neighbors[i][j][k].Down[1];
                // First calculate a->b->c;
                rho1 = sqrt(2.0*(1+dot(PlateInput[ax][ay][k], PlateInput[bx][by][k]))*
                         (1+dot(PlateInput[bx][by][k], PlateInput[cx][cy][k]))*
                         (1+dot(PlateInput[cx][cy][k], PlateInput[ax][ay][k])));
                EXP1 = (1.0+dot(PlateInput[ax][ay][k], PlateInput[bx][by][k])+dot(PlateInput[bx][by][k], PlateInput[cx][cy][k])+dot(PlateInput[cx][cy][k], PlateInput[ax][ay][k])
                     +Complex(0,1.)*dot(PlateInput[ax][ay][k], cross(PlateInput[bx][by][k], PlateInput[cx][cy][k])))/rho1;
                SigmaA1 = 2.0*imag(log(EXP1));
                // Second, a->d->e
                rho2 = sqrt(2.0*(1+dot(PlateInput[ax][ay][k], PlateInput[dx][dy][k]))*
                         (1+dot(PlateInput[dx][dy][k], PlateInput[ex][ey][k]))*
                         (1+dot(PlateInput[ex][ey][k], PlateInput[ax][ay][k])));
                EXP2 = (1.0+dot(PlateInput[ax][ay][k], PlateInput[dx][dy][k])+dot(PlateInput[dx][dy][k], PlateInput[ex][ey][k])+dot(PlateInput[ex][ey][k], PlateInput[ax][ay][k])
                     +Complex(0,1.)*dot(PlateInput[ax][ay][k], cross(PlateInput[dx][dy][k], PlateInput[ex][ey][k])))/rho2;
                SigmaA2 = 2.0*imag(log(EXP2));
                TopologicalChargeDensityOutput[i][j][k] = real(1.0/4.0/PIPI*(SigmaA1+SigmaA2));
                TopologicalCharge += TopologicalChargeDensityOutput[i][j][k];
            }
        }
    }
    return TopologicalCharge;
}


double CalculateContinuumTopologicalCharge(vec** PlateInput, double** TopologicalChargeDensityOutput, NeighborRecord** Neighbors, 
                              int SizeX, int SizeY)
{
    // For each square, we devide the lattice in the following way:
    /*       c
             * 
             * *       
             *   *     
             *     *   
   d * * * * a * * * * b
       *     *
         *   *
           * *
             *
             e 
     The three numbers should be a->b->c and a->d->e
     */
    double TopologicalCharge = 0.0;
    int ax, ay, bx, by, cx, cy, dx, dy, ex, ey;
    Complex rho1, SigmaA1, EXP1;
    Complex rho2, SigmaA2, EXP2;
    for (int i=0; i<SizeX; i++)
    {
        for (int j=0; j<SizeY; j++)
        {
            ax = i;
            ay = j;
            bx = Neighbors[i][j].Right[0];
            by = Neighbors[i][j].Right[1];
            cx = Neighbors[i][j].Up[0];
            cy = Neighbors[i][j].Up[1];
            dx = Neighbors[i][j].Left[0];
            dy = Neighbors[i][j].Left[1];
            ex = Neighbors[i][j].Down[0];
            ey = Neighbors[i][j].Down[1];
            // First calculate a->b->c;
            vec Sa, Sb, Sc, Sd, Se;
            Sa = PlateInput[ax][ay];
            Sb = PlateInput[bx][by];
            Sc = PlateInput[cx][cy];
            Sd = PlateInput[dx][dy];
            Se = PlateInput[ex][ey];
                 
            TopologicalChargeDensityOutput[i][j] = 
                    1/16./PIPI*(
                    dot(Sa, cross(Sb,Sc))+
                    dot(Sa, cross(Sd,Se))+
                    dot(Sa, cross(Se,Sb))+
                    dot(Sa, cross(Sc,Sd)));
            TopologicalCharge += TopologicalChargeDensityOutput[i][j];
        }
    }
    return TopologicalCharge;
}
