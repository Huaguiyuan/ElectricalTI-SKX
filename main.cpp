/* 
 * File:   main.cpp
 * Author: RAF_Jason
 *
 * Created on 2013年10月1日, 上午10:01
 */

#include "Calculations.hpp"
#include "MovieWindow.hpp"
#include <sstream>

int main (int argc, char** argv )
{
    double Temperature = 0.0;
    double Ef = -0.05;
    double t = -1.5;
    double Miu1, Miu2, Miu4;
    double J_exchange = 1.0;
    double alpha = 0.1;
    double DM = 1.0;
    Miu1 = 0.1;
    vec BackgroundField(3);
    BackgroundField << 0.0 << 0.0 << 1.0;
    SpinSystem SpinTexture("input.txt",J_exchange, DM, alpha);
    ElectronSystem Electrons("input.txt", "openBoundaries.txt",t);
    SpinTexture.ReadInElectronSiteIndex("ElectronSiteIndices.txt");
    printf("test.\n");
    vec ax, ay, az;
    ax << 31.0 << 0.0 << 0.0;
    ay << 0.0 << 31.0 << 0.0;
    az << 0.0 << 0.0 << 31.0;
    SpinTexture.GenerateNeighborList(1.00001, 1.00001, true, true, false, ax, ay, az);
    /*std::vector<MagneticNode>::iterator i;
    for (i=SpinTexture.NodeList.begin(); i!=SpinTexture.NodeList.end(); i++)
    {
        printf("Index=%d, Number of neighbours=%lu\n", i->Index, i->ListOfExchangeNeighbours.size());
        for (std::vector<int>::iterator j=i->ListOfExchangeNeighbours.begin(); j!=i->ListOfExchangeNeighbours.end(); j++)
        {
            printf("    %d", *j);
        }
        printf("\n");
    }*/
    Electrons.CalculateGR(Ef);
    SpinTexture.SetBackgroundField(BackgroundField);
    SpinTexture.SetTemperature(Temperature);
    MovieWindow OutputWindow(0.0, 31.0, 0.0, 31.0, SpinTexture.NumSite);
    char Info[256];
    //OutputWindow.UpdateWindow(Haha, Info.str());
    int count = 0;
    for (double Time=0.0; Time<10000000.0; Time += 0.01)
    {
        cx_mat T = Electrons.ThermalAverageTransmission(0.0, Ef);
        Complex P, Q;
        P = T(3,0)*T(1,0) + T(3,1)*T(1,0) + T(3,2)*T(1,0) + T(1,3)*T(3,0);
        Q = T(3,0)*T(1,2) + T(3,1)*T(1,2) + T(3,2)*T(1,2) + T(3,2)*T(1,3);
        Miu2 = real(P/(P+Q)- Q/(P+Q))*Miu1;
        P = T(1,0)*T(3,1) + T(1,0)*T(3,0) + T(1,2)*T(3,0) + T(1,3)*T(3,0);
        Q = T(1,0)*T(3,2) + T(1,2)*T(3,1) + T(1,2)*T(3,2) + T(1,3)*T(3,2);
        Miu4 = real(P/(P+Q) - Q/(P+Q))*Miu1;
        printf("Time=%lf\tMiu2=%lf\tMiu4=%lf\n", Time, Miu2, Miu4);
        SpinTexture.Evolve(0.01, Time, true);
        Electrons.UpdateHamiltonian(SpinTexture);
        Electrons.RenewGR(Ef);
        count++;
        if (count%1 == 0)
        {
            sprintf(Info, "t=%lf", Time);
            OutputWindow.UpdateWindow(SpinTexture, Info);
        }
    }
}
    