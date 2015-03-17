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
    double Temperature = 0.00;
    double Ef = -6.0;
    double t = -1.5;
    double Miu1, Miu2, Miu4;
    double J_Hunds = -1.0;
    double J_exchange = -J_Hunds/100.0;
    double alpha = 0.1;
    double DM = J_exchange*4.0;
    int UpdateTorquePerStep = 100;
    bool CalculateTorque = false;
    double TimeStep = 0.01;
    Miu1 =  0.01;
    
    vec BackgroundField(3);
    BackgroundField << 0.0 << 0.0 << 5.0*J_exchange;
    SpinSystem SpinTexture("input.txt",J_exchange, DM, alpha);
    ElectronSystem Electrons("input.txt", "openBoundaries.txt",t);
    SpinTexture.NodeList[60].Pinned = false;
    if (CalculateTorque == true)
    {
       
        SpinTexture.ReadInElectronSiteIndex("ElectronSiteIndices.txt");
    }
    // Now set up the absorbing boundary
   /* for (int i=0; i<Electrons.ListOfOpenBoundaries.size(); i++)
    {
        for (int j=0; j<Electrons.ListOfOpenBoundaries[i].ListOfBoundarySites.size(); j++)
        {
            int Index = Electrons.ListOfOpenBoundaries[i].ListOfBoundarySites[j].SiteIndex;
            SpinTexture.NodeList[Index].alpha = 20.0;
        }
    }  // now the absorbing boundaries are done.*/
    vec ax, ay, az;
    ax << 11.0 << 0.0 << 0.0;
    ay << 0.0 << 31.0 << 0.0;
    az << 0.0 << 0.0 << 31.0;
    SpinTexture.GenerateNeighborList(1.00001, 1.00001, true, false, false, ax, ay, az);
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
    if (CalculateTorque == true)
        Electrons.CalculateGR(Ef);
    SpinTexture.SetBackgroundField(BackgroundField);
    SpinTexture.SetTemperature(Temperature);
    MovieWindow OutputWindow(0.0, 11.0, 0.0, 11.0, SpinTexture.NumSite);
    char Info[256];
    //OutputWindow.UpdateWindow(Haha, Info.str());
    int count = 0;
    for (double Time=0.0; Time<100000000000; Time += TimeStep)
    {
        if (count%UpdateTorquePerStep == 0 && CalculateTorque == true)
        {
            cx_mat T = Electrons.ThermalAverageTransmission(0.0, Ef);
            /*Complex P, Q;
            P = T(3,0)*T(1,0) + T(3,1)*T(1,0) + T(3,2)*T(1,0) + T(1,3)*T(3,0);
            Q = T(3,0)*T(1,2) + T(3,1)*T(1,2) + T(3,2)*T(1,2) + T(3,2)*T(1,3);
            Miu2 = real(P/(P+Q)- Q/(P+Q))*Miu1;
            P = T(1,0)*T(3,1) + T(1,0)*T(3,0) + T(1,2)*T(3,0) + T(1,3)*T(3,0);
            Q = T(1,0)*T(3,2) + T(1,2)*T(3,1) + T(1,2)*T(3,2) + T(1,3)*T(3,2);
            Miu4 = real(P/(P+Q) - Q/(P+Q))*Miu1;
            printf("Time=%lf\tMiu2=%lf\tMiu4=%lf\n", Time, Miu2, Miu4);
            Electrons.ListOfOpenBoundaries[0].ChemicalPotential = Miu1+Ef;
            Electrons.ListOfOpenBoundaries[1].ChemicalPotential = Miu2+Ef;
            Electrons.ListOfOpenBoundaries[2].ChemicalPotential = -Miu1+Ef;
            Electrons.ListOfOpenBoundaries[3].ChemicalPotential = Miu4+Ef;*/
        //Electrons.OutputElectronSpinMapProFit("debugOutput.txt", Ef);
            Electrons.ListOfOpenBoundaries[0].ChemicalPotential = Miu1+Ef;
            Electrons.ListOfOpenBoundaries[1].ChemicalPotential = -Miu1+Ef;
            SpinTexture.CalculateTorque(Electrons, Ef, J_Hunds);
        }
        SpinTexture.Evolve(TimeStep, Time, true);
        if (count%UpdateTorquePerStep == 0 && CalculateTorque == true)
        {
            Electrons.UpdateHamiltonian(SpinTexture, J_Hunds);
            Electrons.RenewGR(Ef);
        }
        count++;
        
        if (count%(UpdateTorquePerStep*10) == 0)
        {
            sprintf(Info, "t=%lf", Time);
            
            OutputWindow.UpdateWindow(SpinTexture, Info);
            Electrons.OnSiteSpin(60, Ef).print("On-site spin =");
        }
        if (fabs(Time-5000.0)<1.0e-5)
        {
            SpinTexture.OutputEffectiveToProFitTextFile("OutputEffectiveField.txt");
            SpinTexture.OutputTextureToTextFile("OutputTexture.txt");
        }
    }
}
    