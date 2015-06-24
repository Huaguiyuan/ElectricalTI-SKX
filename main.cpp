/* 
 * File:   main.cpp
 * Author: RAF_Jason
 *
 * Created on 2013年10月1日, 上午10:01
 */

#include "Calculations.hpp"
#include "MovieWindow.hpp"
#include <sstream>



void CalculateInjectionBandStructure(ElectronSystem &World)
{
    cx_mat F00 = World.ListOfOpenBoundaries[0].F00;
    cx_mat F01 = World.ListOfOpenBoundaries[0].F01;
    FILE* fp;
    fp = fopen("BandStructure.txt","w");
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



// The following code is used to test the 2-terminal case with a 21by11 ribbon.


int main (int argc, char** argv )
{
    system("rm ./*.gif -f");
    FILE* fpSkyrmionLocationX;
    FILE* fpSkyrmionLocationY;
    FILE* fpRealTimeCurrent;
    FILE* fpEnergyTrack;
    fpSkyrmionLocationX = fopen("SkyrmionLocationX.txt","w");
    fpSkyrmionLocationY = fopen("SkyrmionLocationY.txt","w");
    fpRealTimeCurrent = fopen("RealTimeCurrent.txt","w");
    fpEnergyTrack = fopen("EnergyTrack.txt","w");
    
    
    double Miu1, Miu2, Miu4;
    double J_Hunds = -1.0;
    double t = -0.2*fabs(J_Hunds);
    double Ef = -1.25*fabs(J_Hunds);
    double J_exchange = -J_Hunds/100.0;
    double Temperature = 0.0*J_exchange;
    double H0=6*J_exchange;
    double alpha = 0.1;
    double DM = 4.0*J_exchange;//*4.0;
    double I, Sx, Sy, Sz;
    int UpdateTorquePerStep = 100;
    bool CalculateTorque = false;
    double TimeStep = 0.01;
    Miu1 =  0.004;
    double FixedCurrent = 0.7e-6; //Unit is in A
    vec BackgroundField(3);
    BackgroundField << 0.0 << 0.0 << H0;
    SpinSystem SpinTexture("input.txt",J_exchange, DM, alpha);
    ElectronSystem Electrons("input.txt", "openBoundaries.txt", "BoundaryVirtualShift.txt", t, 1.00001);
    SpinTexture.NodeList[57].Pinned = false;
    double FM_Energy = -2.0*J_exchange*SpinTexture.NodeList.size()-SpinTexture.NodeList.size()*(H0-Temperature);
    if (CalculateTorque == true)
    {  
        SpinTexture.ReadInElectronSiteIndex("ElectronSiteIndices.txt");
    }
    vec ax, ay, az;
    ax << 21.0 << 0.0 << 0.0;
    ay << 0.0 << 31.0 << 0.0;
    az << 0.0 << 0.0 << 51.0;
    SpinTexture.GenerateNeighborList(1.0001, 1.00001, true, false, false, ax, ay, az);
    if (CalculateTorque == true)
    {
        cx_mat Trans(2,2);
        Electrons.CalculateGR(Ef);
        Trans = Electrons.ThermalAverageTransmission(0.0, Ef);
        
        Miu1 = FixedCurrent*4.135667516e-15/1.602176565e-19/real(Trans(0,1))/2.0;
    }
    printf("Miu1=%lf\n",Miu1);
    SpinTexture.SetBackgroundField(BackgroundField);
    SpinTexture.SetTemperature(Temperature);
    /*for (int i=0; i<SpinTexture.NodeList.size(); i++)
    {
        SpinTexture.NodeList[i].Temperature = SpinTexture.NodeList[i].Location(0)*(0.9/99.0)+0.0;
    } This is for the temperature gradient*/
    MovieWindow OutputWindow(0.0, 22, 0.0, 11.0, SpinTexture.NumSite);
    char Info[256];
    int count = 0;
    vec oldLocation(3);
    vec newLocation(3);
    int CrossTimesX = 0;
    int CrossTimesY = 0;
    oldLocation = SpinTexture.SkyrmionLocation();
    for (double Time=0.0; Time<100000000000; Time += TimeStep)
    {
        /*if (Time > 100)
            H0 = 2.0/1000.0*Time; This is for increasing external field */
        BackgroundField(2) = H0*cos(0.0);
        BackgroundField(0) = H0*sin(0.0);
        SpinTexture.UpdateBackgroundField(BackgroundField);
        FM_Energy = -2.0*J_exchange*SpinTexture.NodeList.size()-SpinTexture.NodeList.size()*(H0-Temperature);
        
        if (count%UpdateTorquePerStep == 0 && CalculateTorque == true)
        {
            cx_mat T = Electrons.ThermalAverageTransmission(0.0, Ef);
            Electrons.ListOfOpenBoundaries[0].ChemicalPotential = Miu1+Ef;
            Electrons.ListOfOpenBoundaries[1].ChemicalPotential = -Miu1+Ef;
            SpinTexture.CalculateTorque(Electrons, Ef, J_Hunds);
            if (count == 0)
            {
               // SpinTexture.OutputTorqueFieldToProFitTextFile("TorqueField.txt");
                //SpinTexture.OutputSpinTextureGIF(0,11,0,11,"")
            }
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
            double TotalEnergy;
            TotalEnergy = SpinTexture.TotalEnergy(true);
            SpinTexture.RenormalizeLength();
            double TQ;
            TQ = TopologicalCharge(SpinTexture);
            sprintf(Info, "t=%lf  E=%lf TQ=%lf H0=%lf", Time, (TotalEnergy-FM_Energy)/J_exchange, TQ, H0/fabs(J_exchange));
            fprintf(fpEnergyTrack, "%lf\t%lf\t%1.2lf\t%lf\n", Time, (TotalEnergy-FM_Energy)/J_exchange, TQ, H0/fabs(J_exchange));
            fflush(fpEnergyTrack);
            OutputWindow.UpdateWindow(SpinTexture, Info);
            newLocation = SpinTexture.SkyrmionLocation();
            vec Change = newLocation - oldLocation;
            //Electrons.CalculateTerminalSpinCurrent(Electrons.ListOfOpenBoundaries[0], Ef,I, Sx, Sy, Sz);
            if (fabs(Change(0)) > 0.9 || count/(UpdateTorquePerStep*10)==1)
            {
                //Electrons.CalculateTerminalSpinCurrent(Electrons.ListOfOpenBoundaries[0], Ef,I, Sx, Sy, Sz);
               // printf("% le\t% le\n", Time, Sz);
                if (Change(0) > 1.1)
                    CrossTimesX--;
                if (Change(0) < -1.1)
                    CrossTimesX++;
                
                fprintf(fpSkyrmionLocationX, "%le\t%le\t%le\t%le\n", Time, newLocation(0)+(double)CrossTimesX*ax(0), newLocation(1), newLocation(2));
                fflush(fpSkyrmionLocationX);
                //printf("%le\t%le\t%le\t%le\n", Time, newLocation(0), newLocation(1), newLocation(2));
                oldLocation = newLocation;
                if (CalculateTorque == true)
                {
                    Electrons.CalculateTerminalSpinCurrent(Electrons.ListOfOpenBoundaries[0], Ef,I, Sx, Sy, Sz);
                    fprintf(fpRealTimeCurrent, "%lf\t% le\t% le\n", Time, I*1.602176565e-19, Sz*1.602176565e-19);
                    fflush(fpRealTimeCurrent);
                }
            }
            if (fabs(Change(1)) > 0.9)
            {
                if (Change(1) > 1.1)
                    CrossTimesY--;
                if (Change(1) < -1.1)
                    CrossTimesY++;
                
                fprintf(fpSkyrmionLocationY, "%le\t%le\t%le\t%le\n", Time, newLocation(0), newLocation(1)+(double)CrossTimesY*ay(1), newLocation(2));
                fflush(fpSkyrmionLocationY);
                //printf("%le\t%le\t%le\t%le\n", Time, newLocation(0), newLocation(1), newLocation(2));
                oldLocation = newLocation;
                if (CalculateTorque == true)
                {
                    Electrons.CalculateTerminalSpinCurrent(Electrons.ListOfOpenBoundaries[0], Ef,I, Sx, Sy, Sz);
                    fprintf(fpRealTimeCurrent, "%lf\t% le\t% le\n", Time, I*1.602176565e-19, Sz*1.602176565e-19);
                    fflush(fpRealTimeCurrent);
                }
            }
            if (count%(UpdateTorquePerStep*100) == 0)
            {
                char Buffer[256];
                sprintf(Buffer, "Time = %.2lf H0=%lf", Time, H0/fabs(J_exchange));
                SpinTexture.OutputSpinTextureGIF(0.0, 22.0, 0.0, 11.0, Buffer);
            }
        }
        if (fabs(Time-1000.0)<1.0e-5 )//&& CalculateTorque == true)
        {
            //Electrons.UpdateHamiltonian(SpinTexture, J_Hunds);
            //Electrons.RenewGR(Ef);
            SpinTexture.CalculateTorque(Electrons, Ef, J_Hunds);
            SpinTexture.CalculateEffectiveField(true);
            SpinTexture.OutputOtherTorqueProFit("OutputOtherTorque.txt");
            SpinTexture.OutputEffectiveToProFitTextFile("OutputEffectiveField.txt");
            SpinTexture.OutputTextureToTextFile("OutputTexture.txt");
            //SpinTexture.OutputTorqueFieldToProFitTextFile("OutputTorqueField.txt");
            SpinTexture.OutputTextureToProFeitTextFile("OutputTextureProFit.txt");
           // Electrons.CalculateCurrentDistribution(Ef);
            SpinTexture.OutputSTTProFit("OutputTorqueField.txt");
            //Electrons.OutputSpinCurrentMapProFit("CurrentMapProFit.txt", 0, 21, 0, 11);
            
            if(CalculateTorque == true)
            {
                //Miu1 = 0.0;
                Electrons.ListOfOpenBoundaries[0].CalculateBandStructure("BandStructure.txt");
                //Electrons.OutputSpinCurrentMapProFit("OutputSpinCurrent.txt", 0, 10, 0, 10);
            } 
        }
        if (Time > 100000.0)
            break;
    }
    fclose(fpSkyrmionLocationX);
    fclose(fpSkyrmionLocationY);
    fclose(fpRealTimeCurrent);
    fclose(fpEnergyTrack);
}
    


/*  
 // The following is the code to calculate the NEGF_LLG in the cross-bar case.
int main (int argc, char** argv )
{
    system("rm ./*.gif -f");
    FILE* fpSkyrmionLocation;
    fpSkyrmionLocation = fopen("SkyrmionLocation.txt","w");
    double Temperature = 0.00;
    double Ef = 6.0;
    double t = -1.5;
    double Miu1, Miu2, Miu4;
    double J_Hunds = -1.0;
    double J_exchange = -J_Hunds/100.0;
    double alpha = 0.1;
    double DM = J_exchange*4.0;
    int UpdateTorquePerStep = 1000;
    bool CalculateTorque = true;
    double TimeStep = 0.01;
    Miu1 =  0.01;
    
    vec BackgroundField(3);
    BackgroundField << 0.0 << 0.0 << 5.0*J_exchange;
    SpinSystem SpinTexture("input.txt",J_exchange, DM, alpha);
    ElectronSystem Electrons("input.txt", "openBoundaries.txt", "BoundaryVirtualShift.txt", t, 1.000001);
    SpinTexture.NodeList[60].Pinned = false;
    if (CalculateTorque == true)
    {
       
        SpinTexture.ReadInElectronSiteIndex("ElectronSiteIndices.txt");
    }

    vec ax, ay, az;
    ax << 31.0 << 0.0 << 0.0;
    ay << 0.0 << 31.0 << 0.0;
    az << 0.0 << 0.0 << 51.0;
    SpinTexture.GenerateNeighborList(1.00001, 1.00001, true, true, false, ax, ay, az);
    if (CalculateTorque == true)
        Electrons.CalculateGR(Ef);
    SpinTexture.SetBackgroundField(BackgroundField);
    SpinTexture.SetTemperature(Temperature);
    MovieWindow OutputWindow(0.0, 31, 0.0, 31.0, SpinTexture.NumSite);
    char Info[256];
    //OutputWindow.UpdateWindow(Haha, Info.str());
    int count = 0;
    for (double Time=0.0; Time<100000000000; Time += TimeStep)
    {
        if (count%UpdateTorquePerStep == 0 && CalculateTorque == true)
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
            Electrons.ListOfOpenBoundaries[0].ChemicalPotential = Miu1+Ef;
            Electrons.ListOfOpenBoundaries[1].ChemicalPotential = Miu2+Ef;
            Electrons.ListOfOpenBoundaries[2].ChemicalPotential = -Miu1+Ef;
            Electrons.ListOfOpenBoundaries[3].ChemicalPotential = Miu4+Ef;
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
            vec Location = SpinTexture.SkyrmionLocation();
            fprintf(fpSkyrmionLocation, "%le\t%le\t%le\t%le\n", Time, Location(0), Location(1), Location(2));
            fflush(fpSkyrmionLocation);
            //printf("%le\t%le\t%le\t%le\n", Time, Location(0), Location(1), Location(2));
            SpinTexture.OutputSpinTextureGIF(0.0, 31.0, 0.0, 31.0, "J=1, D=4, H0=5");
           // Electrons.OnSiteSpin(60, Ef).print("On-site spin =");
        }
        if (fabs(Time-5000.0)<1.0e-5)
        {
            SpinTexture.OutputEffectiveToProFitTextFile("OutputEffectiveField.txt");
            SpinTexture.OutputTextureToTextFile("OutputTexture.txt");
            
        }
    }
    fclose(fpSkyrmionLocation);
}
*/