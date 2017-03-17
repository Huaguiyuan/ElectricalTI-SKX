
/* 
 * File:   main.cpp
 * Author: RAF_Jason
 *
 * Created on 2013年10月1日, 上午10:01
 */
#include "Calculations.hpp"
#include "SpinSystem.hpp"


int SkyrmionCalculation(void)
{
    // Parameter Setup
    FILE* fp;
    FILE* test;
    fp = fopen("log","w");
    int SizeX = 50;
    int SizeY = 50;
    int SizeZ = 3;
    int PlotY=SizeY;
    int PlotZ=SizeZ;
    double J = 1.0;
    double D = 0.3;//0.059697; 
    double Radius = 12;
    double CutOff = 20;
    double CenterX = 25;  // Center of the AC field one
    double CenterY = 25;  // Center of the AC field one
    int StepLimit = 999999999;   
    vec h0(3);
    vec h_AC(3);
    h_AC.zeros();
    h0.zeros();
    h0(2) = 0.09;//0.0036;
    h_AC(1) = 0.0;//.225*0.0036;//.225*h0(2);//0.09*0.0;
    double temperature = 0.0;
    double TimeStep = 0.05;
    double alpha = 0.1;
    double rate = 4.0e-5; // unit is in J^2
    double TQ;
    double Energy;
    double FM_Energy = -(double)SizeX*(double)SizeY*(double)(SizeZ-2)*J*2
                       -h0(2)*(double)SizeX*(double)SizeY*(double)(SizeZ-2)
                       -(double)SizeX*(double)SizeY*(double)(SizeZ-3)*J;
    // Now allocate the memory
    vec*** Plate;
    vec*** ExternalField;
    vec*** EffectiveField;
    double*** TopologicalChargeDistribution;
    double*** EnergyDistribution;
    double*** Heisenberg;
    double*** DM;
    double*** External;
    NeighborRecord*** Neighbors;
    float* xp;
    float* yp;  
    float* xv;
    float* yv;
    float* xp1;
    float* yp1;
    float* xv1;
    float* yv1;
    int Created=-1;
    int MomentRecord = 0;
    //Dislin SpinOutput;
    Dislin SpinGraphTop;
    Plate = Allocate_3D_Vector_Plate(SizeX, SizeY, SizeZ);
    ExternalField = Allocate_3D_Vector_Plate(SizeX, SizeY, SizeZ);
    EffectiveField = Allocate_3D_Vector_Plate(SizeX, SizeY, SizeZ);
    TopologicalChargeDistribution = Allocate_3D_Double_Plate(SizeX, SizeY, SizeZ);
    EnergyDistribution = Allocate_3D_Double_Plate(SizeX, SizeY, SizeZ);
    Heisenberg = Allocate_3D_Double_Plate(SizeX, SizeY, SizeZ);
    DM = Allocate_3D_Double_Plate(SizeX, SizeY, SizeZ);
    External = Allocate_3D_Double_Plate(SizeX, SizeY, SizeZ);
    Neighbors = GenerateNeighbors(SizeX, SizeY, SizeZ);
    /// Now get a seed for the random number generator
    int seed = (int)time(0);
    CRandomMersenne RanGen(seed);
    /// Now initialize everything
    //InitializeSpinPlot(SpinGraphBottom, xp, yp, xv, yv, SizeX, SizeY);
    //InitializeSpinPlot(SpinGraphTop, xp, yp, xv, yv, SizeX, SizeY);
    InitializeSpinPlotYZ(SpinGraphTop, xp, yp, xv, yv, SizeY, SizeZ);
    //InitializeSpinPlotYZ(SpinOutput, xp, yp, xv, yv, SizeY, SizeZ);
    InitializeFM_With_Perturbation(Plate, SizeX, SizeY, SizeZ, 0.0);
    InitializeExternalField(ExternalField, SizeX, SizeY, SizeZ, h0);
    /// Now start the simulation
    FILE *fpTQalongZ;
    fpTQalongZ = fopen("TQalongZ.txt","w");
    FILE *fpTQ_R;
    fpTQ_R = fopen("TQ_R.txt","w");
    FILE *fpSpinTrack;
    fpSpinTrack = fopen("SpinTrack.txt","w");
    int StartCount = 0;
    for (int count=0; count<StepLimit; count++)
    {
        double time = (double)count*TimeStep;
        Evolve_PeriodicBoundary_FiniteTemperature(RanGen, temperature, Plate, 
                ExternalField, Neighbors, rate, h_AC, h0, Radius, CutOff, 
                CenterX, CenterY, J, D, time, TimeStep, alpha, SizeX, SizeY, 
                SizeZ);
        //sleep(1);
        SpinNormalization(Plate, SizeX, SizeY, SizeZ);
        TQ = CalculateTopologicalCharge(Plate, TopologicalChargeDistribution, Neighbors, SizeX, SizeY, SizeZ);
        
        if (count % 20 == 0)
        {
            Energy = TotalEnergy(Plate, J, D, SizeX, SizeY, SizeZ, Neighbors, ExternalField, EnergyDistribution, Heisenberg, DM, External);
           // fprintf(fpSpinTrack, "%lf\t% lf\n", time, Energy-FM_Energy);
          //  UpdateSpinPlot(SpinGraphBottom, Plate, rate, count, TimeStep, TQ, xp, yp, xv, yv, SizeX, SizeY, 2, Energy);
            //UpdateSpinPlot(SpinGraphTop, Plate, rate, count, TimeStep, TQ, xp, yp, xv, yv, SizeX, SizeY, 1, Energy-FM_Energy);
            UpdateSpinPlotYZ(SpinGraphTop, Plate, rate, count, TimeStep, TQ, xp, yp, xv, yv, SizeY, SizeZ, 25, Energy-FM_Energy);
            //WriteTextureToFileYZ(Plate, SizeY, SizeZ, 25, PlotY, PlotZ, TQ);
            fprintf(fp, "time= % lf\tHc = % lf\tTQ = % lf\n", time, time*rate, TQ);
            fflush(fp);
           // fflush(fpSpinTrack);
        }
  
        
        if ( fabs(TQ)>(SizeZ-3+0.1) )
        {
            if (rate*time>0.0)
                rate = -4e-5 + 2.0*4e-5*3550./time;
            else
                rate = 0.;
            StartCount++;
            if (count % 100 == 0)
                WriteTextureToFileYZ(Plate, SizeY, SizeZ, 25, PlotY, PlotZ, TQ);
            fprintf(fpTQalongZ, "%lf\t", time);
            for (int k=1; k<SizeZ-1; k++)
            {
                double QLayer = 0.0;
                for (int i=0; i<SizeX; i++)
                {
                    for (int j=0; j<SizeY; j++)
                    {
                        QLayer += TopologicalChargeDistribution[i][j][k];
                    }
                }
                fprintf(fpTQalongZ, "%lf\t", fabs(QLayer));
                
                fflush(fpTQalongZ);
            }
            fprintf(fpTQalongZ, "\n");
            if ( StartCount > 8000 )
            {
                
                break;
            }
        }
        if (time > 3000 && time < 4000)
        {
            vec temp = Plate[25][25][1];
            //if (count%20 == 0)
                //WriteTextureToFileYZ(Plate, SizeY, SizeZ, 25, PlotY, PlotZ, TQ);
            Energy = TotalEnergy(Plate, J, D, SizeX, SizeY, SizeZ, Neighbors, ExternalField, EnergyDistribution, Heisenberg, DM, External);
            fprintf(fpSpinTrack, "%lf\t% lf\t% lf\t% lf\t% lf\n", time, Energy-FM_Energy, temp(0), temp(1), temp(2));
            fprintf(fp, "time= % lf\tHc = % lf\tTQ = % lf\n", time, time*rate, TQ);
            for (int i=0; i<SizeX; i++)
            {
                fprintf(fpTQ_R, "% lf\t% d\t% lf\n",time, i, TopologicalChargeDistribution[i][25][1]);
            }
            
        }
        //if (time > 4000)
           // break;
        /*if (fabs(TQ) > 0.9)
        {
            Created = 1;
            MomentRecord = count;
               //break;
        }
        if (Created > 0 && count > 1000 + MomentRecord)
        {
            fprintf(fp, "time= % lf\tHc = % lf\tTQ = % lf\n", TimeStep*MomentRecord, TimeStep*MomentRecord*rate, TQ);
            fclose(fp);
            break;
        }    */
    }
    fclose(fpTQalongZ);
    fclose(fpSpinTrack);
    fclose(fp);
    fclose(fpTQ_R);
    Delete_3D_Plate(Plate, SizeX, SizeY);
    FinalizeSpinPlot(SpinGraphTop, xp, yp, xv, yv);
    //FinalizeSpinPlot(SpinGraphBottom, xp, yp, xv, yv);
    //FinalizeSpinPlot(SpinGraphTop, xp, yp, xv, yv);
    DeleteNeighbors(Neighbors, SizeX, SizeY);
    
    Delete_3D_Plate(ExternalField, SizeX, SizeY);
    Delete_3D_Plate(EffectiveField, SizeX, SizeY);
    Delete_3D_Double_Plate(TopologicalChargeDistribution, SizeX, SizeY);
    Delete_3D_Double_Plate(EnergyDistribution, SizeX, SizeY);
    Delete_3D_Double_Plate(Heisenberg, SizeX, SizeY);
    Delete_3D_Double_Plate(DM, SizeX, SizeY);
    Delete_3D_Double_Plate(External, SizeX, SizeY);
    
    
    return 1;
    
    
    
}
/*
int main(void)
{
    //SkyrmionCalculation();
    int SizeX=3;
    int SizeY=3;
    int SizeZ=1;
    double ExchangeCutoff = 1.1;
    double StrayCutoff = 99999.;
    double J = 1.0;
    double D=0.3;
    SpinSystem World(SizeX, SizeY, SizeZ, ExchangeCutoff, StrayCutoff, J, D);
    World.CalculateEffectiveField();
    World.EffectiveField[8].print();
    for (int i=0; i<SizeX*SizeY*SizeZ; i++)
    {
        int p=World.NodeList[i].ListOfExchangeNeighbours.size();
        printf("%d\n", p);
    }
    putchar('\n');
    for (int i=0; i<World.EffectiveField.size(); i++)
    {
        World.EffectiveField[i].print();
    }        
}
*/

/*
int main(int argc, char** argv) {
    // Parameter Definition 
    int SizeX = 50;
    int SizeY = 50;
    double J = 1.0;
    double D = 0.3;//0.059697; 
    double Radius = 12;
    double CutOff = 20;
    double Center1X = 25;  // Center of the AC field one
    double Center1Y = 25;  // Center of the AC field one
    double Center2X = 25;  // Center of the AC field two
    double Center2Y = 25;  // Center of the AC field two
    int PlotX = 40;
    int PlotY = 40;
    int StepLimit = 999999999;   

    int N_Record = 1000;
    vec h0(3);
    vec h_AC(3);
    h_AC.zeros(3);
    h0.zeros(3);
    
    //h_AC(2) = 0.0;
    h0(2) = 0.09;//0.0036;
    h0(1) = 0.00;
    h_AC(1) = 0.1;//.225*0.0036;//.225*h0(2);//0.09*0.0;
    double Period = 10000;
    double omega = 2.0*PIPI/Period;
    //double omega = 0.002;
    double temperature = 0.0;
    double TimeStep = 0.3;//1. * 2.0*PIPI/omega/(Period*10.);
    double FFT_SamplingTimeStep = 2.;
    int StepInteger = (int)(FFT_SamplingTimeStep/TimeStep);
    double alpha = 0.1;
    int N_for_average=1;
    cx_vec SpinRecord(N_Record);
    
    // Now allocate the memory;
    vec** Plate;
    vec** ExternalField;
    vec** EffectiveField;
    double** TopologicalChargeDistribution;
    double** EnergyDistribution;
    double** Heisenberg;
    double** DM;
    double** External;
    NeighborRecord** Neighbors;
    float* xp;
    float* yp;  
    float* xv;
    float* yv;
    float* xp1;
    float* yp1;
    float* xv1;
    float* yv1;
    float* ChargePlotPointer;
    Dislin SpinGraph;
    Dislin ChargeGraph;
    Plate = Allocate_Vector_Plate(SizeX, SizeY);
    ExternalField = Allocate_Vector_Plate(SizeX, SizeY);
    Neighbors = GenerateNeighbors(SizeX, SizeY);
    EffectiveField = Allocate_Vector_Plate(SizeX, SizeY);
    FILE* fp;
    fp = fopen("output.txt","w");
    fclose(fp);
    InitializeSpinPlot(SpinGraph, xp, yp, xv, yv, SizeX, SizeY);
    TopologicalChargeDistribution = AllocateDoublePlate(SizeX, SizeY);
    EnergyDistribution = AllocateDoublePlate(SizeX, SizeY);
    Heisenberg = AllocateDoublePlate(SizeX, SizeY);
    DM = AllocateDoublePlate(SizeX, SizeY);
    External = AllocateDoublePlate(SizeX, SizeY);
    
    int seed = (int)time(0);
    CRandomMersenne RanGen(seed);
    
   
    EnforcedSpin* EnforcedList;
    int N_Enforced;
         fp = fopen("output.txt","a");
        fprintf(fp, "%le", h_AC(1));
        fclose(fp);
        double AverageTime=0.0;
        double TQ_previous=0.0;
        double TQ_now;
        double FM_Energy = -(double)SizeX*(double)SizeY*J*2-h0(2)*(double)SizeX*(double)SizeY;
        FILE* fpEnergy;
        fpEnergy = fopen("energy.txt","w");
        FILE* fpCriticalField;
        fpCriticalField = fopen("CriticalField.txt","w");
        fclose(fpCriticalField);

                printf("Current h_AC is %lf\n", h_AC(1));
                InitializeExternalField(ExternalField, SizeX, SizeY, h0);
                //EnforceExternalFieldOnSite(0, 0, 0.1, ExternalField);
                InitializeFM_With_Perturbation(Plate, SizeX, SizeY, 0.0);
                //ReadInTexture(Plate, SizeX, SizeY);
                //InitializeDotSpinWaveExcitation(Plate, 25, 25, 0.0, 1.0, 0.0);
                
                
                //WriteTextureToFile(Plate, SizeX, SizeY);
        
        //InitializeChargePlot(ChargeGraph, ChargePlotPointer ,SizeX, SizeY, PlotX, PlotY);
                
    //InitializePlot(xp, yp, xv, yv, SizeX, SizeY);
   //int count = 0;
                
                double time = 0.0;
                int FFT_Count = 0;
                int FieldFound = 0;
                double TQ;
                double Energy;
                double TotalMinusAC;
                double E_AC;
                double FieldHeisenberg, FieldDM, FieldExtra, FieldTotal, FieldTotalXY, DMXY;
                
                for (int count=0; count<=StepLimit; count++)
                {
                        time = (double)count*TimeStep;
                        if (time > 2000)
                        {
                            //WriteToTextFiles(Plate, SizeX, SizeY, "FinalState.txt");
                            break;
                        }
                        //if (time > Period+1000.)
                        //{
                            //WriteSzToTextFilex(Plate, SizeX, SizeY, "SzFinal.txt");
                            //WriteTheEnergyAlongX(EnergyDistribution, Heisenberg, DM, External, SizeX, 49, "EnergyAlongXFinal.txt", h0, J);
                            //WriteEffectiveFieldAlongX(Plate, Neighbors, ExternalField, SizeX, 49, "HeffAlongXFinal.txt", h0, J, D);
                            //WriteTextureToFile(Plate, SizeX, SizeY, PlotX, PlotY);
                          //  break;
                        //}
                       
       RungeKuttaEvolve_PeriodicBoundary(Plate, ExternalField, Neighbors, omega, h_AC, h0, Radius, CutOff, Center1X, Center1Y, Center2X, Center2Y,
                                         J, D, time, TimeStep, alpha, SizeX, SizeY, EnforcedList, N_Enforced); 
            
                        //Evolve_PeriodicBoundary_FiniteTemperature(RanGen, temperature, Plate, ExternalField, Neighbors, omega, h_AC, h0, Radius, CutOff, Center1X,
                        //Center1Y, Center2X, Center2Y, J, D, time, TimeStep, alpha, SizeX, SizeY);
            
       //CalculateEffectiveField(J, D, Plate, Neighbors, ExternalField, EffectiveField, SizeX, SizeY);
                        SpinNormalization(Plate, SizeX, SizeY);
                        
                        //TQ = CalculateTopologicalCharge(Plate, TopologicalChargeDistribution, Neighbors, SizeX, SizeY);
                        //TQ_now = TQ;
                        //if ( abs(abs(TQ_now-TQ_previous)-1.0)<1.0e-2)
                         //   printf("time is %lf\n", time);
                        //TQ_previous = TQ_now;
                        
                        FieldTotal = OnSiteEffectiveFieldZ(49, 49, Plate, Neighbors, ExternalField, J, D, FieldHeisenberg, FieldDM, FieldExtra, FieldTotalXY, DMXY);
                        //fpCriticalField = fopen("CriticalField.txt","a");
                        //fprintf(fpCriticalField, "% lf\t% lf\t% lf\t% lf\t% lf\t% lf\n", 2.0*omega/PIPI*time*h_AC(1), FieldHeisenberg+FieldExtra, FieldDM, FieldTotal, FieldTotalXY, DMXY);
                        //fclose(fpCriticalField);
                        
                        
                        if (count%10 == 0)
                        {
                            
                                TQ = CalculateTopologicalCharge(Plate, TopologicalChargeDistribution, Neighbors, SizeX, SizeY);
                                UpdateSpinPlot(SpinGraph, Plate, omega, count, TimeStep, TQ, xp, yp, xv, yv, SizeX, SizeY, Energy-FM_Energy);
                                //OutputOnSiteData(Plate, count, 24, 24, 25, 24, 24, 25, 25, 25, 25, 25, 373, 374, fp);
                                Energy = -99999;//TotalEnergy(Plate, J, D, SizeX, SizeY, Neighbors, ExternalField, EnergyDistribution, Heisenberg, DM, External, E_AC, TotalMinusAC);
                                //fprintf(fpEnergy, "%lf\t% le\t% le\t% le\n", time, Energy-FM_Energy, TotalMinusAC-FM_Energy, E_AC);
                                //WriteTextureToFile(Plate, SizeX, SizeY, PlotX, PlotY, TQ);
                                //fflush(fpEnergy);
           
                        }

                        
                }
//                fprintf(fpEnergy, "%lf\t% le\n", (double)Distance, Energy-FM_Energy);
                //fflush(fpEnergy);
                //fclose(fpEnergy);
                //if (TQ < -0.8)
                   // break;
                WriteToTextFiles(Plate, SizeX, SizeY, "FinalState.txt");
       // }
        fclose(fpEnergy);
        // }      
        //}
                
        
        
    
    fclose(fp);
    //cx_vec FFT_Result(N_Record);
    //printf("%lf\n", real(SpinRecord(0)));
    //FFT_Result = FFT_SamplingTimeStep*fft(SpinRecord);
    //fp = fopen("FFT.txt","w");
    //for (int i=0; i<N_Record; i++)
    //{
      //  double CircularFrequency = (double)i*2.0*PIPI/(double)N_Record/FFT_SamplingTimeStep;
       // fprintf(fp, "% lf\t% lf\t% lf\t% lf\n", CircularFrequency, real(FFT_Result(i)), imag(FFT_Result(i)), abs(FFT_Result(i)));
    //}
    //fclose(fp);
    //fp = fopen("TimeDomain.txt","w");
    //for (int i=0; i<N_Record; i++)
    //{
    //        fprintf(fp, "% lf\t% lf\n", (double)i*FFT_SamplingTimeStep, real(SpinRecord(i)));
    //}
    //fclose(fp);
    WriteToTextFiles(Plate, SizeX, SizeY, "FinalState.txt");
    // getchar();
    
    //WriteTheEnergyAlongX(EnergyDistribution, Heisenberg, DM, External, SizeX, 24, "EnergyAlongX.txt", h0, J);
    FinalizeSpinPlot(SpinGraph, xp, yp, xv, yv);
    FinalizeChargePlot(ChargeGraph, ChargePlotPointer);
    Delete_Plate(Plate, SizeX);
    Delete_Plate(ExternalField, SizeX);
    Delete_Plate(EffectiveField, SizeX);
    DeleteDoublePlate(EnergyDistribution, SizeX);
    DeleteDoublePlate(Heisenberg, SizeX);
    DeleteDoublePlate(DM, SizeX);
    DeleteDoublePlate(External, SizeX);
    delete [] EnforcedList;
    DeleteDoublePlate(TopologicalChargeDistribution, SizeX);
    
    return 0;
    
}
///////////////////////////////////////////////////////////////////////

*/

void CalculateCurrentDistribution(double Ef, double Miu1, double Miu2, double Miu4, ElectronSystem &World)
{
    
    
    
    printf("Starting to calculate Current OP\n");
    cx_mat SigmaX, SigmaY, SigmaZ;
    SigmaX << 0.0 << 1.0 << endr
           << 1.0 << 0.0 << endr;
    SigmaY << 0.0 << Complex(0.,-1.) << endr
           << Complex(0.,1.) << 0.0 << endr;
    SigmaZ << 1.0 << 0.0 << endr
           << 0.0 <<-1.0 << endr;   
    double MaxCurrent, MaxSxCurrent, MaxSyCurrent, MaxSzCurrent;
    MaxCurrent = MaxSxCurrent = MaxSyCurrent = MaxSzCurrent = 0.0;
    cx_mat a1, a2, a3, a4;
    a1.set_size(2,2);
    a2.set_size(2,2);
    a3.set_size(2,2);
    a4.set_size(2,2);
    for (int i=0; i<World.NumSite; i++)
    {
        vec Current(3);
        vec SpinXCurrent(3);
        vec SpinYCurrent(3);
        vec SpinZCurrent(3);
        Current.zeros();
        SpinXCurrent.zeros();
        SpinYCurrent.zeros();
        SpinZCurrent.zeros();
        if (World.ListOfSites[i].IsBoundary == true)
        {
            World.ListOfSites[i].Current = Current;
            World.ListOfSites[i].SpinXCurrent = SpinXCurrent;
            World.ListOfSites[i].SpinYCurrent = SpinYCurrent;
            World.ListOfSites[i].SpinZCurrent = SpinZCurrent;
            continue;
        }
         
        
            for (int j=0; j<World.ListOfSites[i].ListOfTightBindingNeighbors.size(); j++)
            {
                int SourceIndex = World.ListOfSites[i].ListOfTightBindingNeighbors[j];
                int I, J;
                I = World.ListOfSites[i].SiteIndex;
                J = World.ListOfSites[SourceIndex].SiteIndex;
                a1 = World.OnSiteSpectralFunctionFromBoundary(I, J, World.ListOfOpenBoundaries[0]);
                a2 = World.OnSiteSpectralFunctionFromBoundary(I, J, World.ListOfOpenBoundaries[1]);
                a3 = World.OnSiteSpectralFunctionFromBoundary(I, J, World.ListOfOpenBoundaries[2]);
                a4 = World.OnSiteSpectralFunctionFromBoundary(I, J, World.ListOfOpenBoundaries[3]);
                cx_mat I_local = -Complex(0.0, 1.0)*(Miu1*(a1-trans(a1)) + Miu2*(a2-trans(a2)) + Miu4*(a4-trans(a4)) - Miu1*(a3-trans(a3)));
                Current = Current + real(trace(I_local))*(World.ListOfSites[i].Location-World.ListOfSites[SourceIndex].Location);
                SpinXCurrent = SpinXCurrent + real(trace(SigmaX*I_local))*(World.ListOfSites[i].Location-World.ListOfSites[SourceIndex].Location);
                SpinYCurrent = SpinYCurrent + real(trace(SigmaY*I_local))*(World.ListOfSites[i].Location-World.ListOfSites[SourceIndex].Location);
                SpinZCurrent = SpinZCurrent + real(trace(SigmaZ*I_local))*(World.ListOfSites[i].Location-World.ListOfSites[SourceIndex].Location);
            }
            World.ListOfSites[i].Current = Current;
            World.ListOfSites[i].SpinXCurrent = SpinXCurrent;
            World.ListOfSites[i].SpinYCurrent = SpinYCurrent;
            World.ListOfSites[i].SpinZCurrent = SpinZCurrent; 
            double Norm = sqrt(dot(Current, Current));
            if (Norm > MaxCurrent)
                MaxCurrent = Norm;
            Norm = sqrt(dot(SpinXCurrent, SpinXCurrent));
            if (Norm > MaxSxCurrent)
                MaxSxCurrent = Norm;
            Norm = sqrt(dot(SpinYCurrent, SpinYCurrent));
            if (Norm > MaxSyCurrent)
                MaxSyCurrent = Norm;
            Norm = sqrt(dot(SpinZCurrent, SpinZCurrent));
            if (Norm > MaxSzCurrent)
                MaxSzCurrent = Norm;
    }
    printf("Current calculation done.\n");
    for (int i=0; i<World.NumSite; i++)
    {
        World.ListOfSites[i].Current = World.ListOfSites[i].Current/MaxCurrent;
        World.ListOfSites[i].SpinXCurrent = World.ListOfSites[i].SpinXCurrent/MaxSxCurrent;
        World.ListOfSites[i].SpinYCurrent = World.ListOfSites[i].SpinYCurrent/MaxSyCurrent;
        World.ListOfSites[i].SpinZCurrent = World.ListOfSites[i].SpinZCurrent/MaxSzCurrent;  
    }
    FILE* fp;
    fp = fopen("LocalMiu.txt","w");
    for (int i=0; i<World.NumSite; i++)
    {
        Complex delta = trace(World.OnSiteSpectralFunctionFromBoundary(i,i,World.ListOfOpenBoundaries[0])*Miu1 +
                World.OnSiteSpectralFunctionFromBoundary(i,i,World.ListOfOpenBoundaries[1])*Miu2 + 
                World.OnSiteSpectralFunctionFromBoundary(i,i,World.ListOfOpenBoundaries[3])*Miu4 - 
                World.OnSiteSpectralFunctionFromBoundary(i,i,World.ListOfOpenBoundaries[2])*Miu1);
        Complex Norm = trace(World.OnSiteSpectralFunctionFromBoundary(i,i,World.ListOfOpenBoundaries[0]) +
                World.OnSiteSpectralFunctionFromBoundary(i,i,World.ListOfOpenBoundaries[1]) + 
                World.OnSiteSpectralFunctionFromBoundary(i,i,World.ListOfOpenBoundaries[3]) + 
                World.OnSiteSpectralFunctionFromBoundary(i,i,World.ListOfOpenBoundaries[2]));
        fprintf(fp, "%lf\t%lf\t%le\n", World.ListOfSites[i].Location(0), World.ListOfSites[i].Location(1), real(delta/Norm));
                
    }
    fclose(fp);
    FILE* fp2;
    
    fp2 = fopen("LocalSz.txt","w");
    printf("Calculating Sz Map.\n");
    for (int i=0; i<World.NumSite; i++)
    {
        cx_mat GN(2,2);
        GN = World.OnSiteCorelationFunctionGn(i,i, Ef);
        Complex sz;
        sz = trace(SigmaZ*GN);
        fprintf(fp2, "% lf\t% lf\t% le\n", World.ListOfSites[i].Location(0), 
                                          World.ListOfSites[i].Location(1),
                                          real(sz));
    }
    fclose(fp2);
    
}    

void CalculateHallEffect(double Temperature, double Ef, double Miu1, 
                         double &Miu2, double &Miu4, double &totalI, double &totalSz,
                         double &UpDos, double &DownDos)
{
    ElectronSystem World("input.txt","openBoundaries.txt","BoundaryVirtualShift.txt", -1.5, 1.0000001);
    //World.PrintBoundaryList();
    
    // debug
    //for (int i=0; i<World.TotalMatrixSize; i++)
    //{
    //    printf("JH = %lf\n", real(World.Hamiltonian(i,i)));
    //}
    // debug
    ofstream outputHamiltonian("Hamiltonian.txt");
    World.Hamiltonian.print(outputHamiltonian);
    outputHamiltonian.close();
    
    printf("calculating Ef=%lf...\n", Ef);
    World.CalculateGR(Ef);
    cx_mat T = World.ThermalAverageTransmission(Temperature, Ef);
    T.print("T = ");
    Complex P, Q;
    P = T(3,0)*T(1,0) + T(3,1)*T(1,0) + T(3,2)*T(1,0) + T(1,3)*T(3,0);
    Q = T(3,0)*T(1,2) + T(3,1)*T(1,2) + T(3,2)*T(1,2) + T(3,2)*T(1,3);
    Miu2 = real(P/(P+Q)- Q/(P+Q))*Miu1;
    P = T(1,0)*T(3,1) + T(1,0)*T(3,0) + T(1,2)*T(3,0) + T(1,3)*T(3,0);
    Q = T(1,0)*T(3,2) + T(1,2)*T(3,1) + T(1,2)*T(3,2) + T(1,3)*T(3,2);
    Miu4 = real(P/(P+Q) - Q/(P+Q))*Miu1;
    printf("%lf\t%lf\n", Miu2, Miu4);
    World.ListOfOpenBoundaries[0].ChemicalPotential = Miu1+Ef;
    World.ListOfOpenBoundaries[1].ChemicalPotential = Miu2+Ef;
    World.ListOfOpenBoundaries[2].ChemicalPotential = -Miu1+Ef;
    World.ListOfOpenBoundaries[3].ChemicalPotential = Miu4+Ef;
    //T = World.ThermalAverageTransmission(0.0, Ef);
    
    CalculateCurrentDistribution(Ef, Miu1, Miu2, Miu4, World);
   
    double Estart = -1.0;
    double Eend = 1.0;
    int Nstep = 100;
    double Estep = fabs((Eend-Estart)/(double)Nstep);
    double I, Sx, Sy, Sz;
    World.CalculateTerminalSpinCurrent(World.ListOfOpenBoundaries[0], Ef, I, Sx, Sy, Sz);
    totalI = I;
    printf("energy = %lf   I=% 8le  Sx=% 8le   Sy=% 8le   Sz=% 8le\n", Ef, I, Sx, Sy, Sz);
    World.CalculateTerminalSpinCurrent(World.ListOfOpenBoundaries[1], Ef, I, Sx, Sy, Sz);
    totalSz = Sz;
    printf("energy = %lf   I=% 8le  Sx=% 8le   Sy=% 8le   Sz=% 8le\n", Ef, I, Sx, Sy, Sz);
    World.CalculateTerminalSpinCurrent(World.ListOfOpenBoundaries[2], Ef, I, Sx, Sy, Sz);
    printf("energy = %lf   I=% 8le  Sx=% 8le   Sy=% 8le   Sz=% 8le\n", Ef, I, Sx, Sy, Sz);
    World.CalculateTerminalSpinCurrent(World.ListOfOpenBoundaries[3], Ef, I, Sx, Sy, Sz);
    printf("energy = %lf   I=% 8le  Sx=% 8le   Sy=% 8le   Sz=% 8le\n", Ef, I, Sx, Sy, Sz);
    double up, down;
    World.InjectionDOS(World.ListOfOpenBoundaries[0], up, down);
    UpDos = up;
    DownDos = down;
     // World.PlotCurrentMap(31,31,31,31);
     // World.PlotSpinXCurrentMap(30,30,30,30);
     // World.PlotSpinYCurrentMap(30,30,30,30);
      World.PlotSpinZCurrentMap(31,31,21,21);
      World.OutputSpinCurrentMapProFit("ProFitSzCurrentMap.txt",5,25,5,25);
      World.OutputSpinTextureProFit("ProFitSpinMap.txt");
      //World.OutputSpinCurrentMapProFit("ProFitSpinCurrentMap.txt");
    /*FILE* fp;
    fp = fopen("TerminalCurrent.txt","w");
    for (double Ef=Estart; Ef<Eend; Ef+=Estep)
    {
        World.CalculateGR(Ef);
        World.CalculateTerminalSpinCurrent(World.ListOfOpenBoundaries[1], I, Sx, Sy, Sz);
        printf("energy = %lf   I=% 8le  Sx=% 8le   Sy=% 8le   Sz=% 8le\n", Ef, I, Sx, Sy, Sz);
        fprintf(fp, "%lf\t% 8le\t% 8le\t% 8le\t% 8le\n", Ef, I, Sx, Sy, Sz);
    }
    fclose(fp);*/
}
/////////
/*void CalculateInjectionBandStructure(double t)
{
    ElectronSystem World("input.txt","openBoundaries.txt", "BoundaryVirtualShift.txt", t, 1.000001);
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
}*/



void CalculateInjectionBandStructure(double t)
{
    ElectronSystem World("input.txt","openBoundaries.txt","BoundaryVirtualShift.txt", t,15.001);
    cx_mat F00 = World.ListOfOpenBoundaries[0].F00;
    cx_mat F01 = World.ListOfOpenBoundaries[0].F01;
    FILE* fp;
    FILE* fp2;
    FILE* fp3;
    fp = fopen("BandStructure.txt","w");
    fp2 = fopen("EigVec.txt","w");
    fp3 = fopen("SzExpectationValue.txt","w");
    cout<<endl<<"size of F00  = "<<F01.n_rows<<"   by "<<F01.n_cols<<endl;
    cx_mat SZ(100,100);
    cx_mat psi, psihc;
    cx_mat SzExpectationValue;
    SZ.eye();
    int j3 = 0;
    
    for(int i3 = 0;i3 <50;i3++)
    {
        j3 = (2*i3+1);
        //cout <<"\n j3 ="<<j3<<"\n";
        SZ(j3,j3) = -1;
    }
    
    
    //cout<<SZ;
    
    for (double ka = 0.0; ka <= PIPI; ka+= PIPI/100.0)
    {
        cx_mat Hk = F00 + F01*exp(Complex(0.0, 1.0)*ka) + trans(F01)*exp(-Complex(0.0, 1.0)*ka);
        //vec eigval = eig_sym(Hk);
        vec eigval;
        cx_mat eigvec;
        eig_sym(eigval,eigvec,Hk);
        fprintf(fp, "% lf\t", ka);
        psi = eigvec.col(50);
        //psihc = trans(psi);
        
        
        
        SzExpectationValue = trans(psi)*SZ*psi;
        fprintf(fp3, "% lf\t  % lf\t  % lf\t % lf\n", ka,abs(SzExpectationValue(0,0)),real(SzExpectationValue(0,0)),imag(SzExpectationValue(0,0)));
        //fprintf(fp3, "\n");
        
        
        if(ka == 0.0)
        {
            //cout<<"\n eigvec no of col"<<eigvec.n_cols<<"\n";
            //cout<<"\n"<<eigvec(0,0)<<"\n";
            //cout<<"\n"<<abs(eigvec(0,0))<<"\n";
            //fprintf(fp2, "% lf\n", ka);
            
            for (int i2=0;i2<eigvec.n_rows;i2++)
            {
                for(int j2=0;j2<eigvec.n_cols;j2++)
                {
                    fprintf(fp2, "% le\t", abs(eigvec(i2,j2)));
                }
                fprintf(fp2, "\n");
            }
            
            
            
            //SzExpectationValue = trans(psi)*SZ*psi;
            //cout<<"\n SzExpectationValue ="<<SzExpectationValue<<"\n";
            //cout<<"\n psi row col "<<psihc.n_rows<<"  "<<psihc.n_cols<<"\n";
            
            
        }
        
        for (int i=0; i<eigval.n_rows; i++)
        {
            fprintf(fp, "% le\t", eigval(i));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    fclose(fp2);
    fclose(fp3);
}




/////////////
double TopologicalCharge(SpinSystem &input)
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
    std::vector<MagneticNode>::iterator i;
    vec Sb(3);
    vec Sc(3);
    vec Sd(3);
    vec Se(3);
    vec Sa(3);
    Sb.zeros();
    Sc.zeros();
    Sd.zeros();
    Se.zeros();
    Complex rho1, SigmaA1, EXP1;
    Complex rho2, SigmaA2, EXP2;
    Complex TQ = 0.0;
    for (i=input.NodeList.begin(); i!=input.NodeList.end(); i++)
    {
        for (int j=0; j<i->ListOfExchangeNeighbours.size(); j++)
        {
            vec temp(3);
            temp = i->ListOfOurwardsVectors[j];
            if (temp(0) > 0.5)
                Sb = input.NodeList[i->ListOfExchangeNeighbours[j]].Spin;
            if (temp(0) < -0.5)
                Sd = input.NodeList[i->ListOfExchangeNeighbours[j]].Spin;
            if (temp(1) > 0.5)
                Sc = input.NodeList[i->ListOfExchangeNeighbours[j]].Spin;
            if (temp(1) < -0.5)
                Se = input.NodeList[i->ListOfExchangeNeighbours[j]].Spin;
        }
        Sa = i->Spin;
        rho1 = sqrt(2.0*(1+dot(Sa, Sb))*
                         (1+dot(Sb, Sc))*
                         (1+dot(Sc, Sa)));
        EXP1 = (1.0+dot(Sa, Sb)+dot(Sb, Sc)+dot(Sc, Sa)
                     +Complex(0,1.)*dot(Sa, cross(Sb, Sc)))/rho1;
        SigmaA1 = 2.0*imag(log(EXP1));
                
        rho2 = sqrt(2.0*(1+dot(Sa, Sd))*
                         (1+dot(Sd, Se))*
                         (1+dot(Se, Sa)));
        EXP2 = (1.0+dot(Sa, Sd)+dot(Sd, Se)+dot(Se, Sa)
                     +Complex(0,1.)*dot(Sa, cross(Sd, Se)))/rho2;
        SigmaA2 = 2.0*imag(log(EXP2));
        TQ += real(1.0/4.0/PIPI*(SigmaA1+SigmaA2));
    }
    return real(TQ);
}



