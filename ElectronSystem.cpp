#include "ElectronSystem.hpp"
#include "SpinSystem.hpp"
#include "SpinSystem.hpp"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

//#include <boost/algorithm/string.hpp>

using namespace arma;
//using namespace boost::algorithm;
typedef complex<double> Complex;

//The following are some string operation functions.
// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}
// trim from both ends
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}
// String operatirons end.

double ElectronSystem::dFdE(double energy, double Ef, double temperature)
{
        //Unit is in eV;
        double kT = 8.6173324e-5*temperature;
        double x = (energy-Ef)/kT;
       	return (-1.0/kT)*1.0/(exp(-x)+2.0+exp(x));
}

void ElectronSystem::ReadInGeometry(const char* filename)
{
    std::ifstream input(filename);
    std::string line;
    //FILE *fp;
    //fp = fopen("filename","r");
    double x, y, z, sx, sy, sz, JH, t, phi;
    int countSite = 0;
    int countTotalSize = 0;
    vec newSpin(3);
    vec newLocation(3);
    getline(input, line);
    do 
    {
        //fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf", &x, &y, &z, &sx, &sy, &sz, &JH, &t, &phi);
        getline(input, line);
        if (line.empty() || input.eof())
            break;
        std::istringstream iss(line);
        iss >> x >> y >> z >> sx >> sy >> sz >> JH >> t >> phi;
        newSpin << sx << sy << sz;
        //newSpin.print();
        newLocation << x << y << z;
        //newLocation.print();
        ElectronSite temp(newSpin, newLocation, JH, phi);
        temp.SiteIndex = countSite;
        countTotalSize += temp.OnSiteBlock.n_cols;
        temp.StartingRowNumber = countTotalSize - temp.OnSiteBlock.n_rows;
        temp.SiteBlockSize = temp.OnSiteBlock.n_rows;
        countSite++;
        this->ListOfSites.push_back(temp);
    } while (!input.eof());
    this->NumSite = countSite;
    this->TotalMatrixSize = countTotalSize;
    input.close();
}
/////////////////////////////////////////
void ElectronSystem::CreateNeighbourList(double CutoffRange, double t)
{
    cx_mat T;
    vec r(3);
    for (int i=0; i<NumSite; i++)
    {
        for (int j=0; j<NumSite; j++)
        {
            r = (ListOfSites[i].Location)-(ListOfSites[j].Location);
            double distance_square = r(0)*r(0)+r(1)*r(1)+r(2)*r(2);
            if (distance_square < CutoffRange*CutoffRange && i!=j)
            {
                ListOfSites[i].ListOfTightBindingNeighbors.push_back(j);
                T.eye(ListOfSites[i].OnSiteBlock.n_rows, ListOfSites[i].OnSiteBlock.n_rows);
                ListOfSites[i].ListOfOutwardsCouplingBlocks.push_back(T*t);
            }
               // if (distance_square < StrayCutoff*StrayCutoff)
                 //   NodeList[i].ListOfStrayFieldNeighbours.push_back(j);
        }
    }
}
/////////////////////////////////////////
void ElectronSystem::GenerateHamiltonian()
{
    this->Hamiltonian.set_size(this->TotalMatrixSize, this->TotalMatrixSize);
    this->Hamiltonian.zeros();
    for (int i=0; i<this->NumSite; i++)
    {
        // First, on-site blocks.
        int OnSiteSize = this->ListOfSites[i].OnSiteBlock.n_rows;
        int I = this->ListOfSites[i].SiteIndex;
        this->Hamiltonian.submat(I*OnSiteSize,I*OnSiteSize,(I+1)*OnSiteSize-1,(I+1)*OnSiteSize-1) = this->ListOfSites[i].OnSiteBlock; 
        // Then coupling blocks.
        for (int j=0; j<ListOfSites[i].ListOfTightBindingNeighbors.size(); j++)
        {
            int J = ListOfSites[i].ListOfTightBindingNeighbors[j];
            this->Hamiltonian.submat(I*OnSiteSize,J*OnSiteSize,(I+1)*OnSiteSize-1,(J+1)*OnSiteSize-1) 
                    = this->ListOfSites[i].ListOfOutwardsCouplingBlocks[j]; 
        }        
    } 
}
/////////////////////////

void ElectronSystem::ReadInOpenBoundaries(const char* filename)
{
    std::ifstream input(filename);
    std::string line;
    int SiteIndex;
    int BoundaryIndex = 0;
    int count = 0;
    do 
    {
        OpenBoundary temp;
        getline(input, line);
        if (line.empty())
            break;
        trim(line);
        temp.BoundaryIndex = BoundaryIndex++;
        std::istringstream iss(line);
        int TotalMatrixSizeCount = 0;
        do
        {
            if (iss.eof())
                break;
            iss >> SiteIndex;
            this->ListOfSites[SiteIndex].IsBoundary = true;
            temp.ListOfBoundarySites.push_back(this->ListOfSites[SiteIndex]);
            TotalMatrixSizeCount += this->ListOfSites[SiteIndex].OnSiteBlock.n_rows;
        }while (!iss.eof());
        temp.TotalBoundaryMatrixSize = TotalMatrixSizeCount;
        ListOfOpenBoundaries.push_back(temp);
        count++;
    } while (!input.eof());
    input.close();
    printf("%d open boundaries detected.\n", count);
    
}
////////////
void ElectronSystem::PrintBoundaryList(void)
{
    for (int i=0; i<ListOfOpenBoundaries.size(); i++)
    {
        printf("Boundary %d; Size=%lu  ; ", i, ListOfOpenBoundaries[i].ListOfBoundarySites.size());
        for (int j=0; j<ListOfOpenBoundaries[i].ListOfBoundarySites.size(); j++)
        {
            printf("%d\t", ListOfOpenBoundaries[i].ListOfBoundarySites[j].SiteIndex);
        }
        printf("\n");
    }
}
/////////////
void ElectronSystem::PrintNeighbourList(void)
{
    for (int i=0; i<NumSite; i++)
    {
        printf("%d    ", i);
        for (int j=0; j<ListOfSites[i].ListOfTightBindingNeighbors.size(); j++)
        {
            printf("%d   ", ListOfSites[i].ListOfTightBindingNeighbors[j]);
        }
        printf("\n");
    }
}
///////////////
void ElectronSystem::CalculateGR(double energy)
{
    cx_mat OpenHamiltonian;
    OpenHamiltonian = this->Hamiltonian;
    //printf("debug.   %d\n", TotalMatrixSize);
    //printf("%d", this->TotalMatrixSize);
    cx_mat I;
    //printf("debug2\n");
    I.eye(TotalMatrixSize, TotalMatrixSize);
    //printf("I generated. %d\n", ListOfOpenBoundaries.size());
    
    for (int i=0; i<this->ListOfOpenBoundaries.size(); i++)
    {
        //printf("Trying to get self\n");
        this->ListOfOpenBoundaries[i].GetSelfEnergy(energy);
        //printf("self generated\n");
        this->ListOfOpenBoundaries[i].AddSelfEnergy(OpenHamiltonian);
        this->ListOfOpenBoundaries[i].GetGamma();
        //printf("added\n");
    }
    OpenHamiltonian = energy*I - OpenHamiltonian;
    //printf("inverting...\n");
    //ofstream test("testdebug.txt");
    //OpenHamiltonian.print(test);
    GR = inv(OpenHamiltonian);
}
//////////////
void ElectronSystem::RenewGR(double energy)
{
    this->GenerateHamiltonian();
    cx_mat OpenHamiltonian;
    OpenHamiltonian = this->Hamiltonian;
    //printf("debug.   %d\n", TotalMatrixSize);
    //printf("%d", this->TotalMatrixSize);
    cx_mat I;
    //printf("debug2\n");
    I.eye(TotalMatrixSize, TotalMatrixSize);
    //printf("I generated. %d\n", ListOfOpenBoundaries.size());
    
    for (int i=0; i<this->ListOfOpenBoundaries.size(); i++)
    {
        //printf("Trying to get self\n");
        //this->ListOfOpenBoundaries[i].GetSelfEnergy(energy);
        //printf("self generated\n");
        this->ListOfOpenBoundaries[i].AddSelfEnergy(OpenHamiltonian);
        this->ListOfOpenBoundaries[i].GetGamma();
        //printf("added\n");
    }
    OpenHamiltonian = energy*I - OpenHamiltonian;
    //printf("inverting...\n");
    //ofstream test("testdebug.txt");
    //OpenHamiltonian.print(test);
    GR = inv(OpenHamiltonian);
}
//////////////
void ElectronSystem::UpdateHamiltonian(SpinSystem &SpinTexture)
{
    for (int i=0; i<SpinTexture.ListOfTorqueSiteIndecies.size(); i++)
    {
        int MagneticIndex = SpinTexture.ListOfTorqueSiteIndecies[i];
        int ElectronicIndex = SpinTexture.NodeList[MagneticIndex].ElectronSiteIndex;
        this->ListOfSites[ElectronicIndex].Spin = SpinTexture.NodeList[MagneticIndex].Spin;
    }
    this->GenerateHamiltonian();
}
/////////////
cx_mat ElectronSystem::ObtainGR_AB(OpenBoundary A, OpenBoundary B)
{
    cx_mat GR_AB;
    GR_AB.zeros();
    GR_AB.set_size(A.TotalBoundaryMatrixSize, B.TotalBoundaryMatrixSize);
    for (int i=0; i<A.ListOfBoundarySites.size(); i++)
    {
        for (int j=0; j<B.ListOfBoundarySites.size(); j++)
        {
            int StartITotalMatrix = A.ListOfBoundarySites[i].StartingRowNumber;
            int StartJTotalMatrix = B.ListOfBoundarySites[j].StartingRowNumber;
            int StartIBoundaryMatrix = A.StartingRowInBoundaryMatrix[i];
            int StartJBoundaryMatrix = B.StartingRowInBoundaryMatrix[j];
            int NumRowI = A.ListOfBoundarySites[i].SiteBlockSize;
            int NumRowJ = B.ListOfBoundarySites[j].SiteBlockSize;
            GR_AB.submat(StartIBoundaryMatrix,StartJBoundaryMatrix,StartIBoundaryMatrix+NumRowI-1, StartJBoundaryMatrix+NumRowJ-1)
            = this->GR.submat(StartITotalMatrix, StartJTotalMatrix, StartITotalMatrix+NumRowI-1, StartJTotalMatrix+NumRowJ-1);
        }
    }
    return GR_AB;
}
////////////
cx_mat ElectronSystem::Transmission(double energy)
{
    //ATTENTION: Before this function, one need to call CalculateGR first.
    //this->CalculateGR(energy);
    //ofstream fileGR("GR.txt");
    //this->GR.print(fileGR);
    //fileGR.close();
    for (int i=0; i<ListOfOpenBoundaries.size(); i++)
    {
        ListOfOpenBoundaries[i].GetGamma();
        //ListOfOpenBoundaries[i].Gamma.print("gamma=");
        //ListOfOpenBoundaries[i].F00.print("F00");
        //ListOfOpenBoundaries[i].F01.print("F01");
    }
    cx_mat T;
    T.set_size(ListOfOpenBoundaries.size(), ListOfOpenBoundaries.size());
    T.zeros();
    for (int i=0; i<T.n_rows; i++)
    {
        for (int j=0; j<T.n_cols; j++)
        {
            cx_mat GR_AB;
            GR_AB = ObtainGR_AB(ListOfOpenBoundaries[i], ListOfOpenBoundaries[j]);
            //printf("i,j = %d, %d   ", i, j);
            //GR_AB.print("GR_AB=");
            T(i,j) = trace(ListOfOpenBoundaries[i].Gamma*GR_AB*ListOfOpenBoundaries[j].Gamma*trans(GR_AB));
        }
    }   
    return T;    
}
///////////
void ElectronSystem::GenerateBoundaryHamiltonians(double CouplingT)
{
    for (int i=0; i<this->ListOfOpenBoundaries.size(); i++)
    {
        this->ListOfOpenBoundaries[i].ConstructF00F01(CouplingT);
    }
    printf("%lu Hamiltonians of boundaries created.\n", this->ListOfOpenBoundaries.size());
}
///////////
ElectronSystem::ElectronSystem(const char* inputFilename, const char* boundaryFilename, double t)
{
    this->ReadInGeometry(inputFilename);
    this->CreateNeighbourList(1.00001, t);
    //this->PrintNeighbourList();
    printf("%d\n", this->TotalMatrixSize);  
    this->GenerateHamiltonian();
    this->ReadInOpenBoundaries(boundaryFilename);
    this->GenerateBoundaryHamiltonians(t);
}
//////////////////////
cx_mat ElectronSystem::ThermalAverageTransmission(double Temperature, double Ef)
{
    if (fabs(Temperature)<1.0e-7)
    {
        // this is for the case of zero temperature;
        return this->Transmission(Ef);
    }
    cx_mat AveragedT;
    AveragedT.set_size(this->ListOfOpenBoundaries.size(), this->ListOfOpenBoundaries.size());
    AveragedT.zeros();
    double kT = 8.6173324e-5*Temperature;
    double Estep = 24.0*kT/20.0;
    int count = 0;
    for (double Energy = Ef-12.0*kT; Energy <= Ef +12.0*kT; Energy += Estep)
    {
        AveragedT += -dFdE(Energy, Ef, Temperature)*Estep*this->Transmission(Energy);
        count ++;
    }
    return AveragedT;
}
///////
cx_mat ElectronSystem::GnBlockFromBoundary(OpenBoundary &Source, double Energy, double Temperature)
{
    cx_mat GnBlock;    
    GnBlock.set_size(Source.TotalBoundaryMatrixSize, Source.TotalBoundaryMatrixSize);
    for (int I=0; I<Source.ListOfBoundarySites.size(); I++)
    {
        for (int J=0; J<Source.ListOfBoundarySites.size(); J++)
        {
            int i=Source.StartingRowInBoundaryMatrix[I];
            int j=Source.StartingRowInBoundaryMatrix[J];
            cx_mat gn(2,2);
            gn.zeros();
            for (int k=0; k<this->ListOfOpenBoundaries.size(); k++)
            {
                cx_mat a;
                a = this->OnSiteSpectralFunctionFromBoundary(Source.ListOfBoundarySites[I].SiteIndex, 
                                                             Source.ListOfBoundarySites[J].SiteIndex, 
                                                             this->ListOfOpenBoundaries[k]);
                gn = gn + a*Fermi_eV(Energy, Temperature, this->ListOfOpenBoundaries[k].ChemicalPotential);    
            }    
            GnBlock.submat(i,j,i+1,j+1) = gn;
        }
    } 
    return GnBlock;
}
///////
cx_mat ElectronSystem::dGnBlockFromBoundary(OpenBoundary& Source, double Ef) 
{
    cx_mat dGnBlock;    
    dGnBlock.set_size(Source.TotalBoundaryMatrixSize, Source.TotalBoundaryMatrixSize);
    for (int I=0; I<Source.ListOfBoundarySites.size(); I++)
    {
        for (int J=0; J<Source.ListOfBoundarySites.size(); J++)
        {
            int i=Source.StartingRowInBoundaryMatrix[I];
            int j=Source.StartingRowInBoundaryMatrix[J];
            cx_mat dgn(2,2);
            dgn.zeros();
            for (int k=0; k<this->ListOfOpenBoundaries.size(); k++)
            {
                cx_mat a;
                a = this->OnSiteSpectralFunctionFromBoundary(Source.ListOfBoundarySites[I].SiteIndex, 
                                                             Source.ListOfBoundarySites[J].SiteIndex, 
                                                             this->ListOfOpenBoundaries[k]);
                dgn = dgn + a*(this->ListOfOpenBoundaries[k].ChemicalPotential-Ef);    
            }    
            dGnBlock.submat(i,j,i+1,j+1) = dgn;
        }
    } 
    return dGnBlock;
}

///////
void ElectronSystem::CalculateTerminalSpinCurrent(OpenBoundary& Source, double Ef,
                                                    double &I, double &Sx, double &Sy, double &Sz)
{
    cx_mat Gr;
    cx_mat dGn;
    double PI = 3.1415926535897932384626;
    double hbar = 6.58211928e-16; //eV*s
    dGn = dGnBlockFromBoundary(Source, Ef);
    Gr.set_size(Source.TotalBoundaryMatrixSize, Source.TotalBoundaryMatrixSize);
    Gr = this->ObtainGR_AB(Source, Source);
    cx_mat I_operator;
    I_operator = Complex(0.0, 1.0)/2.0/PI/hbar*(dGn*trans(Source.SelfEnergy)-Source.SelfEnergy*dGn
                       +Gr*Source.Gamma*(Source.ChemicalPotential-Ef)-Source.Gamma*trans(Gr)*(Source.ChemicalPotential-Ef));
    cx_mat SigmaX, SigmaY, SigmaZ;
    cx_mat x,y,z;
    x << 0.0 << 1.0 << endr
      << 1.0 << 0.0 << endr;
    y << 0.0 << Complex(0.,-1.) << endr
      << Complex(0.,1.) << 0.0 << endr;
    z << 1.000000 << 0.0 << endr
      << 0.0 << -1.000000 << endr; 
    SigmaX.set_size(Source.TotalBoundaryMatrixSize, Source.TotalBoundaryMatrixSize);
    SigmaY.set_size(Source.TotalBoundaryMatrixSize, Source.TotalBoundaryMatrixSize);
    SigmaZ.set_size(Source.TotalBoundaryMatrixSize, Source.TotalBoundaryMatrixSize);
    for (int i=0; i<Source.TotalBoundaryMatrixSize; i+=2)
    {
        SigmaX.submat(i,i,i+1,i+1) = x;
        SigmaY.submat(i,i,i+1,i+1) = y;
        SigmaZ.submat(i,i,i+1,i+1) = z;
    }
    I = real(trace(I_operator));//*1.602176565e-19;
    Sx = real(trace(SigmaX*I_operator));//*hbar/2.0;
    Sy = real(trace(SigmaY*I_operator));//*hbar/2.0;
    Sz = real(trace(SigmaZ*I_operator));//*hbar/2.0;
}

cx_mat ElectronSystem::OnSiteSpectralFunctionFromBoundary(int I, int J, OpenBoundary& Source)
{
    int i = this->ListOfSites[I].StartingRowNumber;
    int j = this->ListOfSites[J].StartingRowNumber;
    cx_mat Aij(2,2);
    Aij.zeros();
    for (int P=0; P<Source.ListOfBoundarySites.size(); P++)
    {
        for (int Q=0; Q<Source.ListOfBoundarySites.size(); Q++)
        {
            cx_mat gr;
            cx_mat ga;
            cx_mat gamma;
            int p = Source.ListOfBoundarySites[P].StartingRowNumber;
            int q = Source.ListOfBoundarySites[Q].StartingRowNumber;
            gr = this->GR.submat(i,p,i+1,p+1);
            ga = trans(this->GR.submat(j,q,j+1,q+1));
            gamma = Source.Gamma.submat(Source.StartingRowInBoundaryMatrix[P], Source.StartingRowInBoundaryMatrix[Q], 
                                        Source.StartingRowInBoundaryMatrix[P]+1, Source.StartingRowInBoundaryMatrix[Q]+1);
            Aij = Aij + gr*gamma*ga; 
        }
    }
    return Aij;
}
///
cx_mat ElectronSystem::OnSiteCorelationFunctionGn(int SiteI, int SiteJ, double Ef)
{
    cx_mat GN;
    GN.set_size(this->ListOfSites[SiteI].SiteBlockSize, this->ListOfSites[SiteJ].SiteBlockSize);
    GN.zeros();
    for (int i=0; i<this->ListOfOpenBoundaries.size(); i++)
    {
        GN = GN + OnSiteSpectralFunctionFromBoundary(SiteI, SiteJ, ListOfOpenBoundaries[i])*
                (ListOfOpenBoundaries[i].ChemicalPotential-Ef);
    }
    return GN;
}
///
void ElectronSystem::InjectionDOS(OpenBoundary& Source, double &Up, double &Down)
{
    Complex UpTemp, DownTemp;
    cx_mat Spectral;
    UpTemp = 0.0;
    DownTemp = 0.0;
    for (int i; i<Source.ListOfBoundarySites.size(); i++)
    {
        UpTemp += -imag(Source.SurfaceGreen(2*i, 2*i));
        DownTemp += -imag(Source.SurfaceGreen(2*i+1, 2*i+1));
    }
    Up = real(UpTemp);
    Down = real(DownTemp);
}

void ElectronSystem::PlotCurrentMap(int SizeX, int SizeY, int PlotX, int PlotY)
{
    double MaxLength = 0.0;
    // To renormalize the arrows within the plot. find the maximum arrow length 
    for (int i=0; i<this->NumSite; i++)
    {
        if (fabs(ListOfSites[i].Location(0)-SizeX/2)<fabs(PlotX/2) &&
            fabs(ListOfSites[i].Location(1)-SizeY/2)<fabs(PlotY/2) )
        {
            double length = sqrt(dot(ListOfSites[i].Current,ListOfSites[i].Current));
            if (length > MaxLength)
                MaxLength = length;
        }        
    }
   float* xp = new float[this->NumSite];
   float* yp = new float[this->NumSite];
   float* xv = new float[this->NumSite];
   float* yv = new float[this->NumSite];
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
  //Graph.scrmod ("black");
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
    sprintf(buffer, "hello.\n");
    for (int i=0; i<this->NumSite; i++)
    {
        xp[i] = (float)(this->ListOfSites[i].Location(0));
        yp[i] = (float)(this->ListOfSites[i].Location(1));
        if (fabs(ListOfSites[i].Location(0)-SizeX/2)<fabs(PlotX/2) &&
            fabs(ListOfSites[i].Location(1)-SizeY/2)<fabs(PlotY/2) )
        {
            xv[i] = (float)(this->ListOfSites[i].Current(0)/MaxLength);
            yv[i] = (float)(this->ListOfSites[i].Current(1)/MaxLength);
        }
        else
        {
            xv[i] = 0.0;
            yv[i] = 0.0;
        }
    }    
    Graph.titlin (buffer, 4);
    Graph.title ();
    Graph.vecfld (xv, yv, xp, yp, this->NumSite, 1901);
    //Graph.crvmat (pointer, PlotX, PlotY, 15, 15);
    
    
    //endgrf();
    Graph.sendbf();
  delete [] xp;
  delete [] yp;
  delete [] xv;
  delete [] yv;
  Graph.disfin();
}
///////////
void ElectronSystem::PlotSpinXCurrentMap(int SizeX, int SizeY, int PlotX, int PlotY)
{
    double MaxLength = 0.0;
    // To renormalize the arrows within the plot. find the maximum arrow length 
    for (int i=0; i<this->NumSite; i++)
    {
        if (fabs(ListOfSites[i].Location(0)-SizeX/2)<fabs(PlotX/2) &&
            fabs(ListOfSites[i].Location(1)-SizeY/2)<fabs(PlotY/2) )
        {
            double length = sqrt(dot(ListOfSites[i].SpinXCurrent,ListOfSites[i].SpinXCurrent));
            if (length > MaxLength)
                MaxLength = length;
        }        
    }
   float* xp = new float[this->NumSite];
   float* yp = new float[this->NumSite];
   float* xv = new float[this->NumSite];
   float* yv = new float[this->NumSite];
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
  //Graph.scrmod ("black");
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
    sprintf(buffer, "hello.\n");
    for (int i=0; i<this->NumSite; i++)
    {
        xp[i] = (float)(this->ListOfSites[i].Location(0));
        yp[i] = (float)(this->ListOfSites[i].Location(1));
        if (fabs(ListOfSites[i].Location(0)-SizeX/2)<fabs(PlotX/2) &&
            fabs(ListOfSites[i].Location(1)-SizeY/2)<fabs(PlotY/2) )
        {
            xv[i] = (float)(this->ListOfSites[i].SpinXCurrent(0)/MaxLength);
            yv[i] = (float)(this->ListOfSites[i].SpinXCurrent(1)/MaxLength);
        }
        else
        {
            xv[i] = 0.0;
            yv[i] = 0.0;
        }
    }    
    Graph.titlin (buffer, 4);
    Graph.title ();
    Graph.vecfld (xv, yv, xp, yp, this->NumSite, 1901);
    //Graph.crvmat (pointer, PlotX, PlotY, 15, 15);
    
    
    //endgrf();
    Graph.sendbf();
  delete [] xp;
  delete [] yp;
  delete [] xv;
  delete [] yv;
  Graph.disfin();
}
/////

void ElectronSystem::PlotSpinYCurrentMap(int SizeX, int SizeY, int PlotX, int PlotY)
{
    double MaxLength = 0.0;
    // To renormalize the arrows within the plot. find the maximum arrow length 
    for (int i=0; i<this->NumSite; i++)
    {
        if (fabs(ListOfSites[i].Location(0)-SizeX/2)<fabs(PlotX/2) &&
            fabs(ListOfSites[i].Location(1)-SizeY/2)<fabs(PlotY/2) )
        {
            double length = sqrt(dot(ListOfSites[i].SpinYCurrent,ListOfSites[i].SpinYCurrent));
            if (length > MaxLength)
                MaxLength = length;
        }        
    }
   float* xp = new float[this->NumSite];
   float* yp = new float[this->NumSite];
   float* xv = new float[this->NumSite];
   float* yv = new float[this->NumSite];
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
  //Graph.scrmod ("black");
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
    sprintf(buffer, "hello.\n");
    for (int i=0; i<this->NumSite; i++)
    {
        xp[i] = (float)(this->ListOfSites[i].Location(0));
        yp[i] = (float)(this->ListOfSites[i].Location(1));
        if (fabs(ListOfSites[i].Location(0)-SizeX/2)<fabs(PlotX/2) &&
            fabs(ListOfSites[i].Location(1)-SizeY/2)<fabs(PlotY/2) )
        {
            xv[i] = (float)(this->ListOfSites[i].SpinYCurrent(0)/MaxLength);
            yv[i] = (float)(this->ListOfSites[i].SpinYCurrent(1)/MaxLength);
        }
        else
        {
            xv[i] = 0.0;
            yv[i] = 0.0;
        }
    }    
    Graph.titlin (buffer, 4);
    Graph.title ();
    Graph.vecfld (xv, yv, xp, yp, this->NumSite, 1901);
    //Graph.crvmat (pointer, PlotX, PlotY, 15, 15);
    
    
    //endgrf();
    Graph.sendbf();
  delete [] xp;
  delete [] yp;
  delete [] xv;
  delete [] yv;
  Graph.disfin();
}
//////
void ElectronSystem::PlotSpinZCurrentMap(int SizeX, int SizeY, int PlotX, int PlotY)
{
    double MaxLength = 0.0;
    // To renormalize the arrows within the plot. find the maximum arrow length 
    for (int i=0; i<this->NumSite; i++)
    {
        if (fabs(ListOfSites[i].Location(0)-SizeX/2)<fabs(PlotX/2) &&
            fabs(ListOfSites[i].Location(1)-SizeY/2)<fabs(PlotY/2) )
        {
            double length = sqrt(dot(ListOfSites[i].SpinZCurrent,ListOfSites[i].SpinZCurrent));
            if (length > MaxLength)
                MaxLength = length;
        }        
    }
   float* xp = new float[this->NumSite];
   float* yp = new float[this->NumSite];
   float* xv = new float[this->NumSite];
   float* yv = new float[this->NumSite];
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
  //Graph.scrmod ("black");
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
    sprintf(buffer, "hello.\n");
    for (int i=0; i<this->NumSite; i++)
    {
        xp[i] = (float)(this->ListOfSites[i].Location(0));
        yp[i] = (float)(this->ListOfSites[i].Location(1));
        if (fabs(ListOfSites[i].Location(0)-SizeX/2)<fabs(PlotX/2) &&
            fabs(ListOfSites[i].Location(1)-SizeY/2)<fabs(PlotY/2) )
        {
            xv[i] = (float)(this->ListOfSites[i].SpinZCurrent(0)/MaxLength);
            yv[i] = (float)(this->ListOfSites[i].SpinZCurrent(1)/MaxLength);
        }
        else
        {
            xv[i] = 0.0;
            yv[i] = 0.0;
        }
    }    
    Graph.titlin (buffer, 4);
    Graph.title ();
    Graph.vecfld (xv, yv, xp, yp, this->NumSite, 1901);
    //Graph.crvmat (pointer, PlotX, PlotY, 15, 15);
    
    
    //endgrf();
    Graph.sendbf();
  delete [] xp;
  delete [] yp;
  delete [] xv;
  delete [] yv;
  Graph.disfin();
}
/////
void ElectronSystem::OutputSpinTextureProFit(const char* filename)
{
    printf("Generating ProFit data for spin plot...\n");
    FILE *fp;
    fp = fopen(filename, "w");
    double x, y, r, theta, sx, sy, sz;
    for (int i=0; i<this->NumSite; i++)
    {
        x = this->ListOfSites[i].Location(0);
        y = this->ListOfSites[i].Location(1);
        sx = this->ListOfSites[i].Spin(0);
        sy = this->ListOfSites[i].Spin(1);
        sz = this->ListOfSites[i].Spin(2);
        r = sqrt(sx*sx+sy*sy);
        theta = atan2(sy, sx);
        fprintf(fp, "% le\t% le\t% le\t% le\t% le\n", x, y, r, theta, sz);
    }
    fclose(fp);
}
/////
void ElectronSystem::OutputSpinCurrentMapProFit(const char* filename, int Xstart,int Xend,int Ystart,int Yend)
{
    printf("Generating ProFit data for current plot...\n");
    FILE *fp;
    fp = fopen(filename, "w");
    double maxR = 0.0;
    for (int i=0; i<this->NumSite; i++)
    {
        double x = this->ListOfSites[i].Location(0);
        double y = this->ListOfSites[i].Location(1);
        if ( x >= (double)Xstart 
             && x <= (double)Xend
             && y >= (double)Ystart
             && y <= (double)Yend)
        {
            double jx = ListOfSites[i].SpinZCurrent(0);
            double jy = ListOfSites[i].SpinZCurrent(1);
            double r = sqrt(jx*jx+jy*jy);
            if (r>maxR)
                maxR = r;
        }
    }
    double x, y, r, theta, Jx, Jy, Jz;
    for (int i=0; i<this->NumSite; i++)
    {
        x = this->ListOfSites[i].Location(0);
        y = this->ListOfSites[i].Location(1);
        Jx = this->ListOfSites[i].SpinZCurrent(0);
        Jy = this->ListOfSites[i].SpinZCurrent(1);
        Jz = this->ListOfSites[i].SpinZCurrent(2);
        r = sqrt(Jx*Jx+Jy*Jy)/maxR;
        theta = atan2(Jy, Jx);
        fprintf(fp, "% le\t% le\t% le\t% le\t% le\n", x, y, r, theta, Jz);
    }
    fclose(fp);
}
