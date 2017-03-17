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
#include <string>

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


void ElectronSystem::CreateSystemVariableAndConstants(void)
{
    this->NumberOfSiteAlongX = 5;
    this->NumberOfSiteAlongY = 5;
    this->delta = 15e-10;
    this->gamma = 5e-10*1;
    this->hv = 3.2955e-10*0.015;
    printf("\n Neighboring site distance = 15A,(width)numsite::AlongX= %d, AlongY= %d \n",this->NumberOfSiteAlongX,this->NumberOfSiteAlongY);
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
        //iss >> x >> y >> z >> sx >> sy >> sz;
        //sx = 0;
        //sy = 0;
        //sz = 1;
        newSpin << sx << sy << sz;
        //newSpin.print();
        newLocation << x << y << z;
        //newLocation.print();
        JH = -(this->JbyT*this->hv)/(this->delta);
        ElectronSite temp(newSpin, newLocation, JH, t);
        temp.SiteIndex = countSite;
        countTotalSize += temp.OnSiteBlock.n_cols;
        temp.StartingRowNumber = countTotalSize - temp.OnSiteBlock.n_rows;
        temp.SiteBlockSize = temp.OnSiteBlock.n_rows;
        countSite++;
        this->ListOfSites.push_back(temp);
    } while (!input.eof());
    this->NumSite = countSite;
    this->TotalMatrixSize = countTotalSize;
    cout<<"JbtT="<<this->JbyT<<"\n";
    input.close();
    
    
   
    
}
/////////////////////////////////////////
/*void ElectronSystem::CreateNeighbourList(double CutoffRange, double t)
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
}*/


void ElectronSystem::CreateNeighbourList(double CutoffRange, double t)
{
    cx_mat T;
    vec r(3);
    //the following variables occur in electronsite, electron system, openboundary
    int width = this->NumberOfSiteAlongX; //number of site along y
    int length = this->NumberOfSiteAlongY; //number of site along x
    double deltaX = this->delta; //interatomic distance along x direction
    double deltaY = this->delta; //intteratomic distance along y direction
    double hv = this->hv; //hbar*VF in electron volt
    double gamma = this->gamma;
    
    cx_mat Pauli_X(2,2);
    cx_mat Pauli_Y(2,2);
    cx_mat Pauli_Z(2,2);
    Pauli_X.zeros();
    Pauli_Y.zeros();
    Pauli_Z.zeros();
    Pauli_X(0,1) = 1.0;
    Pauli_X(1,0) = 1.0;
    Pauli_Y(0,1) = Complex(0.0, -1.0);
    Pauli_Y(1,0) = Complex(0.0, 1.0);
    Pauli_Z(0.0) = 1.0;
    Pauli_Z(1,1) = -1.0;
    
    cx_mat tx(2,2);
    cx_mat txhc(2,2);
    cx_mat ty(2,2);
    cx_mat tyhc(2,2);
    
    tx.zeros();
    txhc.zeros();
    ty.zeros();
    tyhc.zeros();
    
    this->txhc = hv*((Complex(0.0, -1.0)/(2*deltaX))*(Pauli_Y)-(gamma/(deltaX*deltaX))*(Pauli_Z));  //points to the right in the x direction
    txhc = this->txhc;
    tx = trans(txhc);
    this->tyhc = trans(hv*((Complex(0.0, 1.0)/(2*deltaX))*(Pauli_X)-(gamma/(deltaX*deltaX))*(Pauli_Z))); //points up in the y direction
    tyhc = this->tyhc;
    ty = trans(tyhc);

    
    
    for (int i=0; i<NumSite; i++)
    {
        //cout<<"\n for site "<<i<<"neighbour =";
        if(((i+1)%width)!=0)
        {
                ListOfSites[i].ListOfTightBindingNeighbors.push_back(i+1);
                ListOfSites[i].ListOfOutwardsCouplingBlocks.push_back(txhc);
               // cout<<i+1<<"\t";
        }
        if((i%width)!=0)
        {
                ListOfSites[i].ListOfTightBindingNeighbors.push_back(i-1);
                ListOfSites[i].ListOfOutwardsCouplingBlocks.push_back(tx);
                //cout<<i-1<<"\t";
        }

        if((i+width)<=(NumSite-1))
        {
                ListOfSites[i].ListOfTightBindingNeighbors.push_back(i+width);
                ListOfSites[i].ListOfOutwardsCouplingBlocks.push_back(tyhc);
                //cout<<i+width<<"\t";
        }
        
        if((i-width)>=0)
        {
                ListOfSites[i].ListOfTightBindingNeighbors.push_back(i-width);
                ListOfSites[i].ListOfOutwardsCouplingBlocks.push_back(ty);
                //cout<<i-width<<"\t";
        }
        
         
    }
    
    
    printf("create neighborlist done \n"); //tonmoy
    
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



void ElectronSystem::GenerateH01x()       //generates coupling matrix between boundary sites for calculating dispersion in the x direction
{
    //works as long as brillouin zone is rectangular
    this->H01x.set_size(this->TotalMatrixSize, this->TotalMatrixSize);
    this->H01x.zeros();
    int nx = this->NumberOfSiteAlongX;
    int ny = this->NumberOfSiteAlongY;
    int n1 = 0; //index of first member of first row
    int n2= nx-1; //index of last member of first row
    
    for(int i=0;i<ny;i++)
    {
        //cout<<n1<<"  "<<n2<<"\n";
        this->H01x.submat(n2*2,n1*2,n2*2+1,n1*2+1) = this->txhc;  //2 means number of on site basis
        n1 = n1 + nx;
        n2 = n2 + nx;       
        
    }  
 
}



void ElectronSystem::GenerateH01y()       //generates coupling matrix between boundary sites for calculating dispersion in the x direction
{
    //works as long as brillouin zone is rectangular
    this->H01y.set_size(this->TotalMatrixSize, this->TotalMatrixSize);
    this->H01y.zeros();
    int nx = this->NumberOfSiteAlongX;
    int ny = this->NumberOfSiteAlongY;
    int n1 = 0; //index of first member of first column
    int n2= nx*(ny-1) ; //index of last member of first column
    
    for(int i=0;i<nx;i++)
    {
        //cout<<n1<<"  "<<n2<<"\n";
        this->H01y.submat(n2*2,n1*2,n2*2+1,n1*2+1) = this->tyhc;  //2 means number of on site basis
        n1 = n1 + 1;
        n2 = n2 + 1;       
        
    }  
 
}

void ElectronSystem::CalculateDispersion()
{
    
        double pi = 3.14156;
        double kx, ky;
        int berryCurvatureforBand ;
                
        FILE* fp1;
        fp1 = fopen("./OUTPUT/Dispersion.txt","a");
        FILE* fp2;
        fp2 = fopen("./OUTPUT/Psi.txt","w");
        FILE* fp3;
        fp3 = fopen("./OUTPUT/BerryCurvature.txt","w");
        FILE* fp4;
        fp4 = fopen("./OUTPUT/SzExpectation.txt","a");
        cout<<"\n calculating Dispersion \n";
        cx_mat HK;
        
        double nanoribbon = 0.0;
        double psiSquare;
        double BerryCurvature2;
        cx_mat dummy1;
        cx_mat dummy2;
        cx_mat Pauli_Z(2,2);
        Pauli_Z.zeros();
        Pauli_Z(0.0) = 1.0;
        Pauli_Z(1,1) = -1.0;
        cx_mat SZ(50,50);
        SZ.zeros();
        
        for(int i=0;i<50;i=i+2)
        {
            SZ.submat(i,i,i+1,i+1) = Pauli_Z;
        }
        
        //cout<<SZ;
        
        //fprintf(fp7,"Summing Berry Curvature for all Bands for a particular K point and Printing that value \n");
        fprintf(fp1, "% le\t",this->JbyT);
        for(ky = 0;ky<= 6.3;ky = ky + pi/10)  //pi, put 3.2 instead of pi, originall pi/30
        {
            
            for(kx = 0;kx<=6.3 ;kx = kx + pi/10) //pi, put 3.2 instead of pi, originally pi/30
            {
                
                //cout<<"kx ="<<kx<<"ky ="<<ky<<"\n";
                //if(ky!=kx)
                    //continue;
                
                HK = this->Hamiltonian + (this->H01x)*exp(Complex(0.0, 1.0)*kx) + trans(this->H01x)*exp(-Complex(0.0, 1.0)*kx)
                        + (this->H01y)*exp(Complex(0.0, 1.0)*ky)  + trans(this->H01y)*exp(-Complex(0.0, 1.0)*ky);
                
                eig_sym(this->eigval,this->eigvec,HK);
                
            /*    dummy2 = trans(this->eigvec.col(24))*SZ*this->eigvec.col(24);
               // cout<<real(dummy2(0,0));
                fprintf(fp4, "%le\t%le",this->JbyT,real(dummy2(0,0)));
                fprintf(fp4, "\n");   */
                
                //each wavefunction armadillo calculates is indeed normalized
                //ends checking whether wavefunction is normalized
                
                
                fprintf(fp3, "% le\t% le\t",kx,ky);
                for(int n3= 0;n3<=-1;n3++)   //under normal condition, n3=0 to 49
                {
                    fprintf(fp3, "% le\t",BerryCurvature(kx,ky,n3));
                }
                fprintf(fp3, "\n");
                
                  
                
                //fprintf(fp1, "% le\t% le\t",kx,ky);
                //for(int k3=0;k3<this->eigval.n_rows;k3++)
                for(int k3=this->NumberOfSiteAlongX*this->NumberOfSiteAlongY;k3<this->NumberOfSiteAlongX*this->NumberOfSiteAlongY+8;k3++)
                {
                      fprintf(fp1, "% le\t",this->eigval(k3,0));
                      //fprintf(fp1, "\n");
                }
               
                
                
            }
            
            
        }
        

        
        

   /*
        for (int i22=0;i22<this->eigvec.n_rows;i22++)
        {
            for(int j22=0;j22<this->eigvec.n_cols;j22++)
            {
                fprintf(fp2, "% le\t", abs(this->eigvec(i22,j22))*abs(this->eigvec(i22,j22)));
                
            }
                fprintf(fp2, "\n");
    } 
    */    
        fprintf(fp1, "\n");
        fclose(fp1);
        fclose(fp2);
        fclose(fp3);
        fclose(fp4);
}



double ElectronSystem::BerryCurvature(double kx,double ky, int n)
{
    //calculating Berry Curvature for nth band at brillouin zone point kx, ky
    cx_mat dhdkx = (this->H01x)*exp(Complex(0.0, 1.0)*kx) - trans(this->H01x)*exp(-Complex(0.0, 1.0)*kx);
    cx_mat dhdky = (this->H01y)*exp(Complex(0.0, 1.0)*ky)  - trans(this->H01y)*exp(-Complex(0.0, 1.0)*ky);
    
    cx_mat  bc1, bc2, bc3;   //bc contains berry curvature
    double bc = 0.0;
    
    for(int i=0;i<this->TotalMatrixSize;i++)
    {
        if(i==n)
            continue;
        
       
        bc1 = trans((this->eigvec).col(n))*dhdkx*(this->eigvec).col(i);
        bc2 = trans((this->eigvec).col(i))*dhdky*(this->eigvec).col(n);
        bc3 = (bc1*bc2)/((this->eigval(n,0)-this->eigval(i,0))*(this->eigval(n,0)-this->eigval(i,0)));
        
        
        bc = bc + 2*imag(bc3(0,0)); //notice the factor 2 here, since a- conj(a) = i*2Imag(a))

    }
    
    bc = bc * -(this->NumberOfSiteAlongX*this->delta*this->NumberOfSiteAlongY*this->delta);   //that factor = (11*15A0)^2
    
    return bc;
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
        temp.VirtualBoundaryShift.resize(3);
        ListOfOpenBoundaries.push_back(temp);
        
        count++;
    } while (!input.eof());
    input.close();
    printf("%d open boundaries detected.\n", count);
    
}
////////////
void ElectronSystem::ReadInOpenBoundaryVirtialShift(const char* filename)
{
    FILE *fp;
    fp = fopen(filename, "r");
    vec VirtualShiftTemp(3);
    double ax, ay, az;
    for (int i=0; i<this->ListOfOpenBoundaries.size(); i++)
    {
        fscanf(fp, "%le%le%le", &ax, &ay, &az);
        VirtualShiftTemp(0) = ax;
        VirtualShiftTemp(1) = ay;
        VirtualShiftTemp(2) = az;
        this->ListOfOpenBoundaries[i].VirtualBoundaryShift = VirtualShiftTemp;
    }
}
///////////
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
    printf("inverting...\n");
    //ofstream test("testdebug.txt");
    //OpenHamiltonian.print(test);
    GR = inv(OpenHamiltonian);
}
//////////////
void ElectronSystem::RenewGR(double energy) //no need
{
    //this->GenerateHamiltonian();
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
void ElectronSystem::UpdateHamiltonian(SpinSystem &SpinTexture, double JH) //no need
{
    cx_mat Pauli_X(2,2);
    cx_mat Pauli_Y(2,2);
    cx_mat Pauli_Z(2,2);
    vec spinTemp(3);
    Pauli_X.zeros();
    Pauli_Y.zeros();
    Pauli_Z.zeros();
    Pauli_X(0,1) = 1.0;
    Pauli_X(1,0) = 1.0;
    Pauli_Y(0,1) = Complex(0.0, -1.0);
    Pauli_Y(1,0) = Complex(0.0, 1.0);
    Pauli_Z(0.0) = 1.0;
    Pauli_Z(1,1) = -1.0;
    
    
    for (int i=0; i<SpinTexture.ListOfTorqueSiteIndecies.size(); i++)
    {
        int MagneticIndex = SpinTexture.ListOfTorqueSiteIndecies[i];
        int ElectronicIndex = SpinTexture.NodeList[MagneticIndex].ElectronSiteIndex;
        spinTemp = SpinTexture.NodeList[MagneticIndex].Spin;;
        this->ListOfSites[ElectronicIndex].Spin = spinTemp;
        this->ListOfSites[ElectronicIndex].OnSiteBlock =  JH*spinTemp(0)*Pauli_X
                        +JH*spinTemp(1)*Pauli_Y
                        +JH*spinTemp(2)*Pauli_Z;
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
void ElectronSystem::GenerateBoundaryHamiltonians(double CouplingT, double CouplingCutoff)
{
    for (int i=0; i<this->ListOfOpenBoundaries.size(); i++)
    {
        this->ListOfOpenBoundaries[i].ConstructF00F01(CouplingT, CouplingCutoff);
    }
    printf("%lu Hamiltonians of boundaries created.\n", this->ListOfOpenBoundaries.size());
}
///////////
ElectronSystem::ElectronSystem(const char* inputFilename, const char* boundaryFilename,
                               const char* inputBoundaryShiftFilename, double JbyT,double CouplingCutoff)
{
    double t = -1.00;
    this->JbyT = JbyT; 
    this->CreateSystemVariableAndConstants();
    this->ReadInGeometry(inputFilename);
    this->CreateNeighbourList(CouplingCutoff, t);
    //this->PrintNeighbourList();
    printf("\n total matrix size = %d\n", this->TotalMatrixSize);  
    this->GenerateHamiltonian();
    //this->ReadInOpenBoundaries(boundaryFilename);
    //this->ReadInOpenBoundaryVirtialShift(inputBoundaryShiftFilename);
    //this->GenerateBoundaryHamiltonians(t, CouplingCutoff);
    //this->PrintReadInSkyrmionTexture();
    //this->CalculateBoundState();
    //this->CalculateVelocity2();
}


void ElectronSystem::CalculateVelocity(void)
{
    FILE* fp;
    fp = fopen("velocityField.txt","w");
    cout<<"\n Calculating velocity Field for five lowest Hole & Electron state\n";
    int totalNumberofSite = this->NumSite;
    int iPlus1, iMinus1, jPlus1, jMinus1; //neighbour index for each site
    int numberOfSiteAlongX = this->NumberOfSiteAlongX;
    int numberOfSiteAlongY = this->NumberOfSiteAlongY;
    double delta = this->delta; //interatomic distance along x direction
    
    
    for(int k4=2495;k4<2506;k4++) // k4  =  k4 th eigenvalue
    {
        mat a(2,2);
        a.eye();
        cx_mat iPlus = (1/delta)*Complex(0.0, 1.0)*a;
        cx_mat iMinus = (1/delta)*Complex(0.0, -1.0)*a;
        cx_mat vx(this->TotalMatrixSize,this->TotalMatrixSize),vy(this->TotalMatrixSize,this->TotalMatrixSize);
        vx.zeros(), vy.zeros();
        cx_mat avgX, avgY;
        
        
        for(int k10=0;k10<totalNumberofSite;k10++) //k10 = index of each site
        {
            //calculating velocity x component
            if((k10%numberOfSiteAlongX)==0)
            {
                iPlus1 = k10+1;
                iMinus1 = (k10+numberOfSiteAlongX)-1;       
            }
            else if(((k10+1)%numberOfSiteAlongX)==0)
            {
                iPlus1 = k10-(numberOfSiteAlongX-1);
                iMinus1 = k10-1;  
            }
            else
            {
                iPlus1 = k10+1;
                iMinus1 = k10-1;  
            }
            
            vx.submat(k10*2,iPlus1*2,k10*2+1,iPlus1*2+1) = iPlus;
            vx.submat(k10*2,iMinus1*2,k10*2+1,iMinus1*2+1) = iMinus;
            
         
            
            //calculating velocity y component
            if((k10-numberOfSiteAlongX)<0)
            {
                jPlus1 = k10+numberOfSiteAlongX;
                jMinus1 = k10 + numberOfSiteAlongX*(numberOfSiteAlongX-1);       
            }
            else if((k10+numberOfSiteAlongX)>(totalNumberofSite-1))
            {
                jPlus1 = k10 - numberOfSiteAlongX*(numberOfSiteAlongX-1);
                jMinus1 = k10-numberOfSiteAlongX;  
            }
            else
            {
                jPlus1 = k10+numberOfSiteAlongX;
                jMinus1 = k10-numberOfSiteAlongX;  
            }
            
            
            vy.submat(k10*2,jPlus1*2,k10*2+1,jPlus1*2+1) = iPlus;
            vy.submat(k10*2,jMinus1*2,k10*2+1,jMinus1*2+1) = iMinus;
 
        }
        
        
        cx_mat Ix, Iy;
        Ix = vx*this->eigvec.col(k4);
        Iy = vy*this->eigvec.col(k4);
        
        
        
        for(int k10=0;k10<totalNumberofSite;k10++)
        {
            avgX = trans(this->eigvec.submat(k10*2,k4,k10*2+1,k4))*Ix.submat(k10*2,0,k10*2+1,0);
            avgY = trans(this->eigvec.submat(k10*2,k4,k10*2+1,k4))*Iy.submat(k10*2,0,k10*2+1,0);
            fprintf(fp, "% le\t% le\t% le\t% le\t% le\n",ListOfSites[k10].Location(0),ListOfSites[k10].Location(1),ListOfSites[k10].Location(2),real(avgX(0,0)),real(avgY(0,0)));
            //cout<<"\n wave function "<<real(eigvec1.submat())<<"\n";
            //cout<<"\n real part = "<<imag(avgX(0,0))<<"\n";
        }
        //cout<<"\nreal velocity vx = "<<real(Ix(0,0))<<"\t"<<"vy = "<<real(Iy(0,0))<<"\n";
        //cout<<"\nimaginary velocity vx = "<<imag(Ix(0,0))<<"\t"<<"vy = "<<imag(Iy(0,0))<<"\n";
    }
    
    fclose(fp);
    
}

void ElectronSystem::CalculateVelocity2(void)
{
    FILE* fp;
    fp = fopen("velocityField.txt","w");
    int width = this->NumberOfSiteAlongX; //number of site along y
    int length = this->NumberOfSiteAlongY; //number of site along x
    double deltaX = this->delta; //interatomic distance along x direction
    double deltaY = this->delta; //intteratomic distance along y direction
    double hv = 3.2955e-10*0.5; //hbar*VF in electron volt
    double gamma = 5e-10;
    int totalNumberofSite = this->NumSite;
    int numberOfSiteAlongX = this->NumberOfSiteAlongX;
    int numberOfSiteAlongY = this->NumberOfSiteAlongY;
    int iPlus1, iMinus1, jPlus1, jMinus1; //neighbour index for each site
    
    cx_mat Pauli_X(2,2);
    cx_mat Pauli_Y(2,2);
    cx_mat Pauli_Z(2,2);
    Pauli_X.zeros();
    Pauli_Y.zeros();
    Pauli_Z.zeros();
    Pauli_X(0,1) = 1.0;
    Pauli_X(1,0) = 1.0;
    Pauli_Y(0,1) = Complex(0.0, -1.0);
    Pauli_Y(1,0) = Complex(0.0, 1.0);
    Pauli_Z(0.0) = 1.0;
    Pauli_Z(1,1) = -1.0;
    
    cx_mat tx(2,2);
    cx_mat txhc(2,2);
    cx_mat ty(2,2);
    cx_mat tyhc(2,2);
    
    tx.zeros();
    txhc.zeros();
    ty.zeros();
    tyhc.zeros();
    
    txhc = hv*((Complex(0.0, -1.0)/(2*deltaX))*(Pauli_Y)-(gamma/(deltaX*deltaX))*(Pauli_Z));
    tx = trans(txhc);
    tyhc = hv*((Complex(0.0, 1.0)/(2*deltaX))*(Pauli_X)-(gamma/(deltaX*deltaX))*(Pauli_Z));
    ty = trans(tyhc);
    double avgX, avgY;
    cx_mat dummy1, dummy2, dummy3, dummy4;
        
    for(int k4=2495;k4<2506;k4++) // k4  =  k4 th eigenvalue
    {    
        for(int k10=0;k10<totalNumberofSite;k10++) //k10 = index of each site
        {
            //calculating velocity x component
            if((k10%numberOfSiteAlongX)==0)
            {
                iPlus1 = k10+1;
                iMinus1 = (k10+numberOfSiteAlongX)-1;       
            }
            else if(((k10+1)%numberOfSiteAlongX)==0)
            {
                iPlus1 = k10-(numberOfSiteAlongX-1);
                iMinus1 = k10-1;  
            }
            else
            {
                iPlus1 = k10+1;
                iMinus1 = k10-1;  
            }
            

           // avgX = imag(trans(this->eigvec.submat(k10*2,k4,k10*2+1,k4))*txhc*(this->eigvec.submat(iPlus1*2,k4,iPlus1*2+1,k4))) + 
            //       imag(trans(this->eigvec.submat(iMinus1*2,k4,iMinus1*2+1,k4))*txhc*(this->eigvec.submat(k10*2,k4,k10*2+1,k4)));
                   
         
            dummy1=  trans(this->eigvec.submat(k10*2,k4,k10*2+1,k4))*Pauli_Z*txhc*(this->eigvec.submat(iPlus1*2,k4,iPlus1*2+1,k4));
            dummy2 = trans(this->eigvec.submat(iMinus1*2,k4,iMinus1*2+1,k4))*Pauli_Z*txhc*(this->eigvec.submat(k10*2,k4,k10*2+1,k4));
            avgX = imag(dummy1(0,0)) + imag(dummy2(0,0));
            
            
            //For calculating spin current, insert Pauli_Z before txhc, tyhc
            
            //calculating velocity y component
            if((k10-numberOfSiteAlongX)<0)
            {
                jPlus1 = k10+numberOfSiteAlongX;
                jMinus1 = k10 + numberOfSiteAlongX*(numberOfSiteAlongX-1);       
            }
            else if((k10+numberOfSiteAlongX)>(totalNumberofSite-1))
            {
                jPlus1 = k10 - numberOfSiteAlongX*(numberOfSiteAlongX-1);
                jMinus1 = k10-numberOfSiteAlongX;  
            }
            else
            {
                jPlus1 = k10+numberOfSiteAlongX;
                jMinus1 = k10-numberOfSiteAlongX;  
            }
            
            //avgY = imag(trans(this->eigvec.submat(k10*2,k4,k10*2+1,k4))*tyhc*(this->eigvec.submat(jPlus1*2,k4,jPlus1*2+1,k4))) + 
            //       imag(trans(this->eigvec.submat(jMinus1*2,k4,jMinus1*2+1,k4))*tyhc*(this->eigvec.submat(k10*2,k4,k10*2+1,k4)));
            
            
            dummy3 = trans(this->eigvec.submat(k10*2,k4,k10*2+1,k4))*Pauli_Z*tyhc*(this->eigvec.submat(jPlus1*2,k4,jPlus1*2+1,k4));
            dummy4 = trans(this->eigvec.submat(jMinus1*2,k4,jMinus1*2+1,k4))*Pauli_Z*tyhc*(this->eigvec.submat(k10*2,k4,k10*2+1,k4));
            avgY = imag(dummy3(0,0)) + imag(dummy4(0,0));
            
            fprintf(fp, "% le\t% le\t% le\t% le\t% le\n",ListOfSites[k10].Location(0),ListOfSites[k10].Location(1),ListOfSites[k10].Location(2),(avgX),(avgY));

 
        }
    
    }        
            
        fclose(fp);    
}

//////////////////////

void ElectronSystem::CalculateTransmissionForEnergyRange(double Energy)
{
    FILE* fp;
    fp = fopen("Xmission.txt","w");
    for(double e= -Energy;e<= Energy;e=e+0.001)
    {
            this->CalculateGR(e);
            cx_mat tc;
            tc = this->Transmission(e);
            cout<<"\n transmission matrix \n"<<tc<<"\n";
            cout<<"\n transmission = "<<real(tc(0,1))<<"\n";
            fprintf(fp, "% lf\t% lf\n", e,real(tc(0,1)));
    }
    
    
    fclose(fp);
}



void ElectronSystem::CalculateBoundState(void)
{
        FILE* fp1;
        fp1 = fopen("EigenVauleSkx.txt","w");
        FILE* fp2;
        fp2 = fopen("Psi.txt","w");
        cout<<"\n calculating Bound state \n";
        eig_sym(this->eigval,this->eigvec,this->Hamiltonian);
        
        for(int k3=0;k3<this->eigval.n_rows;k3++)
        {
            fprintf(fp1, "% le\n",this->eigval(k3,0));
        }
   
        for (int i22=0;i22<this->eigvec.n_rows;i22++)
        {
            for(int j22=0;j22<this->eigvec.n_cols;j22++)
            {
                fprintf(fp2, "% le\t", abs(this->eigvec(i22,j22))*abs(this->eigvec(i22,j22)));
                
            }
                fprintf(fp2, "\n");
    } 
        
        
        fclose(fp1);
        fclose(fp2);
}


void ElectronSystem::PrintReadInSkyrmionTexture(void)
{
    FILE* fp;
    fp = fopen("SKXtexture.txt","w"); 
    int totalNumberofSite = this->NumSite;
    for(int k14=0;k14<totalNumberofSite;k14++)
    {
        fprintf(fp, "% le\t% le\t% le\t% le\t% le\t% le\n",ListOfSites[k14].Location(0),ListOfSites[k14].Location(1),ListOfSites[k14].Location(2),ListOfSites[k14].Spin(0),ListOfSites[k14].Spin(1),ListOfSites[k14].Spin(2));
    }
    fclose(fp);
}







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


cx_mat ElectronSystem::OnSiteSpectralFunctionA(int SiteI, int SiteJ)
{
    cx_mat A;
    A.set_size(this->ListOfSites[SiteI].SiteBlockSize, this->ListOfSites[SiteJ].SiteBlockSize);
    A.zeros();
    for (int i=0; i<this->ListOfOpenBoundaries.size(); i++)
    {
        A = A + OnSiteSpectralFunctionFromBoundary(SiteI, SiteJ, ListOfOpenBoundaries[i]);
                
    }
    return A;
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
//////////////
vec ElectronSystem::OnSiteSpin(int SiteIndex, double Ef)
{
    cx_mat x,y,z;
    x << 0.0 << 1.0 << endr
      << 1.0 << 0.0 << endr;
    y << 0.0 << Complex(0.,-1.) << endr
      << Complex(0.,1.) << 0.0 << endr;
    z << 1.000000 << 0.0 << endr
      << 0.0 << -1.000000 << endr; 
    vec SpinVector(3);
    cx_mat Gn = this->OnSiteCorelationFunctionGn(SiteIndex, SiteIndex, Ef);
    SpinVector(0) = real(trace(x*Gn));
    SpinVector(1) = real(trace(y*Gn));
    SpinVector(2) = real(trace(z*Gn));
    return SpinVector;
    
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
void ElectronSystem::OutputElectronSpinMapProFit(const char* filename, double Ef)
{
    printf("Generating ProFit data for electron spin map...\n");
    FILE *fp;
    fp = fopen(filename, "w");
    double x, y, r, theta, sx, sy, sz;
    vec tempSpin(3);
    for (int i=0; i<this->NumSite; i++)
    {
        x = this->ListOfSites[i].Location(0);
        y = this->ListOfSites[i].Location(1);
        tempSpin = this->OnSiteSpin(i, Ef);
        sx = tempSpin(0);
        sy = tempSpin(1);
        sz = tempSpin(2);
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
/////////

void ElectronSystem::CalculateCurrentDistribution(double Ef)
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
    cx_mat an; // Spectral function from Terminal n;
    an.set_size(2,2);
    for (int i=0; i<this->NumSite; i++)
    {
        vec Current(3);
        vec SpinXCurrent(3);
        vec SpinYCurrent(3);
        vec SpinZCurrent(3);
        Current.zeros();
        SpinXCurrent.zeros();
        SpinYCurrent.zeros();
        SpinZCurrent.zeros();
        if (this->ListOfSites[i].IsBoundary == true)
        {
            ListOfSites[i].Current = Current;
            ListOfSites[i].SpinXCurrent = SpinXCurrent;
            ListOfSites[i].SpinYCurrent = SpinYCurrent;
            ListOfSites[i].SpinZCurrent = SpinZCurrent;
            continue;
        }
         
        cx_mat I_local;
        I_local.set_size(2,2);
        I_local.zeros();
        double DeltaMiu;
        for (int j=0; j<ListOfSites[i].ListOfTightBindingNeighbors.size(); j++)
        {
            int SourceIndex = ListOfSites[i].ListOfTightBindingNeighbors[j];
            int I, J;
            I = ListOfSites[i].SiteIndex;
            J = ListOfSites[SourceIndex].SiteIndex;
            for (int p=0; p<this->ListOfOpenBoundaries.size(); p++)
            {
                DeltaMiu = ListOfOpenBoundaries[p].ChemicalPotential - Ef;
                an = OnSiteSpectralFunctionFromBoundary(I, J, ListOfOpenBoundaries[p]);
                I_local += -Complex(0.0, 1.0)*(DeltaMiu*(an-trans(an)));
            }
            Current = Current + real(trace(I_local))*(ListOfSites[i].Location-ListOfSites[SourceIndex].Location);
            SpinXCurrent = SpinXCurrent + real(trace(SigmaX*I_local))*(ListOfSites[i].Location-ListOfSites[SourceIndex].Location);
            SpinYCurrent = SpinYCurrent + real(trace(SigmaY*I_local))*(ListOfSites[i].Location-ListOfSites[SourceIndex].Location);
            SpinZCurrent = SpinZCurrent + real(trace(SigmaZ*I_local))*(ListOfSites[i].Location-ListOfSites[SourceIndex].Location);
        }
        ListOfSites[i].Current = Current;
        ListOfSites[i].SpinXCurrent = SpinXCurrent;
        ListOfSites[i].SpinYCurrent = SpinYCurrent;
        ListOfSites[i].SpinZCurrent = SpinZCurrent; 
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
    for (int i=0; i<this->NumSite; i++)
    {
        ListOfSites[i].Current = ListOfSites[i].Current/MaxCurrent;
        ListOfSites[i].SpinXCurrent = ListOfSites[i].SpinXCurrent/MaxSxCurrent;
        ListOfSites[i].SpinYCurrent = ListOfSites[i].SpinYCurrent/MaxSyCurrent;
        ListOfSites[i].SpinZCurrent = ListOfSites[i].SpinZCurrent/MaxSzCurrent;  
    }
    FILE* fp;
    fp = fopen("LocalMiu.txt","w");
    Complex delta;
    Complex Norm;
    delta = 0.0;
    Norm = 0.0;
    double DeltaMiu;
    for (int i=0; i<this->NumSite; i++)
    {
        for (int p=0; p<this->ListOfOpenBoundaries.size(); p++)
        {
            DeltaMiu = this->ListOfOpenBoundaries[p].ChemicalPotential-Ef;
            delta += trace(OnSiteSpectralFunctionFromBoundary(i,i,ListOfOpenBoundaries[p])*DeltaMiu);
            Norm += trace(OnSiteSpectralFunctionFromBoundary(i,i,ListOfOpenBoundaries[p]));
        }
        fprintf(fp, "%lf\t%lf\t%le\n", ListOfSites[i].Location(0), ListOfSites[i].Location(1), real(delta/Norm));
                
    }
    fclose(fp);
}    


