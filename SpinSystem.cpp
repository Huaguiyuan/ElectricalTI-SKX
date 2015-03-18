#include "SpinSystem.hpp"
#include "ElectronSystem.hpp"
//#include "NeighborList.hpp"
#include <armadillo>
#include <math.h>
#define MiuZero 4.0e-7*3.1415926535897932384626 // In SI unit [N/A/A]

using namespace arma;
using namespace std;

/*void SpinSystem::CalculateStrayFieldTensor(void)
{
    // This routine assumes that dx dy and dz are all the same.
    for (int i=0; i<SizeX*SizeY*SizeZ; i++)
    {
        int NumOfNeighbours = NodeList[i].ListOfStrayFieldNeighbours.size();
        for (int j=0; j<NumOfNeighbours; j++)
        {
            vec PQ = NodeList[NodeList[i].ListOfStrayFieldNeighbours[j]].Location
                     -NodeList[i].Location;
            mat temp(3,3);
            temp = Tensor(PQ);
            NodeList[i].ListOfTensors.push_back(temp);
        }        
    }
}
///////////////////////////////////////////////////////////////////////
void SpinSystem::FormatTheSystem(vec S_input)
{
    for (int i=0; i<SizeX*SizeY*SizeZ; i++)
    {
        NodeList[i].Spin = S_input;
    }
}
///////////////////////////////////////////////////////////////////////
mat SpinSystem::Tensor(vec PQ)
{
    mat N(3,3);
    // This routine assumes dx=dy=dz=1
    double I = PQ(0);
    double J = PQ(1);
    double K = PQ(2);
    N.zeros();
    for (int i=0; i<=1; i++)
    {
        for (int j=0; j<=1; j++)
        {
            for (int k=0; k<=1; k++)
            {
                double rijkP = sqrt((I+i-0.5)*(I+i-0.5)+(J+j-0.5)*(J+j-0.5)+(K+k-0.5)*(K+k-0.5));
                N(0,0) += pow(-1.0, i+j+k)*atan(((K+k-0.5)*(J+j-0.5)/rijkP/(I+i-0.5)));
                N(1,1) += pow(-1.0, j+k+i)*atan(((I+i-0.5)*(K+k-0.5)/rijkP/(J+j-0.5)));
                N(2,2) += pow(-1.0, k+i+j)*atan(((J+j-0.5)*(I+i-0.5)/rijkP/(K+k-0.5)));
                N(0,1) += -pow(-1.0, i+j+k)*log((K+k-0.5)+rijkP);
                N(1,2) += -pow(-1.0, j+k+i)*log((I+i-0.5)+rijkP);
                N(2,0) += -pow(-1.0, k+i+j)*log((J+j-0.5)+rijkP);
               
            }
        }
    }
    N(1,0) = N(0,1);
    N(2,1) = N(1,2);
    N(0,2) = N(2,0);
    N = N/4.0/3.14159265;
    return N;
}
///////////////////////////////////////////////////////////////////////
*/
void SpinSystem::GenerateNeighborList(double ExchangeCutoff, double StrayCutoff, 
        bool PeriodicX, bool PeriodicY, bool PeriodicZ, vec ax, vec ay, vec az)
{
    
    for (int i=0; i<this->NumSite; i++)
    {
        for (int j=0; j<this->NumSite; j++)
        {
                vec r = (NodeList[i].Location)-(NodeList[j].Location);
                double distance_square = r(0)*r(0)+r(1)*r(1)+r(2)*r(2);
                if (distance_square < ExchangeCutoff*ExchangeCutoff && i!=j)
                {
                    NodeList[i].ListOfExchangeNeighbours.push_back(j);
                    NodeList[i].ListOfOurwardsVectors.push_back(NodeList[j].Location - NodeList[i].Location);
                }
               // if (distance_square < StrayCutoff*StrayCutoff)
                 //   NodeList[i].ListOfStrayFieldNeighbours.push_back(j);  
        }
    }
    if (PeriodicX == true)
    {
        for (int i=0; i<this->NumSite; i++)
        {
            for (int j=0; j<this->NumSite; j++)
            {
                    vec r = (NodeList[j].Location+ax)-(NodeList[i].Location);
                    double distance_square = r(0)*r(0)+r(1)*r(1)+r(2)*r(2);
                    if (distance_square < ExchangeCutoff*ExchangeCutoff)
                    {
                        NodeList[i].ListOfExchangeNeighbours.push_back(j);
                        NodeList[i].ListOfOurwardsVectors.push_back(NodeList[j].Location+ax - NodeList[i].Location);
                    }
                    r = (NodeList[j].Location-ax)-(NodeList[i].Location);
                    distance_square = r(0)*r(0)+r(1)*r(1)+r(2)*r(2);
                    if (distance_square < ExchangeCutoff*ExchangeCutoff)
                    {
                        NodeList[i].ListOfExchangeNeighbours.push_back(j);
                        NodeList[i].ListOfOurwardsVectors.push_back(NodeList[j].Location-ax - NodeList[i].Location);
                    }
            }
        }
    }
    if (PeriodicY == true)
    {
        for (int i=0; i<this->NumSite; i++)
        {
            for (int j=0; j<this->NumSite; j++)
            {
                    vec r = (NodeList[j].Location+ay)-(NodeList[i].Location);
                    double distance_square = r(0)*r(0)+r(1)*r(1)+r(2)*r(2);
                    if (distance_square < ExchangeCutoff*ExchangeCutoff)
                    {
                        NodeList[i].ListOfExchangeNeighbours.push_back(j);
                        NodeList[i].ListOfOurwardsVectors.push_back(NodeList[j].Location+ay - NodeList[i].Location);
                    }
                    r = (NodeList[j].Location-ay)-(NodeList[i].Location);
                    distance_square = r(0)*r(0)+r(1)*r(1)+r(2)*r(2);
                    if (distance_square < ExchangeCutoff*ExchangeCutoff)
                    {
                        NodeList[i].ListOfExchangeNeighbours.push_back(j);
                        NodeList[i].ListOfOurwardsVectors.push_back(NodeList[j].Location-ay - NodeList[i].Location);
                    }
            }
        }
    }
    if (PeriodicZ == true)
    {
        for (int i=0; i<this->NumSite; i++)
        {
            for (int j=0; j<this->NumSite; j++)
            {
                    vec r = (NodeList[j].Location+az)-(NodeList[i].Location);
                    double distance_square = r(0)*r(0)+r(1)*r(1)+r(2)*r(2);
                    if (distance_square < ExchangeCutoff*ExchangeCutoff)
                    {
                        NodeList[i].ListOfExchangeNeighbours.push_back(j);
                        NodeList[i].ListOfOurwardsVectors.push_back(NodeList[j].Location+az - NodeList[i].Location);
                    }
                    r = (NodeList[j].Location-az)-(NodeList[i].Location);
                    distance_square = r(0)*r(0)+r(1)*r(1)+r(2)*r(2);
                    if (distance_square < ExchangeCutoff*ExchangeCutoff)
                    {
                        NodeList[i].ListOfExchangeNeighbours.push_back(j);
                        NodeList[i].ListOfOurwardsVectors.push_back(NodeList[j].Location-az - NodeList[i].Location);
                    }
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////

SpinSystem::SpinSystem(const char* filename, double J_initial, double D_initial, double New_alpha)
{
    J = J_initial;
    D = D_initial;
    int seed = (int)time(0);
    RanGen = new CRandomMersenne(seed);
    std::ifstream input(filename);
    std::string line;
    double x, y, z, sx, sy, sz;
    int countSite = 0;
    vec newSpin(3);
    vec newLocation(3);
    getline(input, line);
    do 
    {
        getline(input, line);
        if (line.empty() || input.eof())
            break;
        std::istringstream iss(line);
        iss >> x >> y >> z >> sx >> sy >> sz;
        newSpin << sx << sy << sz;
        //newSpin.print();
        newLocation << x << y << z;
        //newLocation.print();
        MagneticNode temp(countSite, newLocation, newSpin);
        countSite++;
        temp.Pinned = false;
        temp.Temperature = 0.0;
        temp.alpha = New_alpha;
        this->NodeList.push_back(temp);
    } while (!input.eof());
    this->NumSite = countSite;
    input.close();
    //GenerateNeighborList(ExchangeCutoff, StrayCutoff, true, true, false);
    vec FieldInitializer(3);
    FieldInitializer.zeros();
    for (int i=0; i<this->NumSite; i++)
    {
        this->ExternalField.push_back(FieldInitializer);
        this->EffectiveField.push_back(FieldInitializer);
        this->RandomField.push_back(FieldInitializer);
        this->BackgroundField.push_back(FieldInitializer);
    }
}
///////////////////////////////////////////////////////////////////////

SpinSystem::~SpinSystem()
{
    delete RanGen;
}

/////////////////////////////////////////////////////////

void SpinSystem::CalculateEffectiveField(bool AddTorque)
{
    vec Heff(3);
    //printf("debug %lu\n", this->ExternalField.size());
    std::vector<MagneticNode>::iterator i;
    for (i=this->NodeList.begin(); i!=this->NodeList.end(); i++)
    {
        Heff.zeros();
        std::vector<int>::iterator j;
        int count=0;
        for (j=i->ListOfExchangeNeighbours.begin(); j!=i->ListOfExchangeNeighbours.end(); j++)
        {
            int J_index = *j;
            vec r = i->ListOfOurwardsVectors[count++];
            Heff = Heff + J*(NodeList[J_index].Spin) - D*cross(NodeList[J_index].Spin, r); 
        }
        Heff = Heff + this->ExternalField[i->Index]; 
        Heff = Heff + this->BackgroundField[i->Index];
        this->EffectiveField[i->Index] = Heff;
       // Heff.print();
    }
    if (AddTorque == true)
    {
        for (int i=0; i<this->ListOfTorqueSiteIndecies.size(); i++)
        {
            int Index = ListOfTorqueSiteIndecies[i];
            Heff = this->EffectiveField[Index];
            this ->EffectiveField[Index] = Heff + this->NodeList[Index].TorqueField;
        }
    }
}
/////////////
void SpinSystem::PrintSpinTextureToTextFile(const char* filename)
{
    FILE* fp;
    fp = fopen(filename, "w");
    for (int i=0; i<NumSite; i++)
    {
        fprintf(fp, "%d    %lf   %lf   %lf   %lf   %lf   %lf\n", NodeList[i].Index,
                                                            NodeList[i].Location(0), 
                                                            NodeList[i].Location(1),
                                                            NodeList[i].Location(2),
                                                            NodeList[i].Spin(0),
                                                            NodeList[i].Spin(1),
                                                            NodeList[i].Spin(2)
                );
    }
    fclose(fp);
}
//////
void SpinSystem::CalculateRandomField(double TimeStep)
{
    vec L(3);
    double TempRand;
    double Lambda;
    double alpha;
    std::vector<MagneticNode>::iterator i;
    for (i=NodeList.begin(); i!=NodeList.end(); i++)
    {
        alpha = i->alpha;
        Lambda= alpha/(1.+alpha*alpha);
        for (int j=0; j<3; j++)
        {
            double TempRand = RanGen->Random() - 0.5;
            L(j) = sqrt(24.0*Lambda*i->Temperature*TimeStep)*TempRand;
        }
        this->RandomField[i->Index] = L;
    }
}
///////
void SpinSystem::Evolve(double TimeStep, double Time, bool AddTorque)
{
    // This function evolves Plate AND ExternalField from time to time+TimeStep
    // Random force is added to the effective field at each time step.
    vec A(3);
    vec BL(3);
    vec A_tilde(3);
    vec BL_tilde(3);
    vec S(3);
    vec L(3);
    vec Heff(3);
    std::vector<vec> S_tilde;
    std::vector<vec> SpinUpdate;
    S_tilde.resize(this->NumSite);
    SpinUpdate.resize(this->NumSite);
    this->CalculateRandomField(TimeStep);
    this->UpdateExternalField(Time);
    this->CalculateEffectiveField(AddTorque);
    for (std::vector<MagneticNode>::iterator i=NodeList.begin(); i!=NodeList.end(); i++)
    {
        S = i->Spin;
        L = RandomField[i->Index];
        Heff = EffectiveField[i->Index];
        double alpha = i->alpha;
        A = -1.0/(1+alpha*alpha)*cross(S, Heff)-alpha/(1.0+alpha*alpha)*(S*(dot(S, Heff)) - Heff);
        BL = -1.0/(1+alpha*alpha)*cross(S, L)-alpha/(1.0+alpha*alpha)*(S*(dot(S, L)) - L);
        S_tilde[i->Index] = S + A*TimeStep + BL;
        SpinUpdate[i->Index] = 0.5*A*TimeStep + 0.5*BL;
    }
    //Now calculate the PlateUpdate and ExternalField to time+TimeStep;
    UpdateExternalField(Time+TimeStep);
    this->CalculateEffectiveField(AddTorque);
    for (std::vector<MagneticNode>::iterator i=NodeList.begin(); i!=NodeList.end(); i++)
    {
        double alpha = i->alpha;
        S = S_tilde[i->Index];
        L = this->RandomField[i->Index];
        Heff = this->EffectiveField[i->Index];
        A = -1.0/(1+alpha*alpha)*cross(S, Heff)-alpha/(1.0+alpha*alpha)*(S*(dot(S, Heff)) - Heff);
        BL = -1.0/(1+alpha*alpha)*cross(S, L)-alpha/(1.0+alpha*alpha)*(S*(dot(S, L)) - L);
        SpinUpdate[i->Index] = SpinUpdate[i->Index] + 0.5*A*TimeStep + 0.5*BL;
        if (i->Pinned == false)
            i->Spin = i->Spin + SpinUpdate[i->Index];
    }
}
/////
void SpinSystem::UpdateExternalField(double Time)
{

}
/////////
void SpinSystem::SetTemperature(double newTemperature)
{
    std::vector<MagneticNode>::iterator i;
    for (i=this->NodeList.begin(); i!=NodeList.end(); i++)
    {
        i->Temperature = newTemperature;
    }
}
//////////
void SpinSystem::ReadInElectronSiteIndex(const char* filename)
{
    FILE* fp;
    fp = fopen(filename, "r");
    int MagneticSiteIndex, ElectronSiteIndex;
    while (!feof(fp))
    {
        fscanf(fp, "%d%d\n", &MagneticSiteIndex, &ElectronSiteIndex);
        this->NodeList[MagneticSiteIndex].ElectronSiteIndex = ElectronSiteIndex;
        this->ListOfTorqueSiteIndecies.push_back(MagneticSiteIndex);
    }
    fclose(fp);
}
//////////
void SpinSystem::SetBackgroundField(vec BackgroundField)
{
    for (int i=0; i<this->NumSite; i++)
    {
        this->BackgroundField[i] = BackgroundField;
    }
}
///////////
void SpinSystem::CalculateTorque(ElectronSystem& Electrons, double Ef, double J_Hunds)
{
    for (int i=0; i<this->ListOfTorqueSiteIndecies.size(); i++)
    {
        int MagneticIndex = this->ListOfTorqueSiteIndecies[i];
        int ElectronIndex = this->NodeList[MagneticIndex].ElectronSiteIndex;
        this->NodeList[MagneticIndex].TorqueField = 2.0*(-J_Hunds)*Electrons.OnSiteSpin(ElectronIndex, Ef);
    }
}
///////////
void SpinSystem::OutputTextureToTextFile(const char* filename)
{
    FILE *fp;
    fp = fopen(filename, "w");
    std::vector<MagneticNode>::iterator i;
    printf("Generating the spin texture to txt file...");
    fprintf(fp, "x    y     z      Sx     Sy     Sz      JH      t     phi\n");
    for (i = this->NodeList.begin(); i!=this->NodeList.end(); i++)
    {
        fprintf(fp, "% lf\t% lf\t% lf\t% le\t% le\t% le\t-1.0000000\t-1.0000000\t0.0000000\n", 
                i->Location(0), 
                i->Location(1), 
                i->Location(2), 
                i->Spin(0), 
                i->Spin(1), 
                i->Spin(2));
    }
    fclose(fp);
    printf("done.\n");
}
//////////
void SpinSystem::OutputTorqueFieldToProFitTextFile(const char* filename)
{
    printf("Generating ProFit data for TorqueField plot...\n");
    FILE *fp;
    fp = fopen(filename, "w");
    double x, y, r, theta, sx, sy, sz;
    for (int i=0; i<this->NumSite; i++)
    {
        x = this->NodeList[i].Location(0);
        y = this->NodeList[i].Location(1);
        sx = this->NodeList[i].TorqueField(0);
        sy = this->NodeList[i].TorqueField(1);
        sz = this->NodeList[i].TorqueField(2);
        r = sqrt(sx*sx+sy*sy);
        theta = atan2(sy, sx);
        fprintf(fp, "% le\t% le\t% le\t% le\t% le\n", x, y, r, theta, sz);
    }
    fclose(fp);
}
////////////
void SpinSystem::OutputEffectiveToProFitTextFile(const char* filename)
{
    printf("Generating ProFit data for EffectiveField plot...\n");
    FILE *fp;
    fp = fopen(filename, "w");
    double x, y, r, theta, sx, sy, sz;
    vec Field(3);
    for (int i=0; i<this->NumSite; i++)
    {
        x = this->NodeList[i].Location(0);
        y = this->NodeList[i].Location(1);
        Field = this->EffectiveField[i];
        sx = Field(0);
        sy = Field(1);
        sz = Field(2);
        r = sqrt(sx*sx+sy*sy);
        theta = atan2(sy, sx);
        fprintf(fp, "% le\t% le\t% le\t% le\t% le\n", x, y, r, theta, sz);
    }
    fclose(fp);
}
/////////////////////////////
vec SpinSystem::SkyrmionLocation(void)
{
    std::vector<MagneticNode>::iterator i;
    vec SkyrmionPosition(3);
    SkyrmionPosition.zeros();
    double MinSz= 1.0;
    for (i=this->NodeList.begin(); i!=this->NodeList.end(); i++)
    {
        if (i->Spin(2) < MinSz)
        {
            MinSz = i->Spin(2);
            SkyrmionPosition = i->Location;
        }
    }
    return SkyrmionPosition;
}
/////////////////////////
void SpinSystem::OutputSpinTextureGIF(double Xmin, double Xmax, double Ymin, double Ymax, const char* title)
{
    Dislin TheWindow;
    TheWindow.errmod("all", "off");
    //char titleBuffer[256];
    float* xp;
    float* yp;
    float* xv;
    float* yv;
    xp = new float[NumSite];
    xv = new float[NumSite];
    yp = new float[NumSite];
    yv = new float[NumSite];
    //TheWindow.winsiz ((PlotX/2), PlotX/2);
    double YoverX_Ratio = (Ymax-Ymin)/(Xmax-Xmin);
    TheWindow.winsiz(600, (int)(600*YoverX_Ratio)+200);
    TheWindow.page (1800, (int)(1800*YoverX_Ratio)+500);
    TheWindow.sclmod ("full");
    TheWindow.scrmod ("revers");
    //Graph.scrmod ("black");
    TheWindow.metafl ("gif");
    TheWindow.x11mod("nostore");
    TheWindow.disini ();
    TheWindow.pagera ();
    TheWindow.hwfont ();
    TheWindow.axspos (400, -300);
    TheWindow.axslen (1200, (int)(1200*YoverX_Ratio));
    TheWindow.name ("X-axis", "x");
    TheWindow.name ("Y-axis", "y");
    TheWindow.vecopt(1.0, "scale");
    TheWindow.vecopt(1.0-0.618, "size");
    TheWindow.vecopt(1.0, "length");
    TheWindow.vecopt(18.0, "angle");
    //TheWindow.selwin(1);
    TheWindow.graf (Xmin-1, Xmax, 0, 10, Ymin-1, Ymax, 0, 10);
    TheWindow.height (50);
    TheWindow.vecclr (-2);   
    
    std::vector<MagneticNode>::iterator i;
    for (i=this->NodeList.begin(); i!=this->NodeList.end(); i++)
    {
        xp[i->Index] = (float)(i->Location(0));
        yp[i->Index] = (float)(i->Location(1));
        xv[i->Index] = (float)(i->Spin(0));
        yv[i->Index] = (float)(i->Spin(1));
    }    
    TheWindow.titlin (title, 4);
    TheWindow.title ();
    TheWindow.vecfld (xv, yv, xp, yp, this->NumSite, 1901);
    TheWindow.disfin();
}