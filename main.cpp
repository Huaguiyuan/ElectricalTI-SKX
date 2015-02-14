/* 
 * File:   main.cpp
 * Author: RAF_Jason
 *
 * Created on 2013年10月1日, 上午10:01
 */

#include "Calculations.hpp"

int main (int argc, char** argv )
{
    SpinSystem Haha("input.txt",1.0, 1.0);
    printf("test.\n");
    vec ax, ay, az;
    ax << 31.0 << 0.0 << 0.0;
    ay << 0.0 << 31.0 << 0.0;
    az << 0.0 << 0.0 << 31.0;
    Haha.GenerateNeighborList(1.00001, 1.00001, true, false, false, ax, ay, az);
    std::vector<MagneticNode>::iterator i;
    for (i=Haha.NodeList.begin(); i!=Haha.NodeList.end(); i++)
    {
        printf("Index=%d, Number of neighbours=%lu\n", i->Index, i->ListOfExchangeNeighbours.size());
        for (std::vector<int>::iterator j=i->ListOfExchangeNeighbours.begin(); j!=i->ListOfExchangeNeighbours.end(); j++)
        {
            printf("    %d", *j);
        }
        printf("\n");
    }
    vec BackgroundField(3);
    BackgroundField << 0.0 << 0.0 << 1.0;
    for (int i=0; i<Haha.NumSite; i++)
    {
        Haha.ExternalField[i] = Haha.ExternalField[i] + BackgroundField;
    }
    Haha.CalculateEffectiveField();
    for (int i=0; i<Haha.NumSite; i++)
    {
        printf("%lf\n", (Haha.EffectiveField[i])(2));  
        
    }
}
    