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
    for (int i=0; i<Haha.NumSite; i++)
    {
        printf("Index=%d, Number of neighbours=%lu\n", Haha.NodeList[i].Index, Haha.NodeList[i].ListOfExchangeNeighbours.size());
        for (int j=0; j<Haha.NodeList[i].ListOfExchangeNeighbours.size(); j++)
        {
            printf("    %d", Haha.NodeList[i].ListOfExchangeNeighbours[j]);
        }
        printf("\n");
    }
}
    