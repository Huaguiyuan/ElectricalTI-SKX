/* 
 * File:   NeighborList.hpp
 * Author: RAF_Jason
 *
 * Created on 2013年10月3日, 下午4:06
 */

#ifndef NEIGHBORLIST_HPP
#define	NEIGHBORLIST_HPP
#include <armadillo>
struct NeighborRecord
{
    int Right[3];
    int Left[3];
    int Up[3];
    int Down[3];
    int Top[3];
    int Bottom[3];
};
using namespace arma;

NeighborRecord*** GenerateNeighbors(int SizeX, int SizeY, int SizeZ);
void DeleteNeighbors(NeighborRecord*** p, int SizeX, int SizeY);
void NeighborTest(void);

#endif	/* NEIGHBORLIST_HPP */
