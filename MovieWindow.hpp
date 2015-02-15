/* 
 * File:   MovieWindow.hpp
 * Author: jason
 *
 * Created on February 14, 2015, 6:03 PM
 */

#ifndef MOVIEWINDOW_HPP
#define	MOVIEWINDOW_HPP

#include <armadillo>
#include <iostream>
#include <string>
#include "discpp.h"
#include "SpinSystem.hpp"
#include "discpp.h"

using namespace arma;
using namespace std;

class MovieWindow
{
public:
    Dislin TheWindow;
    float* xp;
    float* yp;
    float* xv;
    float* yv;
    int PlotX;
    int PlotY;
    string TheTitle;
    MovieWindow(double Xmin, double Xmax, double Ymin, double Ymax, int NumSite);
    ~MovieWindow();
    void UpdateWindow(SpinSystem &SpinTexture, string TitleUpdate);
};

#endif	/* MOVIEWINDOW_HPP */

