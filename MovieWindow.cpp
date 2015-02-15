#include "MovieWindow.hpp"

MovieWindow::MovieWindow(double Xmin, double Xmax, double Ymin, double Ymax, int NumSite)
{
    this->PlotX = (int)(Xmax-Xmin);
    this->PlotY = (int)(Ymax-Ymin);
    this->xp = new float[NumSite];
    this->xv = new float[NumSite];
    this->yp = new float[NumSite];
    this->yv = new float[NumSite];
    this->TheTitle = "The movie window";
    //TheWindow.winsiz ((PlotX/2), PlotX/2);
    TheWindow.winsiz(600, 600);
    TheWindow.page (1800, 1800);
    TheWindow.sclmod ("full");
    TheWindow.scrmod ("revers");
    //Graph.scrmod ("black");
    TheWindow.metafl ("xwin");
    TheWindow.x11mod("nostore");
    TheWindow.disini ();
    TheWindow.pagera ();
    TheWindow.hwfont ();
    TheWindow.axspos (100, -100);
    TheWindow.axslen (1200, 1200);
    TheWindow.name ("X-axis", "x");
    TheWindow.name ("Y-axis", "y");
    TheWindow.vecopt(1.0, "scale");
    TheWindow.vecopt(1.0-0.618, "size");
    TheWindow.vecopt(1.0, "length");
    TheWindow.vecopt(18.0, "angle");
    TheWindow.selwin(1);
    TheWindow.graf (Xmin-1, Xmax, 0, 10, Ymin-1, Ymax, 0, 10);
    TheWindow.height (50);
    TheWindow.vecclr (-2);   
}
//////////////////////////////////////
MovieWindow::~MovieWindow()
{
    delete [] xp;
    delete [] xv;
    delete [] yp;
    delete [] yv;
    TheWindow.disfin();
}
////////////////////////////////////
void MovieWindow::UpdateWindow(SpinSystem& SpinTexture, string TitleUpdate)
{
    TheWindow.erase();
    const char* buffer = TitleUpdate.c_str();
    std::vector<MagneticNode>::iterator i;
    for (i=SpinTexture.NodeList.begin(); i!=SpinTexture.NodeList.end(); i++)
    {
        xp[i->Index] = (float)(i->Location(0));
        yp[i->Index] = (float)(i->Location(1));
        xv[i->Index] = (float)(i->Spin(0));
        yv[i->Index] = (float)(i->Spin(1));
    }    
    TheWindow.titlin (buffer, 4);
    TheWindow.title ();
    TheWindow.vecfld (xv, yv, xp, yp, SpinTexture.NumSite, 1901);
}
