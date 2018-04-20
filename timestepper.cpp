#include "timestepper.h"
#include "cassert"

TimeStepper::TimeStepper()
{
    gridF = NULL;
    gridVx = NULL;
    gridVy = NULL;
}

void TimeStepper::setGridProperties(double lx, double ly, int N, int M)
{
    assert(lx>0.0);
    lenx = lx;
    assert(ly>0.0);
    leny = ly;
    assert(N>=10);
    Nsteps = N;
    assert(M>=10);
    Msteps = M;

    this->updateGrid();
}
