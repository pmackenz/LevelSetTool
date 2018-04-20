#include "finitedifference.h"

FiniteDifference::FiniteDifference()
    : TimeStepper()
{

}

void FiniteDifference::stepOne(double deltaT)
{
    double dF;

    double lx = lenx / Nsteps;
    double ly = leny / Msteps;

    // x-direction transport stage
    for (int j=0; j <= Msteps; j++ )
    {
        for (int i=0; i<=Nsteps; i++)
        {
            double vx = (*gridVx)[i][j];

            if (vx > 0.0)
            {
                if (i>0)
                    { dF = (*gridF)[i][j] - (*gridF)[i-1][j]; }
                else
                    { dF = (*gridF)[i][j] - (3.*(*gridF)[0][j] - 3.*(*gridF)[1][j]   + 1.*(*gridF)[2][j]); }
            }
            else
            {
                if (i<Nsteps)
                    { dF = (*gridF)[i+1][j] - (*gridF)[i][j]; }
                else
                    { dF = (3.*(*gridF)[i][j] - 3.*(*gridF)[i-1][j] + 1.*(*gridF)[i-2][j]) - (*gridF)[i][j]; }
            }
            (*gridF)[i][j] -= (deltaT * vx / lx) * dF;
        }
    }

    // y-direction transport stage
    for (int i=0; i<=Nsteps; i++)
    {
        //
        // now take care of the interior
        for (int j=0; j<=Msteps; j++)
        {
            double vy = (*gridVy)[i][j];

            if (vy > 0.0)
            {
                if (j>0)
                    { dF = (*gridF)[i][j] - (*gridF)[i][j-1]; }
                else
                    { dF = (*gridF)[i][j] - (3.*(*gridF)[i][0] - 3.*(*gridF)[i][1]   + 1.0*(*gridF)[i][2]); }
            }
            else
            {
                if (j<Msteps)
                    { dF = (*gridF)[i][j+1] - (*gridF)[i][j]; }
                else
                    { dF = (3.*(*gridF)[i][j] - 3.*(*gridF)[i][j-1] + 1.0*(*gridF)[i][j-2]) - (*gridF)[i][j]; }
            }
            (*gridF)[i][j] -= (deltaT * vy / ly) * dF;
        }
    }
}
