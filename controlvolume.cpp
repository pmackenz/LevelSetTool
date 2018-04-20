#include "controlvolume.h"
#include <QVector>
#include <iostream>

ControlVolume::ControlVolume()
    : TimeStepper()
{

}

void ControlVolume::updateGrid()
{
    foreach (QVector<double> vec, dF)
    {
        vec.clear();
    }
    dF.clear();

    dF.resize(Nsteps+1);
    foreach (QVector<double> vec, dF) {
        vec.resize(Msteps);
        vec.fill(0.0);
    }

    foreach (QVector<double> vec, dF) {
        std::cout << vec.size() << "..";
    }
    std::cout << std::endl;
}

void ControlVolume::stepOne(double deltaT)
{
    foreach (QVector<double> vec, dF) {
        vec.resize(Msteps);
        vec.fill(0.0);
    }

    double lx = lenx / Nsteps;
    double ly = leny / Msteps;

    // x-direction transport stage
    for (int j=0; j<=Msteps; j++)
    {

        // handle the left boundary first
        int n = Nsteps;    //  + 1 for number of points, - 1 for indexing 0..n-1

        // x-flux
        double FexpandLeft  = 1.875*(*gridF)[0][j] - 1.25*(*gridF)[1][j]   + 0.375*(*gridF)[2][j];
        double FexpandRight = 1.875*(*gridF)[n][j] - 1.25*(*gridF)[n-1][j] + 0.375*(*gridF)[n-2][j];

        // now take care of the interior
        for (int i=0; i<=Nsteps; i++)
        {
            double vx;
            double flxX;

            if (i<n)
                { vx = 0.50 * ((*gridVx)[i][j] + (*gridVx)[i+1][j]); }
            else
                { vx = (*gridVx)[n][j]; }

            if (vx >= 0.0)
                { flxX = ly * (*gridF)[i][j] * vx; }
            else
            {
                if (i<n)
                    { flxX = ly * (*gridF)[i+1][j] * vx; }
                else
                    { flxX = ly * FexpandRight * (*gridVx)[n][j]; }
            }

            dF[i][j] -= flxX;
            if (i==0) { dF[i][j] += ly * FexpandLeft * (*gridVx)[0][j]; }
            if (i<n)  { dF[i+1][j] += flxX; }
        }
    }
    // y-direction transport stage
    for (int i=0; i<=Nsteps; i++)
    {
        //
        // handle the left boundary first
        int m = Msteps;    //  + 1 for number of points, - 1 for indexing 0..n-1

        // y-flux
        double FexpandBottom  = 1.875*(*gridF)[i][0] - 1.25*(*gridF)[i][1]   + 0.375*(*gridF)[i][2];
        double FexpandTop     = 1.875*(*gridF)[i][m] - 1.25*(*gridF)[i][m-1] + 0.375*(*gridF)[i][m-2];

        // now take care of the interior
        for (int j=0; j<=Msteps; j++)
        {
            double vy;
            double flxY;

            if (j<m)
                { vy = 0.50 * ((*gridVy)[i][j] + (*gridVy)[i][j+1]); }
            else
                { vy = (*gridVy)[i][m];}

            if (vy >= 0.0)
                { flxY = lx * (*gridF)[i][j] * vy; }
            else
            {
                if (j<m)
                    { flxY = lx * (*gridF)[i][j+1] * vy; }
                else
                    { flxY = lx * FexpandTop * vy; }
            }

            dF[i][j] -= flxY;
            if (j==0)  { dF[i][j] += lx * FexpandBottom * (*gridVy)[i][0]; }
            if (j<m)   { dF[i][j+1] += flxY;}
        }
    }
    double mult = deltaT / lx / ly;

    for (int i=0; i<=Nsteps; i++)
    {
        for (int j=0; j<=Msteps; j++)
        {
            gridF[i][j] += mult * dF[i][j];
        }
    }
}
