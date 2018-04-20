#ifndef CONTROLVOLUME_H
#define CONTROLVOLUME_H

#include "timestepper.h"

class ControlVolume : public TimeStepper
{
public:
    ControlVolume();
    void stepOne(double deltaT = 1.0);

protected:
    void updateGrid();

    QVector<QVector<double> > dF;
};

#endif // CONTROLVOLUME_H
