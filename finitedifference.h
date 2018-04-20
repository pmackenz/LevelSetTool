#ifndef FINITEDIFFERENCE_H
#define FINITEDIFFERENCE_H

#include "timestepper.h"

class FiniteDifference : public TimeStepper
{
public:
    FiniteDifference();
    void stepOne(double deltaT = 1.0);
};

#endif // FINITEDIFFERENCE_H
