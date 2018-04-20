#ifndef TIMESTEPPER_H
#define TIMESTEPPER_H

#include <QVector>

enum class TimeStepperType { ControlVolume, FiniteDifference };

class TimeStepper
{
public:
    TimeStepper();
    virtual void stepOne(double deltaT = 1.0) = 0;
    virtual void setGridProperties(double lx, double ly, int N, int M);
    virtual void setF ( QVector<QVector<double> > &newF ) { gridF  = &newF;  };
    virtual void setVx( QVector<QVector<double> > &newVx) { gridVx = &newVx; };
    virtual void setVy( QVector<QVector<double> > &newVy) { gridVy = &newVy; };

    QVector<QVector<double> > * getF() { return gridF; };

protected:
    virtual void updateGrid() {};

    QVector<QVector<double> > * gridF;
    QVector<QVector<double> > * gridVx;
    QVector<QVector<double> > * gridVy;

    double lenx;
    double leny;
    int Msteps;
    int Nsteps;
};

#endif // TIMESTEPPER_H
