#ifndef LEVELSET_H
#define LEVELSET_H

#include <QList>
#include <QVector>
#include "timestepper.h"
class QPointF;

struct LINE_INFO {
    double  dl;
    QPointF dlvec;
    QPointF nvec;
};

class LevelSet
{
public:
    LevelSet(double lx=10., double ly=10., int n=10, int m=10);
    ~LevelSet();
    void setAlgorithm(TimeStepperType type);
    void setGridSize(double lx, double ly);
    void setGridSize(double lx, double ly, int n, int m);
    void printGrid(void);
    void initGrid(void);
    void setVelocityField(QVector<QVector<double> > &vField);
    void setVelocityType(int type=0);
    void setShape(int type=0);
    void initF(void);
    LINE_INFO getBase4line(QPointF xi, QPointF xj);
    void stepOne(double deltaT);
    void reset(void);
    QVector<QVector<double> > & getF(void);
    QList<QVector<double> > * getGrid(void);
    QList<QVector<double> > * getShape(void);
    int stepsX(void) { return Nsteps; };
    int stepsY(void) { return Msteps; };
    double getCFL(void);
    double *getSize() { return dim; };

private:
    double len;   // side width of grid
    int Nsteps;   // number of grid cells per side
    int Msteps;   // number of grid cells per side
    double dim[2];
    int type;

    QList<QPointF> shape;
    QVector<QVector<double>> gridF;
    QVector<QVector<double>> gridX;
    QVector<QVector<double>> gridY;
    QVector<QVector<double>> gridVx;
    QVector<QVector<double>> gridVy;

    TimeStepper * timeStepper;

};

#endif // LEVELSET_H
