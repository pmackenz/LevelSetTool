#include "levelset.h"
#include "cmath"
#include <iostream>

#include "timestepper.h"
#include "finitedifference.h"
#include "controlvolume.h"

#include "algorithm"
#include "limits"

LevelSet::LevelSet(double lx, double ly, int n, int m)
{
    if (lx > 1.0) dim[0] = lx; else dim[0] = 1.0;     // side width of grid
    if (ly > 1.0) dim[1] = ly; else dim[1] = 1.0;     // side width of grid
    if (n > 10) Nsteps = n; else Nsteps = 10;      // number of grid cells per side
    if (m > 10) Msteps = m; else Msteps = 10;      // number of grid cells per side

    type   = 0;         // undefined

    shape.clear();

    timeStepper = new FiniteDifference();
    initGrid();
    setVelocityType(0);
}

LevelSet::~LevelSet()
{
    if (timeStepper) delete timeStepper;
}

void LevelSet::setAlgorithm(TimeStepperType type)
{
    if (timeStepper)
    {
        delete timeStepper;
        timeStepper = NULL;
    }

    switch (type)
    {
    case TimeStepperType::ControlVolume:
        timeStepper = new ControlVolume();
        break;
    case TimeStepperType::FiniteDifference:
        timeStepper = new FiniteDifference();
        break;
    }

    if (timeStepper)
    {
        timeStepper->setGridProperties(dim[0], dim[1], Nsteps, Msteps);
    }
}

void LevelSet::setGridSize(double lx, double ly)
{
    if (lx > 1.0) dim[0] = lx; else dim[0] = 1.0;     // side width of grid
    if (ly > 1.0) dim[1] = ly; else dim[1] = 1.0;     // side width of grid

    if (Nsteps < 10) Nsteps = 10;      // number of grid cells per side
    if (Msteps < 10) Msteps = 10;      // number of grid cells per side

    initGrid();
}

void LevelSet::setGridSize(double lx, double ly, int n, int m)
{
    if (lx > 1.0) dim[0] = lx; else dim[0] = 1.0;     // side width of grid
    if (ly > 1.0) dim[1] = ly; else dim[1] = 1.0;     // side width of grid

    if (n>10 && m>10)
    {
        Nsteps = n;
        Msteps = m;
    }

    initGrid();
}

void LevelSet::printGrid(void)
{

}

void LevelSet::initGrid(void)
{
    double dx = dim[0]/Nsteps;
    double dy = dim[1]/Msteps;
    double xi = -dim[0]/2.;
    double yj = -dim[1]/2.;

    QVector<double> xvals;
    for (int i=0; i<=Nsteps; i++) {
        xvals.append(xi);
        xi += dx;
    }
    gridX.clear();
    for (int j=0; j<=Msteps; j++) {
        gridX.append(xvals);
    }

    gridY.clear();
    for (int j=0; j<=Msteps; j++) {
         gridY.append(QVector<double>(Nsteps+1, yj));
         yj += dy;
    }
}

void LevelSet::setVelocityField(QVector<QVector<double> > &vField)
{

}

void LevelSet::setVelocityType(int type)
{
    if (type<0 || type>6) return;

    gridVx.clear();
    gridVy.clear();

    switch (type)
    {
    case 0:    // to the right
        for (int i=0; i<=Nsteps; i++) {
            gridVx.append(QVector<double>(Msteps, 1.0));
            gridVy.append(QVector<double>(Msteps, 0.0));
        }
        break;
    case 1:    // to the top
        for (int i=0; i<=Nsteps; i++) {
            gridVx.append(QVector<double>(Msteps, 0.0));
            gridVy.append(QVector<double>(Msteps, 1.0));
        }
        break;
    case 2:    // 45 degrees translation to upper right
        for (int i=0; i<=Nsteps; i++) {
            gridVx.append(QVector<double>(Msteps, 0.7071));
            gridVy.append(QVector<double>(Msteps, 0.7071));
        }
        break;
    case 3:    // 45 degrees translation to lower left
        for (int i=0; i<=Nsteps; i++) {
            gridVx.append(QVector<double>(Msteps,-0.7071));
            gridVy.append(QVector<double>(Msteps,-0.7071));
        }
        break;
    case 4:     // rigid rotation
        gridVx = QVector<QVector<double>>(Nsteps+1, QVector<double>(Msteps+1, 0.0));
        gridVy = QVector<QVector<double>>(Nsteps+1, QVector<double>(Msteps+1, 0.0));

        for (int i=0; i<=Nsteps; i++) {
            for (int j=0; j<=Msteps; j++) {
                double x = gridX[i][j];
                double y = gridY[i][j];
                double r = sqrt(x*x + y*y);
                // prevent division by zero
                if (r < 0.001)  { r = 0.001; }
                gridVx[i][j] = -0.10 * y;
                gridVy[i][j] =  0.10 * x;
            }
        }
        break;
    case 5:     // rigid rotation
        gridVx = QVector<QVector<double>>(Nsteps+1, QVector<double>(Msteps+1, 0.0));
        gridVy = QVector<QVector<double>>(Nsteps+1, QVector<double>(Msteps+1, 0.0));

        for (int i=0; i<=Nsteps; i++) {
            for (int j=0; j<=Msteps; j++) {
                double x = gridX[i][j];
                double y = gridY[i][j];
                double r = sqrt(x*x + y*y);
                // prevent division by zero
                if (r < 0.001)  { r = 0.001; }
                gridVx[i][j] =  0.10 * y;
                gridVy[i][j] = -0.10 * x;
            }
        }
        break;
    case 6:   // swirling rotation
        gridVx = QVector<QVector<double>>(Nsteps+1, QVector<double>(Msteps+1, 0.0));
        gridVy = QVector<QVector<double>>(Nsteps+1, QVector<double>(Msteps+1, 0.0));

        for (int i=0; i<=Nsteps; i++) {
            for (int j=0; j<=Msteps; j++) {
                double x = gridX[i][j];
                double y = gridY[i][j];
                double r = sqrt(x*x + y*y);
                // prevent division by zero
                if (r < 0.001)  { r = 0.001; }
                gridVx[i][j] = -y/r;
                gridVy[i][j] =  x/r;
            }
        }
        break;
    default:
        std::cout << "unknown velocity field code: " << type << std::endl;
    }
}

void LevelSet::setShape(int newType)
{
    double radius = 2.0;
    double x, y;

    if (newType < 0 ||newType > 4) return;   // not a known type -- do nothing

    if (newType != type || shape.size() < 2)
    {
        type = newType;

        shape.clear();

        switch (newType)
        {
        case 0:  // circle
            radius = 2.;
            for (int i=0; i<36; i++)
            {
                x = radius * cos(2.*M_PI*i/36.);
                y = radius * sin(2.*M_PI*i/36.);
                shape.append(QPointF(x,y));
            }
            shape.append(shape[0]);
            break;
        case 1: // square
            shape.append(QPointF( 0., 0.));
            shape.append(QPointF(-4., 0.));
            shape.append(QPointF(-4.,-4.));
            shape.append(QPointF( 0.,-4.));
            shape.append(QPointF( 0., 0.));
            break;
        case 2:  // star 1
            shape.append(QPointF( 0., 0.));
            shape.append(QPointF(-2.,-.5));
            shape.append(QPointF(-4., 0.));
            shape.append(QPointF(-3.5,-2.));
            shape.append(QPointF(-4.,-4.));
            shape.append(QPointF(-2.,-3.5));
            shape.append(QPointF( 0.,-4.));
            shape.append(QPointF(-0.5,-2.));
            shape.append(QPointF( 0., 0.));
            break;
        case 3:   // star 2
            shape.append(QPointF( 0., 0.));
            shape.append(QPointF(-4.,-1.));
            shape.append(QPointF(-8., 0.));
            shape.append(QPointF(-7.,-4.));
            shape.append(QPointF(-8.,-8.));
            shape.append(QPointF(-4.,-7.));
            shape.append(QPointF( 0.,-8.));
            shape.append(QPointF(-1.,-4.));
            shape.append(QPointF( 0., 0.));
            break;
        case 4:   // irregular
            shape.append(QPointF( 0., 0.));
            shape.append(QPointF(-4., 1.));
            shape.append(QPointF(-8.,-4.));
            shape.append(QPointF(-8.,-8.));
            shape.append(QPointF(-4.,-8.));
            shape.append(QPointF(-4.,-4.5));
            shape.append(QPointF( 1.,-4.5));
            shape.append(QPointF( 0., 0.));
            break;
        default:
            std::cout << "unknown shape type" << std::endl;
        }

        this->initF();
    }
}

void LevelSet::initF(void)
{
    if (shape.size() < 2) return;

    for (int k; k<=Nsteps; k++)
    {
        for (int l=0; l<=Msteps; l++)
        {
            gridF[k][l] = std::numeric_limits<double>::max();
        }
    }

    int npts = shape.size();

    QPointF xi = shape[0];
    QPointF xk;
    QPointF xj;

    for (int j=1; j<npts; j++)
    {
        QPointF xj = shape[j];

        // check if there is another point
        if (j + 1 < npts)
        {
            // we have another segment
            xk = shape[j + 1];
        }
        else if ( shape[0] == xj && npts > 2)
        {
            // we have a closed shape and the first segment is the next segment
            xk = shape[1];
        }
        else
        {
            // we have an open shape and reached the end.  Expand the shape
            // ... simply use d[m]
            xk = xj;
        }

        LINE_INFO val0 = getBase4line(xi, xj);
        double  dl0    = val0.dl;
        QPointF dlvec0 = val0.dlvec;
        QPointF nvec0  = val0.nvec;

        if (dl0 < 0.0) continue;

        LINE_INFO val1 = getBase4line(xj, xk);
        double  dl1    = val1.dl;
        QPointF dlvec1 = val1.dlvec;
        QPointF nvec1  = val1.nvec;

        double metric[2][2] = {
            QPointF::dotProduct(nvec1,nvec1), -QPointF::dotProduct(nvec0,nvec1),
            -QPointF::dotProduct(nvec1,nvec0),  QPointF::dotProduct(nvec0,nvec0)
        };
        double detMetric = metric[0][0]*metric[1][1] - metric[1][0]*metric[0][1];
        if (detMetric > 0.000001)
        {
            metric[0][0] /= detMetric;
            metric[0][1] /= detMetric;
            metric[1][0] /= detMetric;
            metric[1][1] /= detMetric;
        }
        else
        {
            metric[0][0] = 1.0; metric[1][1] = 1.0; metric[0][1] = 0.0; metric[1][0] = 0.0;
            std::cout << "escape singular metric at j = " << j << std::endl;
        }

        for (int k; k<=Nsteps; k++)
        {
            for (int l=0; l<=Msteps; l++)
            {
                double F = std::numeric_limits<double>::max();
                QPointF pt(gridX[k][l],gridY[k][l]);

                QPointF R = pt - xi;
                double d = QPointF::dotProduct(nvec0, R);
                double s = QPointF::dotProduct(dlvec0, R);

                if (s < -1.e-12) continue;  // numeric zero; outside segment

                if (s > 1.0)
                {
                    // we are sitting in a corner region
                    //
                    if (detMetric > 0.000001)
                    {
                        // we have a next segment
                        QPointF vec  = pt - xj;

                        double ds0 = QPointF::dotProduct(vec, nvec0);
                        double ds1 = QPointF::dotProduct(vec, nvec1);

                        double l0 = metric[0][0]*ds0 + metric[0][1]*ds1;
                        double l1 = metric[1][0]*ds0 + metric[1][1]*ds1;

                        double eps = 1.0e-5;

                        if ( l0 < eps and l1 < eps )
                        {
                            d = -sqrt(QPointF::dotProduct(vec,vec));
                        }
                        else if ( l0 > -eps and l1 > -eps )
                        {
                            d = sqrt(QPointF::dotProduct(vec,vec));
                        }
                        else
                        {
                            // we should not get here, or do we?
                            continue;
                        }
                    }
                }
                if (std::abs(d) < std::abs(gridF[k][l])) { gridF[k][l] = d; }
            }
        }

        xi = xj;
    }
}

LINE_INFO LevelSet::getBase4line(QPointF xi, QPointF xj)
{
    double dl;

    QPointF dlvec = xj - xi;
    QPointF nvec;

    double dl2 = QPointF::dotProduct(dlvec,dlvec);

    if (dl2 >= 0.0001)
    {
        dl = sqrt(dl2);
        nvec = QPointF(dlvec.y()/dl, -dlvec.x()/dl);
        dlvec /= dl2;
    }
    else
    {
        dl = -1.;
        nvec  = QPointF(0.0,0.0);
        dlvec = nvec;
    }

    LINE_INFO val;
    val.dl    = dl;
    val.dlvec = dlvec;
    val.nvec  = nvec;

    return val;
}

void LevelSet::stepOne(double deltaT)
{
    timeStepper->stepOne(deltaT);
}

void LevelSet::reset(void)
{
    setShape(type);
    initF();
}

QVector<QVector<double> > * LevelSet::getF(void)
{
    return &gridF;
}

QList<QVector<double> > * LevelSet::getGrid(void)
{
    QList<QVector<double>> *lists = new QList<QVector<double>>;
    QVector<double> x;
    QVector<double> y;

    x.clear();
    y.clear();

    for (int i=0; i<=Nsteps; i++)  {
        x += gridX[i];
        y += gridY[i];
    }

    lists->append(x);
    lists->append(y);

    return lists;
}

QList<QVector<double> > * LevelSet::getShape(void)
{
    QList<QVector<double>> *lists = new QList<QVector<double>>;
    QVector<double> x;
    QVector<double> y;

    x.clear();
    y.clear();

    foreach (QPointF pt, shape)
    {
        x.append(pt.x());
        y.append(pt.y());
    }

    lists->append(x);
    lists->append(y);

    return lists;
}

double LevelSet::getCFL(void)
{
    double cfl = 1.0;

    double lx = dim[0]/Nsteps;
    double ly = dim[1]/Msteps;

    double maxVx = 0.0001;
    foreach (QVector<double> vec, gridVx)
    {
        double max = *std::max_element(vec.constBegin(), vec.constEnd());
        double min = *std::min_element(vec.constBegin(), vec.constEnd());
        if ( max>maxVx) maxVx =  max;
        if (-min>maxVx) maxVx = -min;
    }

    double maxVy = 0.0001;
    foreach (QVector<double> vec, gridVy)
    {
        double max = *std::max_element(vec.constBegin(), vec.constEnd());
        double min = *std::min_element(vec.constBegin(), vec.constEnd());
        if ( max>maxVy) maxVy =  max;
        if (-min>maxVy) maxVy = -min;
    }

    if (maxVx > 0.00001 and lx/maxVx < cfl) cfl = lx/maxVx;
    if (maxVy > 0.00001 and ly/maxVy < cfl) cfl = ly/maxVy;

    return cfl;
}
