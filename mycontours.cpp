#include "mycontours.h"
#include <QVector>
#include <QColor>
#include <QDebug>
#include "qcustomplot/qcustomplot.h"

MyContours::MyContours(QVector<QVector<double> > &X, QVector<QVector<double> > &Y, QVector<QVector<double> > &F, QCustomPlot *plot)
{
    mPlot = plot;

    mCurves.clear();
    //new QCPCurve();

    Nsteps = X.size();
    if (Nsteps>1)
        Msteps = X[0].size();
    else
        Msteps = 0;
    if (Nsteps==0 || Msteps==0)
    {
        // friendly error message and die
        qWarning() << "too few grid points for contour plot";
    }
    else
    {
        mX = X;
        mY = Y;
        mF = F;
        // add size checking for mX, mY, and mF
    }
}

MyContours::~MyContours(void)
{

}

void MyContours::setLevels(QVector<double> &levels)
{
    if (levels.size() > 1)
    {
        //mLevels.clear();
        mLevels = levels;
        mColors = QVector<QColor>(levels.size(), QColor(0,0,0));
    }
}

void MyContours::setLevels(QVector<double> &levels, QVector<QColor> &colors)
{
    if (levels.size()>0)
        mLevels = levels;
    else
        mLevels.clear();

    if (colors.size()>0)
    {
        mColors.clear();
        if (colors.size() == levels.size())
            mColors = colors;
        if (colors.size() > levels.size())
        {
            mColors.resize(levels.size());
        }
        if (colors.size() < levels.size())
        {
            int j=0;
            for (int k=0; k<levels.size(); k++)
            {
                mColors[k] = colors[j];
                if (++j == colors.size()) j=0;
            }
        }
    }
    else
    {
        mColors = QVector<QColor>(levels.size(), QColor(0,0,0));
    }
}

void MyContours::setLabels(bool show=true)
{
    if (show != showLabels)
    {
        showLabels = show;
        refresh();
    }
}

void MyContours::refresh(void)
{

}

void MyContours::findContours(void)
{
    if (mLevels.size()>0)
    {
        foreach (QCPCurve *curve, mCurves) {
            if (curve)
            {
                delete curve;
                curve = NULL;
            }
        }
        mCurves.clear();

        // loop through cells
        for (int i=0; i<Nsteps; i++)
        {
            for (int j=0; j<Msteps; j++)
            {
                //
            }
        }
    }
}
