#ifndef MYCONTOURS_H
#define MYCONTOURS_H

#include <QVector>
#include <QColor>

class QCustomPlot;
class QCPCurve;

class MyContours
{
public:
    MyContours(QVector<QVector<double> > &X, QVector<QVector<double> > &Y, QVector<QVector<double> > &F, QCustomPlot *);
    ~MyContours(void);
    void setLevels(QVector<double> &levels);
    void setLevels(QVector<double> &levels, QVector<QColor> &colors);
    void setLabels(bool );
    void refresh(void);

private:
    void findContours(void);

protected:
    QCustomPlot *mPlot;

    bool showLabels = false;
    QVector<double> mLevels;
    QVector<QColor> mColors;
    QList<QCPCurve *> mCurves;

    QVector< QVector< double > > mX;
    QVector< QVector< double > > mY;
    QVector< QVector< double > > mF;

    int Nsteps;
    int Msteps;

};

#endif // MYCONTOURS_H
