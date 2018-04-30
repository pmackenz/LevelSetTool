#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "levelset.h"
#include "qcustomplot/qcustomplot.h"
#include <QVector>
#include "mycontours.h"

#ifdef USE_QWT

#include <qwt.h>
#include <qwt_plot.h>
#include <qwt_plot_canvas.h>
#include <qwt_plot_item.h>
#include <qwt_plot_curve.h>
#include <qwt_symbol.h>
#include <qwt_plot_spectrogram.h>
#include <qwt_matrix_raster_data.h>
#include <qwt_raster_data.h>

#endif



MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    driver = NULL;
    contours = NULL;

    ui->setupUi(this);

    ui->cbx_selectShape->addItem("Circle");
    ui->cbx_selectShape->addItem("Square");
    ui->cbx_selectShape->addItem("Star 1");
    ui->cbx_selectShape->addItem("Star 2");
    ui->cbx_selectShape->addItem("Irregular");

    ui->cbx_selectVelocityField->addItem("To the right");
    ui->cbx_selectVelocityField->addItem("To the top");
    ui->cbx_selectVelocityField->addItem("To upper right");
    ui->cbx_selectVelocityField->addItem("To bottom left");
    ui->cbx_selectVelocityField->addItem("Rotate CC");
    ui->cbx_selectVelocityField->addItem("Rotate CW");
    ui->cbx_selectVelocityField->addItem("Swirl");

    driver = new LevelSet();
    if (driver)
    {
        driver->setGridSize(GRIDSIZE, GRIDSIZE, NUMGRIDPOINTS,NUMGRIDPOINTS);
        //driver->setAlgorithm(TimeStepperType::FiniteDifference);
        driver->setAlgorithm(TimeStepperType::ControlVolume);
        driver->setShape(ui->cbx_selectShape->currentIndex());
        driver->setVelocityType(ui->cbx_selectVelocityField->currentIndex());

        this->setTimeStep(driver->getCFL());
    }

    ui->cbx_selectShape->setCurrentIndex(0);
    ui->cbx_selectVelocityField->setCurrentIndex(0);

    ui->chk_showGrid->setChecked(showGrid);
    ui->chk_showLabels->setChecked(showLabels);
    ui->chk_showLevelSet->setChecked(showLevelSet);

#ifdef USE_QWT
    plot = new QwtPlot(ui->plotWidget);
    plot->setCanvasBackground(QBrush(Qt::white));
#else
    plot = new QCustomPlot(ui->plotWidget);
#endif
    QLayout *lyt = ui->plotWidget->layout();
    lyt->addWidget(plot);

    refreshUI();
}

MainWindow::~MainWindow()
{
    delete ui;
    delete driver;
}

void MainWindow::refreshUI(void)
{
#ifdef USE_QWT
    this->QwtRefreshUI();
#else
    this->QCPRefreshUI();
#endif
}

#ifndef USE_QWT
void MainWindow::QCPRefreshUI(void)
{
    if (!driver) return;

    double *dim = driver->getSize();
    double lenx = dim[0];
    double leny = dim[1];

    QCustomPlot *customPlot = plot;

    customPlot->clearPlottables();
    customPlot->clearItems();

    customPlot->setLocale(QLocale(QLocale::English, QLocale::UnitedKingdom)); // period as decimal separator and comma as thousand separator
    customPlot->legend->setVisible(true);
    QFont legendFont = font();  // start out with MainWindow's font..
    legendFont.setPointSize(9); // and make a bit smaller for legend
    customPlot->legend->setFont(legendFont);
    customPlot->legend->setBrush(QBrush(QColor(255,255,255,230)));
    // by default, the legend is in the inset layout of the main axis rect. So this is how we access it to change legend placement:
    customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);

    // setup for graph 0: key axis left, value axis bottom
    customPlot->addGraph(customPlot->yAxis, customPlot->xAxis);
    customPlot->graph(0)->setPen(QPen(QColor(100, 100, 255)));
    // customPlot->graph(0)->setBrush(QBrush(QPixmap("./balboa.jpg"))); // fill with texture of specified image
    customPlot->graph(0)->setLineStyle(QCPGraph::lsNone);
    customPlot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssPlus, 5));
    customPlot->graph(0)->setName("Grid");

    // setup for graph 1: key axis bottom, value axis left (those are the default axes)
    // will contain bottom maxwell-like function with error bars

    //customPlot->addGraph();
    //customPlot->graph(1)->setPen(QPen(Qt::red));
    //customPlot->graph(1)->setBrush(QBrush(QPixmap("./balboa.jpg"))); // same fill as we used for graph 0
    //customPlot->graph(1)->setLineStyle(QCPGraph::lsLine);
    //customPlot->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssNone));
    //customPlot->graph(1)->setName("The Curve");
    //QCPErrorBars *errorBars = new QCPErrorBars(customPlot->xAxis, customPlot->yAxis);
    //errorBars->removeFromLegend();
    //errorBars->setDataPlottable(customPlot->graph(1));

    // plot the surface
    QList<QVector<double>> *curve = driver->getShape();

    if (curve)
    {
        QVector<double> x = (*curve)[0];
        QVector<double> y = (*curve)[1];

        QCPCurve* newCurve = new QCPCurve(customPlot->xAxis, customPlot->yAxis);
        newCurve->setName("initial curve");
        newCurve->addToLegend();
        newCurve->setPen(QPen(QColor(255,0,0)));
        newCurve->setData(x, y);
    }

    // plot the grid
    if (showGrid)
    {
        QList<QVector<double>> *grid = driver->getGrid();

        if (grid)
        {
            QVector<double> x = (*grid)[0];
            QVector<double> y = (*grid)[1];

            customPlot->graph(0)->setData(x, y);
        }
    }

    // plot level set
    // contours for level set function
    QVector<QVector<double> > *F = driver->getF();
    int mxSteps = driver->stepsX();
    if (driver->stepsY() > mxSteps) mxSteps = driver->stepsY();

    QVector<double> levels(4*mxSteps + 1);

    levels[0] = -2.*mxSteps;
    for (int i=1; i<= 4*mxSteps; i++)
        levels[i] = levels[i-1] + 1.0;

    if (showLevelSet)
    {
        if (contours) delete contours;
        //contours = new MyContours(X, Y, F, ui->plotWidget);
        //contours->setLevels(levels);
        if (showLabels)
        {
        //    contours->setLabels(true);
        }
    }

    //axes.contour(X, Y, F, [0.0], colors=('r'));

    // set ranges appropriate to show data:
    customPlot->xAxis->setRange(-lenx/2., lenx/2.);
    customPlot->yAxis->setRange(-leny/2., leny/2.);

    // set labels:
    customPlot->xAxis->setLabel("x");
    customPlot->yAxis->setLabel("y");

    // make ticks on bottom axis go outward:
    customPlot->xAxis->setTickLength(0, 5);
    customPlot->xAxis->setSubTickLength(0, 3);

    plot->replot();
}
#endif

#ifdef USE_QWT
void MainWindow::QwtRefreshUI(void)
{
    if (!driver) return;

    double *dim = driver->getSize();
    double lenx = dim[0];
    double leny = dim[1];

    QwtPlot *customPlot = plot;

    QwtPlotItemList list = customPlot->itemList();
    foreach (QwtPlotItem *item, list)
    {
        item->detach();
        delete item;
    }

    customPlot->setLocale(QLocale(QLocale::English, QLocale::UnitedStates)); // period as decimal separator and comma as thousand separator
    customPlot->legend();
    QFont legendFont = font();  // start out with MainWindow's font..
    legendFont.setPointSize(9); // and make a bit smaller for legend
    //customPlot->legend()->setFont(legendFont);
    //customPlot->legend()->setBrush(QBrush(QColor(255,255,255,230)));

    // by default, the legend is in the inset layout of the main axis rect. So this is how we access it to change legend placement:
    //customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);

    // setup for graph 0: key axis left, value axis bottom
    QwtPlotCurve *grid = new QwtPlotCurve();
    grid->setStyle(QwtPlotCurve::CurveStyle::NoCurve);
    QwtSymbol *theSymbol = new QwtSymbol(QwtSymbol::Cross);
    theSymbol->setSize(5);
    theSymbol->setColor(QColor(100, 100, 255));
    grid->setSymbol(theSymbol);
    grid->setTitle("Grid");

    grid->attach(customPlot);


    // setup for graph 1: key axis bottom, value axis left (those are the default axes)
    // will contain bottom maxwell-like function with error bars

    // plot the surface
    QList<QVector<double>> *curve = driver->getShape();

    if (curve)
    {
        QwtPlotCurve* newCurve = new QwtPlotCurve();
        //customPlot->xAxis, customPlot->yAxis);
        newCurve->setTitle("initial curve");

        newCurve->setPen(QPen(QColor(255,0,0)));
        newCurve->setSamples((*curve)[0], (*curve)[1]);

        newCurve->setItemAttribute(QwtPlotItem::Legend, false);
        newCurve->attach(customPlot);
    }

    // plot the grid
    if (showGrid)
    {
        QList<QVector<double>> *gridPoints = driver->getGrid();

        if (gridPoints)
        {
            QVector<double> x = (*gridPoints)[0];
            QVector<double> y = (*gridPoints)[1];

            grid->setSamples(x, y);
        }
    }

    if (showLevelSet)
    {
        // plot level set
        // contours for level set function
        QVector<QVector<double> > *F = driver->getF();

        QwtPlotSpectrogram *spectrogram = new QwtPlotSpectrogram(tr("Level set function"));
        spectrogram->setDisplayMode( QwtPlotSpectrogram::ContourMode, true );
        spectrogram->setDisplayMode( QwtPlotSpectrogram::ImageMode, false );

        int mxSteps = driver->stepsX();
        if (driver->stepsY() > mxSteps) mxSteps = driver->stepsY();

        QList<double> levels;
        levels.clear();
        levels.append(-2.*mxSteps);
        for (int i=1; i<= 4*mxSteps; i++)  { levels.append(levels.last() + 1.0); }
        spectrogram->setContourLevels( levels );

        QwtMatrixRasterData *specData = new QwtMatrixRasterData();
        int numCols = (*F)[0].length();
        int numRows = (*F).length();
        QVector<double> dataMatrix;
        foreach (QVector<double> vec, *F)
        {
            dataMatrix += vec;
        }
        specData->setValueMatrix(dataMatrix, numCols);
        specData->setInterval(Qt::Axis::XAxis, QwtInterval(-lenx/2, lenx/2));
        specData->setInterval(Qt::Axis::YAxis, QwtInterval(-leny/2, leny/2));
        specData->setResampleMode(QwtMatrixRasterData::BilinearInterpolation);

        spectrogram->setData(specData);

        spectrogram->attach(plot);

        //if (contours) delete contours;
        //contours = new MyContours(X, Y, F, ui->plotWidget);
        //contours->setLevels(levels);
        if (showLabels)
        {
            //    contours->setLabels(true);
        }
    }

    //axes.contour(X, Y, F, [0.0], colors=('r'));

    // set ranges appropriate to show data:

    plot->setAxisScale(QwtPlot::xBottom, -lenx/2., lenx/2.);
    plot->setAxisScale(QwtPlot::yLeft,   -leny/2., leny/2.);

    // set labels:

    plot->setAxisTitle(QwtPlot::xBottom, "x");
    plot->setAxisTitle(QwtPlot::yLeft,  "y");

    plot->replot();
}
#endif

/* *** slots *** */

void MainWindow::on_cbx_selectShape_currentIndexChanged(int index)
{
    if (driver)
    {
        driver->setShape(index);
        refreshUI();
    }
}

void MainWindow::on_cbx_selectVelocityField_currentIndexChanged(int index)
{
    if (driver)
    {
        driver->setVelocityType(index);
        refreshUI();
    }
}

void MainWindow::on_chk_showGrid_clicked(bool checked)
{
    showGrid = checked;
    refreshUI();
}

void MainWindow::on_chk_showLevelSet_clicked(bool checked)
{
    showLevelSet = checked;
    refreshUI();
}

void MainWindow::on_chk_showLabels_clicked(bool checked)
{
    showLabels = checked;
    refreshUI();
}

void MainWindow::on_theTimeStep_editingFinished()
{
    try {
        double deltaT = ui->theTimeStep->text().toDouble();
        if (deltaT < 0.0) setTimeStep(driver->getCFL());
        setTimeStep(deltaT);
    }
    catch (...) {
        setTimeStep(driver->getCFL());
    }
}

void MainWindow::setTimeStep(double dt)
{
    if (dt >= 0.0)
    {
        ui->theTimeStep->setText( QString("%1").arg(dt,7,'f',5) );
    }
}

void MainWindow::on_btn_setCFL100_clicked()
{
    setTimeStep(driver->getCFL());
}

void MainWindow::on_btn_setCFL050_clicked()
{
    setTimeStep(0.50*driver->getCFL());
}

void MainWindow::on_btn_setCFL025_clicked()
{
    setTimeStep(0.25*driver->getCFL());
}

void MainWindow::on_btn_singleStep_clicked()
{
    double deltaT = ui->theTimeStep->text().toDouble();
    driver->stepOne(deltaT);
    refreshUI();
}

void MainWindow::on_btn_10steps_clicked()
{
    double deltaT = ui->theTimeStep->text().toDouble();
    for (int i=0; i<10; i++)
    {
        driver->stepOne(deltaT);
        refreshUI();
    }
}

void MainWindow::on_btn_reset_clicked()
{
    driver->reset();
    setTimeStep(driver->getCFL());
    refreshUI();
}
