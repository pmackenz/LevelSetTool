#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#define GRIDSIZE 20.0
#define NUMGRIDPOINTS 40

class LevelSet;
class MyContours;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    void refreshUI(void);

private slots:
    void on_cbx_selectShape_currentIndexChanged(int index);
    void on_cbx_selectVelocityField_currentIndexChanged(int index);
    void on_chk_showGrid_clicked(bool checked);
    void on_chk_showLevelSet_clicked(bool checked);
    void on_chk_showLabels_clicked(bool checked);
    void on_theTimeStep_editingFinished();
    void on_btn_setCFL100_clicked();
    void on_btn_setCFL050_clicked();
    void on_btn_setCFL025_clicked();
    void on_btn_singleStep_clicked();
    void on_btn_10steps_clicked();
    void on_btn_reset_clicked();

private:
    void setTimeStep(double dt);
    Ui::MainWindow *ui;
    LevelSet *driver;
    MyContours *contours;

    bool showGrid     = true;
    bool showLabels   = true;
    bool showLevelSet = true;
};

#endif // MAINWINDOW_H
