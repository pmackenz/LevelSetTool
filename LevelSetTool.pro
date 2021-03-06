#-------------------------------------------------
#
# Project created by QtCreator 2018-03-16T10:15:03
#
#-------------------------------------------------

QT += core gui printsupport

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = LevelSetTool
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

INCLUDEPATH += /Library/Frameworks/qwt.framework/Headers/

#include(/Library/Frameworks/qwt.framework/qwt.prl)
include( /usr/local/qwt-6.2.0-svn/features/qwtconfig.pri )
include( /usr/local/qwt-6.2.0-svn/features/qwt.prf )


SOURCES += \
        main.cpp \
        mainwindow.cpp \
    levelset.cpp \
    timestepper.cpp \
    finitedifference.cpp \
    controlvolume.cpp \
    qcustomplot/qcustomplot.cpp \
    mycontours.cpp

HEADERS += \
        mainwindow.h \
    levelset.h \
    timestepper.h \
    finitedifference.h \
    controlvolume.h \
    qcustomplot/qcustomplot.h \
    mycontours.h

FORMS += \
        mainwindow.ui
