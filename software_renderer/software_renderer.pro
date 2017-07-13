#-------------------------------------------------
#
# Project created by QtCreator 2017-07-13T17:45:37
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = software_renderer
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

INCLUDEPATH += ../include

SOURCES += main.cpp\
        mainwindow.cpp \
    camera.cpp \
    ekran.cpp \
    light_source.cpp \
    matrix.cpp \
    model3d.cpp \
    pixbuf.cpp \
    renderer3d.cpp \
    structures.cpp \
    texture.cpp \
    vertex.cpp \
    world.cpp

HEADERS  += mainwindow.h \
    camera.h \
    ekran.h \
    light_source.h \
    matrix.h \
    model3d.h \
    pixbuf.h \
    renderer3d.h \
    structures.h \
    texture.h \
    ui_mainwindow.h \
    vertex.h \
    world.h

FORMS    += mainwindow.ui

QMAKE_CXXFLAGS += -O2 -larmadillo

unix|win32: LIBS += -L$$PWD/../libs/ -lliblapack

INCLUDEPATH += $$PWD/../libs
DEPENDPATH += $$PWD/../libs

unix|win32: LIBS += -L$$PWD/../libs/ -llibblas

INCLUDEPATH += $$PWD/../libs
DEPENDPATH += $$PWD/../libs
