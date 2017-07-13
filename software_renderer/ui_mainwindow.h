/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.8.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QSlider>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>
#include "ekran.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    Ekran *ekran;
    QGroupBox *groupBox;
    QSlider *trX;
    QSlider *trY;
    QSlider *trZ;
    QSlider *rotX;
    QSlider *rotY;
    QSlider *rotZ;
    QSlider *shX;
    QSlider *shY;
    QSlider *shZ;
    QLabel *label;
    QLabel *label_2;
    QLabel *label_3;
    QLabel *label_4;
    QSlider *scX;
    QSlider *scY;
    QSlider *scZ;
    QSlider *nrSlider;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QStringLiteral("MainWindow"));
        MainWindow->resize(1599, 795);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        ekran = new Ekran(centralWidget);
        ekran->setObjectName(QStringLiteral("ekran"));
        ekran->setGeometry(QRect(10, 10, 1280, 720));
        groupBox = new QGroupBox(centralWidget);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        groupBox->setGeometry(QRect(1310, 10, 271, 721));
        trX = new QSlider(groupBox);
        trX->setObjectName(QStringLiteral("trX"));
        trX->setGeometry(QRect(10, 70, 251, 22));
        trX->setMinimum(-200);
        trX->setMaximum(200);
        trX->setOrientation(Qt::Horizontal);
        trY = new QSlider(groupBox);
        trY->setObjectName(QStringLiteral("trY"));
        trY->setGeometry(QRect(10, 110, 251, 22));
        trY->setMinimum(-200);
        trY->setMaximum(200);
        trY->setOrientation(Qt::Horizontal);
        trZ = new QSlider(groupBox);
        trZ->setObjectName(QStringLiteral("trZ"));
        trZ->setGeometry(QRect(10, 150, 251, 22));
        trZ->setMinimum(-200);
        trZ->setMaximum(200);
        trZ->setOrientation(Qt::Horizontal);
        rotX = new QSlider(groupBox);
        rotX->setObjectName(QStringLiteral("rotX"));
        rotX->setGeometry(QRect(10, 230, 251, 22));
        rotX->setMaximum(359);
        rotX->setOrientation(Qt::Horizontal);
        rotY = new QSlider(groupBox);
        rotY->setObjectName(QStringLiteral("rotY"));
        rotY->setGeometry(QRect(9, 270, 251, 22));
        rotY->setMaximum(359);
        rotY->setOrientation(Qt::Horizontal);
        rotZ = new QSlider(groupBox);
        rotZ->setObjectName(QStringLiteral("rotZ"));
        rotZ->setGeometry(QRect(9, 310, 251, 22));
        rotZ->setMaximum(359);
        rotZ->setOrientation(Qt::Horizontal);
        shX = new QSlider(groupBox);
        shX->setObjectName(QStringLiteral("shX"));
        shX->setGeometry(QRect(50, 390, 160, 22));
        shX->setMinimum(-100);
        shX->setMaximum(100);
        shX->setOrientation(Qt::Horizontal);
        shY = new QSlider(groupBox);
        shY->setObjectName(QStringLiteral("shY"));
        shY->setGeometry(QRect(50, 430, 160, 22));
        shY->setMinimum(-100);
        shY->setMaximum(100);
        shY->setOrientation(Qt::Horizontal);
        shZ = new QSlider(groupBox);
        shZ->setObjectName(QStringLiteral("shZ"));
        shZ->setGeometry(QRect(50, 470, 160, 22));
        shZ->setMinimum(-100);
        shZ->setMaximum(100);
        shZ->setOrientation(Qt::Horizontal);
        label = new QLabel(groupBox);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(50, 40, 47, 13));
        label_2 = new QLabel(groupBox);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setGeometry(QRect(50, 200, 47, 13));
        label_3 = new QLabel(groupBox);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setGeometry(QRect(50, 360, 47, 13));
        label_4 = new QLabel(groupBox);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setGeometry(QRect(50, 530, 47, 13));
        scX = new QSlider(groupBox);
        scX->setObjectName(QStringLiteral("scX"));
        scX->setGeometry(QRect(50, 560, 160, 22));
        scX->setMinimum(-100);
        scX->setMaximum(100);
        scX->setOrientation(Qt::Horizontal);
        scY = new QSlider(groupBox);
        scY->setObjectName(QStringLiteral("scY"));
        scY->setGeometry(QRect(50, 600, 160, 22));
        scY->setMinimum(-100);
        scY->setMaximum(100);
        scY->setOrientation(Qt::Horizontal);
        scZ = new QSlider(groupBox);
        scZ->setObjectName(QStringLiteral("scZ"));
        scZ->setGeometry(QRect(50, 640, 160, 22));
        scZ->setMinimum(-100);
        scZ->setMaximum(100);
        scZ->setOrientation(Qt::Horizontal);
        nrSlider = new QSlider(groupBox);
        nrSlider->setObjectName(QStringLiteral("nrSlider"));
        nrSlider->setGeometry(QRect(10, 10, 251, 22));
        nrSlider->setMaximum(0);
        nrSlider->setOrientation(Qt::Horizontal);
        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1599, 21));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        MainWindow->setStatusBar(statusBar);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", Q_NULLPTR));
        groupBox->setTitle(QString());
        label->setText(QApplication::translate("MainWindow", "Tranformation", Q_NULLPTR));
        label_2->setText(QApplication::translate("MainWindow", "Rotation", Q_NULLPTR));
        label_3->setText(QApplication::translate("MainWindow", "Shearnig", Q_NULLPTR));
        label_4->setText(QApplication::translate("MainWindow", "Scale", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
