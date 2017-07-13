#include "mainwindow.h"
#include <QApplication>

#include <iostream>

int main(int argc, char *argv[])
{
    //std::ios::sync_with_stdio(0);
    QApplication a(argc, argv);
    World world;
    world.load_assets();

    MainWindow w(&world);
    w.show();
    //w.loadAssets();

    return a.exec();
}
