#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(World* world, QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    world(world)
{
    ui->setupUi(this);

    connect((ui->nrSlider), &QSlider::valueChanged, [=](){ set_model_nr(ui->nrSlider->value()); });

    connect((ui->trX), &QSlider::valueChanged, [=](){ set_tr_x(ui->trX->value()); });
    connect((ui->trY), &QSlider::valueChanged, [=](){ set_tr_y(ui->trY->value()); });
    connect((ui->trZ), &QSlider::valueChanged, [=](){ set_tr_z(ui->trZ->value()); });

    connect((ui->rotX), &QSlider::valueChanged, [=](){ set_rot_x(ui->rotX->value()); });
    connect((ui->rotY), &QSlider::valueChanged, [=](){ set_rot_y(ui->rotY->value()); });
    connect((ui->rotZ), &QSlider::valueChanged, [=](){ set_rot_z(ui->rotZ->value()); });

    connect((ui->shX), &QSlider::valueChanged, [=](){ set_sh_x(ui->shX->value()); });
    connect((ui->shY), &QSlider::valueChanged, [=](){ set_sh_y(ui->shY->value()); });
    connect((ui->shZ), &QSlider::valueChanged, [=](){ set_sh_z(ui->shZ->value()); });

    connect((ui->scX), &QSlider::valueChanged, [=](){ set_sc_x(ui->scX->value()); });
    connect((ui->scY), &QSlider::valueChanged, [=](){ set_sc_y(ui->scY->value()); });
    connect((ui->scZ), &QSlider::valueChanged, [=](){ set_sc_z(ui->scZ->value()); });

    connect((ui->lightA), &QSlider::valueChanged, [=](){ set_l_a(ui->lightA->value()); });
    connect((ui->lightB), &QSlider::valueChanged, [=](){ set_l_b(ui->lightB->value()); });
    connect((ui->lightC), &QSlider::valueChanged, [=](){ set_l_c(ui->lightC->value()); });

    model_nr = 0;
    ui->nrSlider->setMaximum(world->models.size());

    ui->ekran->setWorld(world);
}

MainWindow::~MainWindow() {
    delete ui;
}

//void MainWindow::loadAssets()
//{
//    ui->ekran->loadAssets();
//    setModelCount();
//}

void MainWindow::keyPressEvent(QKeyEvent* event) {
    ui->ekran->handleKeyPress(event);
}

void MainWindow::keyReleaseEvent(QKeyEvent* event) {
    ui->ekran->handleKeyRelease(event);
}

//void MainWindow::setModelCount() {
//    ui->nrSlider->setMaximum(ui->ekran->getModelCount());
//}
