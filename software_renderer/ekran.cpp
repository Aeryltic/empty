#include "ekran.h"

#include <QMouseEvent>
#include <QPainter>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stack>

using namespace std;
Ekran::Ekran(QWidget *parent) : QWidget(parent), width(1280), height(720), renderer(width, height, 0xff000000) {
    pressed = false;
    world = nullptr;
}

Ekran::~Ekran() {

}

void Ekran::setWorld(World *world) {
    this->world = world;
}

void Ekran::paintEvent(QPaintEvent *event) {
    QPainter p(this);
    //p.fillRect(0, 0, width(), height(), Qt::white);
    if(world)
        renderer.render(world);
    QImage im((uchar*)renderer.get_buffer(), width, height, QImage::Format_ARGB32);
    p.drawImage(0,0,im);
}

void Ekran::mousePressEvent(QMouseEvent *event) {
    if(event->buttons()==Qt::LeftButton) {
        pressed = true;
        mx = event->x();
        my = event->y();
    }
//    if(event->buttons()==Qt::RightButton){
//        world->cam.v = 100;
//    }
    update();
}

void Ekran::mouseMoveEvent(QMouseEvent *event) {
    if(event->buttons()==Qt::LeftButton) {
        if(pressed) {
            world->cam.v = 0;
            int nx = event->x(),
                ny = event->y();
            world->cam.ry -= double(nx - mx) / width * M_PI;
            world->cam.rx += double(ny - my) / height * M_PI;
            mx = nx;
            my = ny;
        }
    }
    update();
}

void Ekran::mouseReleaseEvent(QMouseEvent *event) {
//    world->cam.v = 0;
    pressed = false;
    //if(event->buttons()==Qt::LeftButton) pressed = false;
//    if(event->buttons()==Qt::RightButton){
//        world->cam.v = 0;
//    }
    update();
}

void Ekran::wheelEvent(QWheelEvent *event) {
    world->cam.pz += event->delta() / 20.0;
    update();
}


void Ekran::handleKeyPress(QKeyEvent *event) {
//    if(event->key()==Qt::Key_W){
//        world->cam.v = 100;
//    }

//    world->cam.update(10);
}

void Ekran::handleKeyRelease(QKeyEvent *event) {
//    if(event->key()==Qt::Key_W){
//        world->cam.v = 0;
//    }

    //    world->cam.update(10);
}
