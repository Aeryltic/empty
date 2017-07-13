#ifndef EKRAN_H
#define EKRAN_H

#include <QWidget>
#include <mainwindow.h>

#include "world.h"
#include "renderer3d.h"

using namespace std;
class Ekran : public QWidget
{
    Q_OBJECT
public:
    explicit Ekran(QWidget *parent = 0);
    ~Ekran();

    void setWorld(World* world);

    void handleKeyPress(QKeyEvent* event);
    void handleKeyRelease(QKeyEvent* event);

protected:
    void paintEvent(QPaintEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);

private:
    int width, height;
    Renderer3D renderer;

    bool pressed;
    int mx, my;

    World* world;

signals:

public slots:

};

#endif // EKRAN_H
