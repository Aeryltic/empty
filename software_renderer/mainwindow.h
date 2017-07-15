#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "world.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(World* world, QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    World* world;
    int model_nr;

    void set_tr_x(float v){ if(world->models.size()) world->models[model_nr - 1].tx = v; update(); }
    void set_tr_y(float v){ if(world->models.size()) world->models[model_nr - 1].ty = v; update(); }
    void set_tr_z(float v){ if(world->models.size()) world->models[model_nr - 1].tz = v; update(); }

    void set_rot_x(float v){ if(world->models.size()) world->models[model_nr - 1].rx = v / 180 * M_PI; update(); }
    void set_rot_y(float v){ if(world->models.size()) world->models[model_nr - 1].ry = v / 180 * M_PI; update(); }
    void set_rot_z(float v){ if(world->models.size()) world->models[model_nr - 1].rz = v / 180 * M_PI; update(); }

    void set_sh_x(float v){ if(world->models.size()) world->models[model_nr - 1].shx = v / 10; update(); }
    void set_sh_y(float v){ if(world->models.size()) world->models[model_nr - 1].shy = v / 10; update(); }
    void set_sh_z(float v){ if(world->models.size()) world->models[model_nr - 1].shz = v / 10; update(); }

    void set_sc_x(float v){ if(world->models.size()) world->models[model_nr - 1].scx = 1 + v/20; update(); }
    void set_sc_y(float v){ if(world->models.size()) world->models[model_nr - 1].scy = 1 + v/20; update(); }
    void set_sc_z(float v){ if(world->models.size()) world->models[model_nr - 1].scz = 1 + v/20; update(); }

    void set_l_a(float v){ if(world->lights.size()) world->lights[0].att_a = 1 + v/20; update(); }
    void set_l_b(float v){ if(world->lights.size()) world->lights[0].att_b = std::pow(1.0 / v, 2); update(); }
    void set_l_c(float v){ if(world->lights.size()) world->lights[0].att_c = std::pow(1.0 / v, 3); update(); }

    void set_model_nr(int n){ this->model_nr = n; update(); }

    void keyPressEvent(QKeyEvent* event);
    void keyReleaseEvent(QKeyEvent* event);
};

#endif // MAINWINDOW_H
