#include "matrix.h"

namespace matrix {

arma::mat trans(float tx, float ty, float tz) {
    return arma::mat{
        {1, 0, 0, tx},
        {0, 1, 0, ty},
        {0, 0, 1, tz},
        {0, 0, 0, 1}
    };
}

arma::mat rotOX(float a) {
    return arma::mat{
        {1, 0, 0, 0},
        {0, cos(a), -sin(a), 0},
        {0, sin(a), cos(a), 0},
        {0, 0, 0, 1}
    };
}

arma::mat rotOY(float a) {
    return arma::mat{
        {cos(a), 0, sin(a), 0},
        {0, 1, 0, 0},
        {-sin(a), 0, cos(a), 0},
        {0, 0, 0, 1}
    };
}

arma::mat rotOZ(float a) {
    return arma::mat{
        {cos(a), -sin(a), 0, 0},
        {sin(a), cos(a), 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
}

arma::mat scale(float scx, float scy, float scz)
{
    return arma::mat{
        {scx, 0, 0, 0},
        {0, scy, 0, 0},
        {0, 0, scz, 0},
        {0, 0, 0, 1}
    };
}

arma::mat shearOX(float shy, float shz)
{
    return arma::mat{
        {1, 0, 0, 0},
        {shy, 1, 0, 0},
        {shz, 0, 1, 0},
        {0, 0, 0, 1}
    };
}

arma::mat shearOY(float shx, float shz){
    return arma::mat{
        {1, shx, 0, 0},
        {0, 1, 0, 0},
        {0, shz, 1, 0},
        {0, 0, 0, 1}
    };
}

arma::mat shearOZ(float shx, float shy){
    return arma::mat{
        {1, 0, shx, 0},
        {0, 1, shy, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
}

arma::mat head_camera(const camera& cam){
    return rotOZ(-cam.rz) * rotOX(-cam.rx) * rotOY(-cam.ry) * trans(-cam.px, -cam.py, -cam.pz);
}

arma::mat flying_camera(const camera& cam){
    return rotOY(-cam.ry) * rotOX(-cam.rx) * rotOZ(-cam.rz) * trans(-cam.px, -cam.py, -cam.pz);
}

arma::mat local_transform(model3d& model) {
    return    trans(model.tx, model.ty, model.tz)
            * rotOY(model.ry) * rotOX(model.rx) * rotOZ(model.rz)
            * shearOX(model.shy, model.shz) * shearOY(model.shx, model.shz) * shearOZ(model.shx, model.shy)
            * scale(model.scx, model.scy, model.scz);
}

arma::mat clip_matrix(float aspect_ratio, float fov, float near, float far)
{
    return arma::mat{
        {fov*aspect_ratio, 0, 0, 0},
        {0, fov, 0, 0},
        {0, 0, (far+near)/(far-near), 1},
        {0, 0, (2*near*far)/(near-far), 0}
    };
}

arma::mat anchored_camera(const camera &cam)
{
    return trans(-cam.px, -cam.py, -cam.pz) * rotOZ(-cam.rz) * rotOX(-cam.rx) * rotOY(-cam.ry);
}


}
