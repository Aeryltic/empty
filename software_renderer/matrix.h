#ifndef MATRIX_H
#define MATRIX_H

#include <armadillo>

#include "camera.h"
#include "model3d.h"

namespace matrix{
arma::mat trans(float tx, float ty, float tz);
arma::mat rotOX(float a);
arma::mat rotOY(float a);
arma::mat rotOZ(float a);
arma::mat scale(float scx, float scy, float scz);
arma::mat shearOX(float shy, float shz);
arma::mat shearOY(float shx, float shz);
arma::mat shearOZ(float shx, float shy);

arma::mat head_camera(const camera& cam);
arma::mat flying_camera(const camera& cam);
arma::mat anchored_camera(const camera& cam);


arma::mat local_transform(model3d& model);

arma::mat clip_matrix(float aspect_ratio, float fov, float near, float far);
}


#endif // MATRIX_H
