#include "pixbuf.h"

#include <algorithm>
//#include <limits>

#include "matrix.h"

//---------------------------------------------
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
//---------------------------------------------

PixBuf::PixBuf(int width, int height, uint32_t color) {
    this->width = width;
    this->height = height;

    buf = new uint32_t[height * width];
    z_buf = new float[height * width];
    this->baseColor = color;
}

PixBuf::~PixBuf() {
    delete [] buf;
    delete [] z_buf;
}

bool PixBuf::resize(int new_width, int new_height) {
    width = new_width;
    height = new_height;
    delete [] buf;
    buf = new uint32_t[height * width];
    return true;
}

void PixBuf::clear() {
    for(int y = 0; y < height; y++)
        for(int x = 0; x < width; x++)
            buf[y*width+x] = baseColor;
}

void PixBuf::clear_zbuf(float far) {
    //float value = std::numeric_limits<float>::max();
    for(int i=0; i<height*width; i++) {
        z_buf[i] = far;
    }
}

//uint32_t PixBuf::unpackColor(const QColor &color) {
//    int c = (color.alpha() << 24) | (color.red() << 16) | (color.green() << 8) | color.blue();
//    return c;
//}

//void PixBuf::copy_buffer(const PixBuf& other) {
//    int tw = (other.width < width ? other.width : width);
//    int th = (other.buf_h < buf_h ? other.buf_h : buf_h);
//    clear();
//    for(int y = 0; y < th; y++)
//        for(int x = 0; x < tw; x++)
//            buf[y*width+x] = other.buf[y*other.width+x];
//}

//float dist(arma::vec v) {
//    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
//}

//float PixBuf::calc_zb(float z) {
//    return ((far + near) / 2 * (far - near) + (-far * near)/(far - near)/z + 1.5)/2;
//}
/*
void PixBuf::render(model_base &models, camera viewport) {
    //qDebug("rendering");
    clear();
    clear_zbuf();

    arma::mat viewport_transform = //matrix::clip_matrix(aspect, fov, near, far) *
            //matrix::head_camera(viewport);
            matrix::anchored_camera(viewport);

    for(unsigned i=0; i<models.lights.size(); i++) { // obliczenie pozycji świateł
        arma::vec pos = viewport_transform * arma::vec{{models.lights[i].real_x}, {models.lights[i].real_y}, {models.lights[i].real_z}, {1}};
        models.lights[i].x = pos[0];
        models.lights[i].y = pos[1];
        models.lights[i].z = pos[2];
    }

    // właściwy render
    unsigned model_count = models.size();
    for(unsigned m=0; m<model_count; m++) {
        model3d& model = models[m];
        arma::mat global_transform = viewport_transform * matrix::local_transform(model);

        unsigned total_vertices = model.vs.size();
        arma::vec* v = new arma::vec[total_vertices];
        float* x = new float[total_vertices];
        float* y = new float[total_vertices];
        float* z = new float[total_vertices];

        for(unsigned i=1; i<total_vertices; i++) {
            v[i] = global_transform * arma::vec{{model.vs[i].x}, {model.vs[i].y}, {model.vs[i].z}, {1}};
            v[i].shed_row(3); // do back-face cullingu
            // punkt na ekranie
            x[i] = (-(v[i][0] / v[i][2])) * center_x;
            y[i] = (-(v[i][1] / v[i][2])) * center_y;
            if(aspect > 1.0) {
                x[i] = x[i] / aspect;
            } else {
                y[i] = y[i] * aspect;
            }
            x[i] += center_x;
            y[i] += center_y;

            z[i] = v[i][2];//dist(v[i]);// calc_zb(v[i][2]);//;
        }

        unsigned total_faces = model.fs.size();
        for(unsigned f=0; f<total_faces; f++) {
            face& fc = model.fs[f];
            material& mtl = model.mtls[fc.mtl_id];
            unsigned vertex_count = fc.size();
            unsigned v0 = fc[0];
            for(unsigned i=1; i<vertex_count-1; i++) {
                unsigned v1 = fc[i], v2 = fc[i+1];
                if(v[v0][2] <= near || v[v1][2] <= near || v[v2][2] <= near) continue;
                arma::vec normal = arma::normalise(arma::cross((v[v1]-v[v0]), (v[v2]-v[v0])));
                if(arma::dot(-v[v0], normal) >= 0) { // back-face culling
                    render_point rpts[3] = {
                        {x[v0], y[v0], z[v0], model.vts[fc.vts_id[0]], model.vns[fc.vns_id[0]]},
                        {x[v1], y[v1], z[v1], model.vts[fc.vts_id[i]], model.vns[fc.vns_id[i]]},
                        {x[v2], y[v2], z[v2], model.vts[fc.vts_id[i+1]], model.vns[fc.vns_id[i+1]]}
                    };
                    draw_textured_triangle(rpts, mtl, normal);
                }
            }
        }

        delete [] v;
        delete [] x;
        delete [] y;
        delete [] z;
    }
}
*/
/*
void PixBuf::draw_textured_triangle(render_point v[3], const material& mtl, const arma::vec& normal) {
    if(v[0].y < v[1].y){
        std::swap(v[0], v[1]);
    }
    if(v[0].y < v[2].y){
        std::swap(v[0], v[2]);
    }
    if(v[1].y < v[2].y){
        std::swap(v[1], v[2]);
    }

    float a = (v[1].y - v[2].y) / (v[0].y - v[2].y);

    float rx = v[2].x + (v[0].x - v[2].x) * a;
    float rz = v[2].z + (v[0].z - v[2].z) * a;

    float ru = v[2].vt.u + (v[0].vt.u - v[2].vt.u) * a;
    float rv = v[2].vt.v + (v[0].vt.v - v[2].vt.v) * a;

    float rnx = v[2].vn.x + (v[0].vn.x - v[2].vn.x) * a;
    float rny = v[2].vn.y + (v[0].vn.y - v[2].vn.y) * a;
    float rnz = v[2].vn.z + (v[0].vn.z - v[2].vn.z) * a;

    render_point r = {rx, v[1].y, rz, {ru, rv}, {rnx, rny, rnz}};

    render_point top[3] = {v[0], v[1], r};
    if(top[1].x > top[2].x) std::swap(top[1], top[2]);
    render_point bottom[3] = {v[2], v[1], r};
    if(bottom[1].x > bottom[2].x) std::swap(bottom[1], bottom[2]);

    texture_top_tri(top, mtl, normal);
    texture_bottom_tri(bottom, mtl, normal);

//    bresenham2d(v[0].x, v[0].y, v[1].x, v[1].y, 0xffffffff);
//    bresenham2d(v[0].x, v[0].y, v[2].x, v[2].y, 0xffffffff);
//    bresenham2d(v[2].x, v[2].y, v[1].x, v[1].y, 0xffffffff);
}

void PixBuf::texture_top_tri(render_point v[3], const material& mtl, const arma::vec& normal) {

    render_point top = v[0];
    render_point left = v[1];
    render_point right = v[2];

    float dy = top.y - left.y;

    float lx = left.x,
          rx = right.x;

    float ldx = (top.x - left.x) / dy,
          rdx = (top.x - right.x) / dy;

    float lz = left.z,
          rz = right.z;

    float ldz = (top.z - left.z) / dy,
          rdz = (top.z - right.z) / dy;

    int middle_bound = left.y;

    float top_zr = 1.0f / top.z;
    float left_zr = 1.0f / left.z;
    float right_zr = 1.0f / right.z;

    for(int y = middle_bound; y <= top.y; y++) {

        float a = (y - middle_bound) / dy;
        //if(a<0 || a>1)continue; // tylko testy

        // affine texture mapping (works okey-ish)
        texture_vertex ltv = top.vt * a + left.vt * (1.0 - a);
        texture_vertex rtv = top.vt * a + right.vt * (1.0 - a);

        // perspective correct mapping (doesn't work... or does?)
        //float one_m_a = 1.0 - a;
//        texture_vertex ltv = (top.vt * top_zr * (a) + left.vt * left_zr * (1.0 - a))
//                            /(top_zr * (a) + left_zr * (1.0 - a));
//        texture_vertex rtv = (top.vt * (1.0 - a) * top_zr + right.vt * right_zr * a)
//                            /(top_zr * (a) + right_zr * (1.0 - a));

        int left = lx - 0.5,
            right = rx + 0.5;
        float reciprocal = 1.0 / float(right - left);
        for(int x = left; x <= right; x++) {
            float alpha = (x - left) * reciprocal;
            float z = lz * (1.0 - alpha) + rz * alpha;
            if(is_inside(x, y)) {
                int i = y * width + x;
                if((z < z_buf[i]) && (z >= near) && (z <= far)) {
                    texture_vertex vt = ltv * (1.0 - alpha) +  rtv * alpha;
                    if(vt.u < 0)vt.u = 0;
                    if(vt.u > 1)vt.u = 1;
                    if(vt.v < 0)vt.v = 0;
                    if(vt.v > 1)vt.v = 1;
                    uint32_t color = calc_color(x, y, z, vt.u, vt.v, normal, arma::vec{ mtl);//mtl.tex.get(vt.u, vt.v);// 0xffffffff;//
                    buf[i] = color;
                    z_buf[i] = z;
                    //set_pixel(x, y, z, color);
                }
            }
        }

        lx += ldx;
        rx += rdx;
        lz += ldz;
        rz += rdz;
    }
}

void PixBuf::texture_bottom_tri(render_point v[3], const material& mtl, const arma::vec& normal) {

    render_point bottom = v[0];
    render_point left = v[1];
    render_point right = v[2];

    float dy = left.y + 0.5 - bottom.y;

    float lx = left.x,
          rx = right.x;

    float ldx = (bottom.x - left.x) / dy,
          rdx = (bottom.x - right.x) / dy;

    float lz = left.z,
          rz = right.z;

    float ldz = (bottom.z - left.z) / dy,
          rdz = (bottom.z - right.z) / dy;

    int middle_bound = left.y + 0.5;
    float zr[3] = {1.0f / bottom.z, 1.0f / left.z, 1.0f / right.z};

    for(int y = middle_bound; y >= bottom.y; y--) {

        float a = (middle_bound - y) / dy;
       // if(a>1 || a<0)continue;

        texture_vertex ltv = left.vt * (1 - a) + bottom.vt * a;
        texture_vertex rtv = right.vt * (1 - a) + bottom.vt * a;

//        texture_vertex ltv = (left.vt * (1 - a) * zr[1] + bottom.vt * (a) * zr[0])
//                            /((1 - a) * zr[1] + (a) * zr[0]);
//        texture_vertex rtv = (right.vt * (1 - a) * zr[2] + bottom.vt * (a) * zr[0])
//                            /((1 - a) * zr[2] + (a) * zr[0]);

        int left = lx - 0.5,
            right = rx + 0.5;
        float reciprocal = 1.0 / (right - left);
        for(int x = left; x <= right; x++) {
            float alpha = (x - left) * reciprocal;
            float z = lz * (1 - alpha) + rz * alpha;
            texture_vertex vt = ltv * (1 - alpha) +  rtv * alpha;
            if(vt.u < 0)vt.u = 0;
            if(vt.u > 1)vt.u = 1;
            if(vt.v < 0)vt.v = 0;
            if(vt.v > 1)vt.v = 1;
            uint32_t color = mtl.tex.get(vt.u, vt.v);// 0xffffffff;//
//            uchar c = 0xff * (1.0 - z / far);
//            uint32_t color = 0xff << 24 | c << 16 | c << 8 | c;
            set_pixel(x, y, z, color);
        }

        lx += ldx;
        rx += rdx;
        lz += ldz;
        rz += rdz;

    }
}

uint32_t treshold(uint32_t v) {
    return v > 255 ? 255 : v;
}

uint32_t PixBuf::calc_color(float x, float y, float z, float tu, float tv, const arma::vec& normal, const arma::vec& view, const material &mtl, std::vector<light_source> &lights) {
    // wektory muszą(?) być znormalizowane

//        L_a - ambient,
//        L_d - diffuse,
//        L_s - specular
//        a*x + b*y + c*z + d = 0
//        normal = (a, b, c)
//    // AMBIENT
//        k_a - surface’s coefficient of ambient reflection
//        I_a - ambient component of illumination
//    #   I_a = k_a * L_a;
//    // DIFFUSE
//        k_d - surface’s coefficient of diffuse reflection
//    #   I_d = k_d * max(0, dot(n, l) * L_d
//    // SPECULAR
//        k_s - surface’s coefficient of specular reflection
//        s - shininess (Ns ?)
//        r = 2 * (dot(normal, light_v)) * normal - light_v
//    #1  I_s = k_s * dot(r, view_v)^s * L_s
//        h = normalize(light_v + view_v)
//    #2  I_s = k_s * dot(n, h)^s * L_s
//    // TOTAL
//        I = I_e + I_a + (I_d + I_s) / (a + b*d + c*d*d)

    if(lights.size()) { // na razie max 1 światło
        light_source& light = lights[0];
        // light and reflection vectors
        arma::vec l{{light.x - x}, {light.y - y}, {light.z - z}};
        arma::vec r = 2 * arma::dot(normal, l) * normal - l;
        // obliczenia światła...
        rgb ambient = mtl.ambient * light.ambient;
        rgb diffuse = mtl.diffuse * std::max(0.0, arma::dot(normal, l)) * light.diffuse;
        rgb specular = mtl.specular * std::pow(arma::dot(r, view), mtl.ns) * light.specular;

        float a = normal[0], b = normal[1], c = normal[2];
        float d = std::sqrt(l[0] * l[0] + l[1] * l[1] + l[2]* l[2]);

        rgb illumination = ambient + (diffuse + specular) / (a + b * d + c * d * d); // czy przy takich liczbach coś się jeszcze ostanie z diffuse i specular?
        // kolory
        uint32_t color = mtl.tex.get(tu, tv); // zamiast tu i tv brać po prostu kolor?
        uint32_t lr = (color >> 16) & 0xff;
        uint32_t lg = (color >> 8) & 0xff;
        uint32_t lb = (color) & 0xff;
        lr = treshold(lr * illumination.r);
        lg = treshold(lg * illumination.g);
        lb = treshold(lb * illumination.b);
        return (color & (0xff << 24)) | lr << 16 | lg << 8 | lb;
    }
    return 0xff000000;
}
*/
/*
void PixBuf::bresenham_tri(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3) {
    float tx = x1, ty = y1, tz = z1,
          mx = x2, my = y2, mz = z2,
          bx = x3, by = y3, bz = z3;

    if(ty < my){
        std::swap(tx, mx);
        std::swap(ty, my);
        std::swap(tz, mz);
    }
    if(ty < by){
        std::swap(tx, bx);
        std::swap(ty, by);
        std::swap(tz, bz);
    }
    if(my < by){
        std::swap(mx, bx);
        std::swap(my, by);
        std::swap(mz, bz);
    }
    int rx   = bx + (tx - bx) / (ty - by) * (my - by);
    float rz = bz + (tz - bz) / (ty - by) * (my - by);

    fill_top_tri(tx+0.5, ty+0.5, mx+0.5, rx+0.5, my+0.5, tz, mz, rz);
    fill_bottom_tri(bx+0.5, by+0.5, mx+0.5, rx+0.5, my+0.5, bz, mz, rz);

//    bresenham3d(x1, y1, z1, x2, y2, z2);
//    bresenham3d(x1, y1, z1, x3, y3, z3);
//    bresenham3d(x3, y3, z3, x2, y2, z2);

//    bresenham2d(x1, y1, x2, y2, 0xffffffff);
//    bresenham2d(x1, y1, x3, y3, 0xffffffff);
//    bresenham2d(x3, y3, x2, y2, 0xffffffff);
}

void PixBuf::bresenham2d(int x1, int y1, int x2, int y2, uint32_t color) {
    if(!cohen_sutherland_clipping(x1, y1, x2, y2)) return;
    bool changed = false;
    int x = x1;
    int y = y1;

    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);

    int signx = sgn(x2 - x1);
    int signy = sgn(y2 - y1);

    if(dy > dx) {
        std::swap(dx, dy);
        changed = true;
    }

    int e = 2 * dy - dx;
    for(int i = 1; i <= dx; i++) {
        set_pixel(x, y, color);
        while (e >= 0) {
            if(changed) {
                x += signx;
            } else {
                y += signy;
            }
            e -= 2 * dx;
        }
        if(changed) {
            y += signy;
        } else {
            x += signx;
        }
        e += 2 * dy;
    }
}

void PixBuf::bresenham3d(int x1, int y1, float z1, int x2, int y2, float z2, uint32_t color) {

    bool changed = false;
    int x = x1;
    int y = y1;

    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);

    int signx = sgn(x2 - x1);
    int signy = sgn(y2 - y1);

    if(dy > dx) {
        std::swap(dx, dy);
        changed = true;
    }

    int e = 2 * dy - dx;
    //---------------------------------
    float dz = (z2 - z1) / dx;
    float z = z1;
    //---------------------------------
    for(int i = 1; i <= dx; i++) {
        set_pixel(x, y, z, color);
        while (e >= 0) {
            if(changed) {
                x += signx;
            } else {
                y += signy;
            }
            e -= 2 * dx;
        }
        if(changed) {
            y += signy;
        } else {
            x += signx;
        }
        e += 2 * dy;
        //---------------------------------
        z += dz;
        //---------------------------------
    }
}
*/
inline void bresenham_next_pixel(int& x, int& y, int& dx, int& dy, int& signx, int& signy, int& e, bool changed) {
    while (e >= 0) {
        if(changed) {
            x += signx;
        } else {
            y += signy;
        }
        e -= 2 * dx;
    }
    if(changed) {
        y += signy;
    } else {
        x += signx;
    }
    e += 2 * dy;
}
/*
void PixBuf::bresenham3d(int x1, int y1, float z1, int x2, int y2, float z2) {

    bool changed = false;
    int x = x1;
    int y = y1;

    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);

    int signx = sgn(x2 - x1);
    int signy = sgn(y2 - y1);

    if(dy > dx) {
        std::swap(dx, dy);
        changed = true;
    }

    int e = 2 * dy - dx;
    //---------------------------------
    float dz = (z2 - z1) / dx;
    float z = z1;
    //---------------------------------
    for(int i = 1; i <= dx; i++) {
        //---------------------------------
        uchar c = 0xff * (1.0 - z / far);
        //if(c<0)c=0; if(c>255)c=255;
        uint32_t color = 0xff << 24 | c << 16 | c << 8 | c;
        set_pixel(x, y, z, color);
        bresenham_next_pixel(x, y, dx, dy, signx, signy, e, changed);
        //---------------------------------
        z += dz;
        //---------------------------------
    }
}


void PixBuf::fill_top_tri(int tx, int ty, int ax, int bx, int my, float tz, float az, float bz) {

    bool lchanged = false;
    bool rchanged = false;

    int lx = tx;
    int ly = ty;
    int rx = tx;
    int ry = ty;

    int ldx = abs(ax - tx);
    int ldy = abs(my - ty);
    int rdx = abs(bx - tx);
    int rdy = abs(my - ty);

    int lsignx = sgn(ax - tx);
    int rsignx = sgn(bx - tx);
    int signy = sgn(my - ty);

    if(ldy > ldx) {
        std::swap(ldx, ldy);
        lchanged = true;
    }
    if(rdy > rdx) {
        std::swap(rdx, rdy);
        rchanged = true;
    }

    int le = 2 * ldy - ldx;
    int re = 2 * rdy - rdx;
    //---------------------------------
    float ldz = float(az - tz) / ldx;
    float lz = tz;
    float rdz = float(bz - tz) / rdx;
    float rz = tz;
    //---------------------------------
    int li=0, ri=0;
    //int y = ty;
    int w = ly;
    while(w > my) {
        int y = ly;
        while((y==ly)&&(li<ldx)){
            uchar c = 0xff * (1.0 - lz / far);
            uint32_t color = 0xff << 24 | c << 16 | c << 8 | c;
            set_pixel(lx, ly, lz, color);
            bresenham_next_pixel(lx, ly, ldx, ldy, lsignx, signy, le, lchanged);
            //qDebug("l: %d, %d", lx, ly);
            lz += ldz;
            li++;
        }
        y = ry;
        while((y==ry)&&(ri<rdx)) {
            uchar c = 0xff * (1.0 - rz / far);
            uint32_t color = 0xff << 24 | c << 16 | c << 8 | c;
            set_pixel(rx, ry, rz, color);
            bresenham_next_pixel(rx, ry, rdx, rdy, rsignx, signy, re, rchanged);
            //qDebug("r: %d, %d", lx, ly);
            ri++;
            rz += rdz;
        }
        bresenham3d(lx, ly, lz, rx, ry, rz);
        w--;
    }
}


void PixBuf::fill_bottom_tri(int btx, int bty, int ax, int bx, int my, float btz, float az, float bz) {

    bool lchanged = false;
    bool rchanged = false;

    int lx = btx;
    int ly = bty;
    int rx = btx;
    int ry = bty;

    int ldx = abs(ax - btx);
    int ldy = abs(my - bty);
    int rdx = abs(bx - btx);
    int rdy = abs(my - bty);

    int lsignx = sgn(ax - btx);
    int rsignx = sgn(bx - btx);
    int signy = sgn(my - bty);

    if(ldy > ldx) {
        std::swap(ldx, ldy);
        lchanged = true;
    }
    if(rdy > rdx) {
        std::swap(rdx, rdy);
        rchanged = true;
    }

    int le = 2 * ldy - ldx;
    int re = 2 * rdy - rdx;
    //---------------------------------
    float ldz = (az - btz) / ldx;
    float lz = btz;
    float rdz = (bz - btz) / rdx;
    float rz = btz;
    //---------------------------------
    int w = bty;
    int li=0, ri=0;
    while(w <= my) {
        int y = ly;
        while((y==ly)&&(li<ldx)){
            uchar c = 0xff * (1.0 - lz / far);
            uint32_t color = 0xff << 24 | c << 16 | c << 8 | c;
            set_pixel(lx, ly, lz, color);
            bresenham_next_pixel(lx, ly, ldx, ldy, lsignx, signy, le, lchanged);
            lz += ldz;

            li++;
        }
        y = ry;
        while((y==ry)&&(ri<rdx)) {
            uchar c = 0xff * (1.0 - rz / far);
            uint32_t color = 0xff << 24 | c << 16 | c << 8 | c;
            set_pixel(rx, ry, rz, color);
            bresenham_next_pixel(rx, ry, rdx, rdy, rsignx, signy, re, rchanged);
            rz += rdz;

            ri++;
        }
        bresenham3d(lx, ly, lz, rx, ry, rz); // draw_horizontal
        w++;
    }
}


char PixBuf::cohen_sutherland_clipping(int& x0, int& y0, int& x1, int& y1) {
    char b0 =   (y0 < 0 ? 4 : (y0 >= height ? 8 : 0)) |
                (x0 < 0 ? 1 : (x0 >= width ? 2 : 0));
    char b1 =   (y1 < 0 ? 4 : (y1 >= height ? 8 : 0)) |
                (x1 < 0 ? 1 : (x1 >= width ? 2 : 0));
    if((!b0) && (!b1)) return -1; // true
    if(b0 & b1) return 0; // false
    //if(x0 == x1 && y0 == y1) return 0;
    if(b0) {
        if(b0 & 1) { //x < 0
            if(x1 - x0)
                y0 += (y1 - y0) * (-x0) / (x1 - x0);
            x0 = 0;
        } else if(b0 & 2) { // x >= buf_w
            if(x1 - x0)
                y0 += (y1 - y0) * (width - 1 - x0) / (x1 - x0);
            x0 = width - 1;
        }
        if(b0 & 4) {
            if(y1 - y0)
                x0 += (x1 - x0) * (-y0) / (y1 - y0);
            y0 = 0;
        } else if(b0 & 8) {
            if(y1 - y0)
                x0 += (x1 - x0) * (height - 1 - y0) / (y1 - y0);
            y0 = height - 1;
        }
    }
    if(b1) {
        if(b1 & 1) { //x < 0
            if(x0 - x1)
                y1 += (y0 - y1) * (-x1) / (x0 - x1);
            x1 = 0;
        } else if(b1 & 2) { // x >= buf_w
            if(x0 - x1)
                y1 += (y0 - y1) * (width - 1 - x1) / (x0 - x1);
            x1 = width - 1;
        }
        if(b1 & 4) {
            if(y0 - y1)
                x1 += (x0 - x1) * (-y1) / (y0 - y1);
            y1 = 0;
        } else if(b1 & 8) {
            if(y0 - y1)
                x1 += (x0 - x1) * (height - 1 - y1) / (y0 - y1);
            y1 = height - 1;
        }
    }
    return -1;
}


void PixBuf::draw_horizontal(int x0, int x1, int y, uint32_t color)
{
    if(y>=0 && y<height) {
        if(x0 < 0) x0 = 0;
        if(x0 >= width) x0 = width - 1;
        if(x1 < 0) x1 = 0;
        if(x1 >= width) x1 = width - 1;
        // to można zoptymalizować jeszcze, bo to są kolejne komórki pamięci
        for(int x=x0; x<=x1; x++) {
            buf[y*width+x] = color;
        }
    }
}


void PixBuf::draw_vertical(int x, int y0, int y1, uint32_t color)
{
    if(x>=0 && x<width) {
        if(y0 < 0) y0 = 0;
        if(y0 >= height) y0 = height - 1;
        if(y1 < 0) y1 = 0;
        if(y1 >= height) y1 = height - 1;

        for(int y=y0; y<=y1; y++) {
            buf[y*width+x] = color;
        }
    }
}


void PixBuf::draw_diagonal(int x0, int y0, int x1, int y1, uint32_t color)
{
    int d = y1 > y0 ? 1 : -1;
    for(int x=x0, y = y0; x<=x1; x++, y+=d){
        buf[y*width+x] = color;
    }
}
*/
/*
bool PixBuf::load_image(const QString& filepath, const char *format, bool force_resize) {
    QImage img;
    if(!img.load(filepath, format)) return false;

    const uchar* bits = img.constBits();
    int img_h = img.height(),
        img_w = img.width();

    if(force_resize)
        resize(img_w, img_h);
    clear();

    for(int y = 0; y < img_h; y++)
    {
        for(int x = 0; x < img_w; x++)
        {
            int i = (y * img_w + x) * 4;
            uint32_t a = bits[i+3],
                     r = bits[i+2],
                     g = bits[i+1],
                     b = bits[i];
            uint32_t color = a << 24 | r << 16 | g << 8 | b;
            set_pixel_nc(x, y, color);
        }
    }
    return true;
}
*/

/*
void PixBuf::draw_grid(int horizontal, int vertical, float grid_size, const arma::mat& viewport_transform) {

//    for(int i=0; i<horizontal; i++) {
//        arma::vec begin = viewport_transform * arma::vec{{(i-horizontal/2)*grid_size},{0}, {-vertical/2*grid_size}, {1}};
//        arma::vec end = viewport_transform * arma::vec{{(i-horizontal/2)*grid_size},{0}, {vertical/2*grid_size}, {1}};
//        float bx = ((begin[0] / begin[2]) + 1) * center_x;
//        float by = (-(begin[1] / begin[2]) + 1) * center_y;
//        float ex = ((end[0] / end[2]) + 1) * center_x;
//        float ey = (-(end[1] / end[2]) + 1) * center_y;
//        float bz = dist(begin);
//        float ez = dist(end);
//        bresenham3d(bx, by, bz, ex, ey, ez);
//    }

//    for(int i=0; i<vertical; i++){
//        arma::vec begin = viewport_transform * arma::vec{{-horizontal/2*grid_size},{0}, {(i-vertical/2)*grid_size}, {1}};
//        arma::vec end = viewport_transform * arma::vec{{horizontal/2*grid_size},{0}, {(i+vertical/2)*grid_size}, {1}};
//        float bx = ((begin[0] / begin[2]) + 1) * center_x;
//        float by = (-(begin[1] / begin[2]) + 1) * center_y;
//        float ex = ((end[0] / end[2]) + 1) * center_x;
//        float ey = (-(end[1] / end[2]) + 1) * center_y;
//        float bz = dist(begin);
//        float ez = dist(end);
//        bresenham3d(bx, by, bz, ex, ey, ez);
//    }

}
*/
