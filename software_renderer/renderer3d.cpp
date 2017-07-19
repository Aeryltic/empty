#include "renderer3d.h"

#include "matrix.h"

//---------------------------------------------
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
//---------------------------------------------

Renderer3D::Renderer3D(int width, int height, uint32_t base_color) : width(width), height(height), pixBuf(width, height, base_color) {
    center_x = width / 2.0;
    center_y = height / 2.0;
    aspect = float(width) / height;

    near = 10.0;
    far = 5000.0;

    float angle = 90.0;
    fov = 1 / tan(angle/2.0);

    x_mult = - center_x * (aspect > 1 ? 1 / aspect : 1);
    y_mult = - center_y * (aspect > 1 ? 1 : aspect);
}


struct Comp {
    float* z;
    Comp(float* fz) : z(fz) { }
    bool operator()(unsigned i, unsigned j) {
        return z[i] < z[j];
    }
};

void Renderer3D::render_v2(World* world) {
    //qDebug("rendering");

    //calls = 0;
    this->world = world;

    pixBuf.clear();
    pixBuf.clear_zbuf(far);

    arma::mat viewport_transform = //matrix::clip_matrix(aspect, fov, near, far) *
            //matrix::head_camera(world->cam);
            matrix::anchored_camera(world->cam);

     // obliczenie pozycji świateł
    for(unsigned i=0; i<world->lights.size(); i++) {
        light_source& light = world->lights[i];
        arma::mat transform = viewport_transform;
        if(light.parent) transform *= matrix::local_transform(*light.parent);
        arma::vec pos = transform * arma::vec{{light.real_x}, {light.real_y}, {light.real_z}, {1}};
        world->lights[i].x = pos[0];
        world->lights[i].y = pos[1];
        world->lights[i].z = pos[2];
    }

    unsigned model_count = world->models.size();

    // malarz
    unsigned* model_order = new unsigned[model_count];
    float* model_z = new float[model_count];
    for(unsigned i=0; i < model_count; i++) {
        model_order[i] = i;
        model_z[i] = arma::vec(viewport_transform * arma::vec{{world->models[i].tx}, {world->models[i].ty}, {world->models[i].tz}, {1}})[2];
    }
    std::sort(model_order, model_order + model_count, Comp(model_z));
    delete [] model_z;

    for(unsigned m=0; m<model_count; m++) {

        model3d& model = world->models[model_order[m]];
        arma::mat global_transform = viewport_transform * matrix::local_transform(model);

        unsigned total_vertices = model.vs.size();
        arma::vec* v = new arma::vec[total_vertices];
        float* x = new float[total_vertices];
        float* y = new float[total_vertices];

        // transformacje
        for(unsigned i=1; i<total_vertices; i++) {
            v[i] = global_transform * arma::vec{{model.vs[i].x}, {model.vs[i].y}, {model.vs[i].z}, {1}};
            v[i].shed_row(3);
            // punkt na ekranie
            x[i] = ((v[i][0] / v[i][2])) * x_mult;
            y[i] = ((v[i][1] / v[i][2])) * y_mult;

            x[i] += center_x;
            y[i] += center_y;
        }

        unsigned total_faces = model.fs.size();

        // "malarz"
//        float* fz = new float[total_faces];
//        unsigned* order = new unsigned[total_faces];
//        for(unsigned f=0; f<total_faces; f++) {
//            order[f] = f;
//            face& fc = model.fs[f];
//            unsigned vertex_count = fc.size();
//            fz[f] = 0;
//            for(unsigned i=0; i<vertex_count; i++) {
//                fz[f] += v[fc.vs_id[i]][2];
//            }
//            fz[f] /= vertex_count;
//        }
//        std::sort(order, order + total_faces, Comp(fz));
//        delete [] fz;

        // rysowanie
        for(unsigned f=0; f<total_faces; f++) {
            face& fc = model.fs[f];
            material& mtl = model.mtls[fc.mtl_id];
            unsigned vertex_count = fc.size();

            // v1
            unsigned v0 = fc.vs_id[0];
            for(unsigned i=1; i<vertex_count-1; i++) {
                unsigned v1 = fc.vs_id[i], v2 = fc.vs_id[i+1];
                if(v[v0][2] <= near || v[v1][2] <= near || v[v2][2] <= near) continue;
                arma::vec normal = arma::normalise(arma::cross((v[v0]-v[v1]), (v[v1]-v[v2])));
                if(arma::dot(-v[v0], normal) >= 0) { // back-face culling
                    render_point rpts[3] = {
                        {x[v0] + 0.5, y[v0] + 0.5, v[v0], model.vts[fc.vts_id[0]], model.vns[fc.vns_id[0]]},
                        {x[v1] + 0.5, y[v1] + 0.5, v[v1], model.vts[fc.vts_id[i]], model.vns[fc.vns_id[i]]},
                        {x[v2] + 0.5, y[v2] + 0.5, v[v2], model.vts[fc.vts_id[i+1]], model.vns[fc.vns_id[i+1]]}
                    };
                    draw_textured_triangle(rpts, mtl, normal);
                }
            }

            // v2
//            std::vector<render_point> poly;
//            for(unsigned i=0; i<vertex_count; i++) {
//                unsigned v_id = fc.vs_id[i],
//                         vt_id = fc.vts_id[i],
//                         vn_id = fc.vns_id[i];
//                poly.push_back({x[v_id], y[v_id], v[v_id], model.vts[vt_id], model.vns[vn_id]});
//            }
//            for(int i=0; i<poly.size(); i++) { // do interpolacji
//                poly[i].v[2] = 1.0 / poly[i].v[2];
//                poly[i].vt.u *= poly[i].v[2];
//                poly[i].vt.v *= poly[i].v[2];
//            }
//            //clip_poly(poly);
//            int poly_size = poly.size();
//            //if(poly[0].v[2] <= near) continue;
//            for(int i = 1; i < poly_size - 1; i++) {
//                //if(poly[i].v[2] <= near || poly[i].v[i+1] <= near) continue;
//                arma::vec normal = arma::normalise(arma::cross((poly[0].v - poly[i].v), (poly[i].v - poly[i+1].v)));
//                if(arma::dot(-poly[0].v, normal) >= 0) { // back-face culling
//                    render_point rpts[3] = {poly[0], poly[i], poly[i+1]};
//                    draw_textured_triangle(rpts, mtl, normal);
//                }
//            }
        }

        delete [] v;
        delete [] x;
        delete [] y;
//        delete [] order;
    }
//    delete [] model_order;

    // okręgi oznaczające źródło światła
    for(unsigned i=0; i < world->lights.size(); i++) {
        light_source& light = world->lights[i];
        if(light.z < near || light.z > far) continue;

        float x = ((light.x / light.z)) * x_mult;
        float y = ((light.y / light.z)) * y_mult;

        x += center_x;
        y += center_y;

        float r = 20 * (1 - (light.z - near) / float(far - near));
        draw_ellipse(x, y, r, r, 0, 0xffffffff, 16);
    }

    //qDebug("calls: %u", calls);
}

void Renderer3D::draw_textured_triangle(render_point v[3], const material& mtl, const arma::vec& normal) {

    for(int i=0; i<3; i++) { // do interpolacji
        v[i].v[2] = 1.0 / v[i].v[2];
        v[i].vt.u *= v[i].v[2];
        v[i].vt.v *= v[i].v[2];
    }

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

    arma::vec rv{
        {v[2].v[0] + (v[0].v[0] - v[2].v[0]) * a},
        {v[1].v[1]},
        {v[2].v[2] + (v[0].v[2] - v[2].v[2]) * a}
    };

    vertex_texture rvt =  v[2].vt + (v[0].vt - v[2].vt) * a;

    vertex_normal rvn = v[2].vn + (v[0].vn - v[2].vn) * a; // ?

    render_point r = {rx, v[1].y, rv, rvt, rvn};

    render_point top[3] = {v[0], v[1], r};
    if(top[1].x > top[2].x) std::swap(top[1], top[2]);
    render_point bottom[3] = {v[2], v[1], r};
    if(bottom[1].x > bottom[2].x) std::swap(bottom[1], bottom[2]);

    texture_top_tri(top, mtl, normal);
    texture_bottom_tri(bottom, mtl, normal);

//    bresenham2d(v[0].x, v[0].y, v[1].x, v[1].y, 0xffffffff);
//    bresenham2d(v[0].x, v[0].y, v[2].x, v[2].y, 0xffffffff);
//    bresenham2d(v[2].x, v[2].y, v[1].x, v[1].y, 0xffffffff);
//    bresenham3d(v[0], v[1], mtl, normal);
//    bresenham3d(v[0], v[2], mtl, normal);
//    bresenham3d(v[2], v[1], mtl, normal);
}

void Renderer3D::texture_top_tri(render_point v[3], const material& mtl, const arma::vec& normal) {

    render_point top = v[0];
    render_point left = v[1];
    render_point right = v[2];

    float dy = top.y - left.y;

    float lx = left.x,
          rx = right.x;

    float ldx = (top.x - left.x) / dy,
          rdx = (top.x - right.x) / dy;

    float lz = left.v[2],
          rz = right.v[2];

    float ldz = (top.v[2] - left.v[2]) / dy,
          rdz = (top.v[2] - right.v[2]) / dy;


    float middle_bound = (left.y < 0 ? 0 : left.y);
    float top_bound = (top.y >= height-0.6 ? height-1.6 : top.y);

    if(top_bound < 0 || middle_bound >= height) return;

    float temp = middle_bound - left.y;
    lx += (temp) * ldx;
    rx += (temp) * rdx;
    lz += (temp) * ldz;
    rz += (temp) * rdz;

    for(float y = middle_bound; y <= top_bound; y++) {

        float y_a = (y - left.y) / dy;

        vertex_texture ltv = top.vt * y_a + left.vt * (1.0 - y_a);
        vertex_texture rtv = top.vt * y_a + right.vt * (1.0 - y_a);

        float left = (lx < 0 ? 0 : lx),
              right = (rx >= width-0.6 ? width-1.6 : rx) + 0.5;
        float x_rec = 1.0 / (rx - lx);
        int iy = y + 0.5;
        for(float x = left; x <= right; x++) {
            float x_a = (x - lx) * x_rec;
            float z = 1.0 / (lz * (1.0 - x_a) + rz * x_a);
            int ix = x + 0.5;
            int i = iy * width + ix;
            if((z < pixBuf.z_buf[i]) && (z >= near)) {
                vertex_texture vt = ltv * (1.0 - x_a) +  rtv * x_a;
                float u = vt.u * z;
                float v = vt.v * z;
                if(u < 0)u = 0;
                if(u > 1)u = 1;
                if(v < 0)v = 0;
                if(v > 1)v = 1;
                pixBuf.buf[i] = calc_color(ix, iy, z, u, v, normal, mtl);
                pixBuf.z_buf[i] = z;
            }
        }

        lx += ldx;
        rx += rdx;
        lz += ldz;
        rz += rdz;
    }
}

void Renderer3D::texture_bottom_tri(render_point v[3], const material& mtl, const arma::vec& normal) {

    render_point bottom = v[0];
    render_point left = v[1];
    render_point right = v[2];

    float dy = left.y - bottom.y;

    float lx = left.x,
          rx = right.x;

    float ldx = (bottom.x - left.x) / dy,
          rdx = (bottom.x - right.x) / dy;

    float lz = left.v[2],
          rz = right.v[2];

    float ldz = (bottom.v[2] - left.v[2]) / dy,
          rdz = (bottom.v[2] - right.v[2]) / dy;

    float middle_bound = (left.y >= height-0.6 ? height-1.6 : left.y);
    float bottom_bound = (bottom.y < 0 ? 0 : bottom.y);

    if(bottom_bound >= height-0.6 || middle_bound < 0) return;

    float temp = left.y - middle_bound;
    lx += (temp) * ldx;
    rx += (temp) * rdx;
    lz += (temp) * ldz;
    rz += (temp) * rdz;

    for(float y = middle_bound; y >= bottom_bound; y--) {

        float a = (left.y - y) / dy;

        vertex_texture ltv = left.vt * (1 - a) + bottom.vt * a;
        vertex_texture rtv = right.vt * (1 - a) + bottom.vt * a;

        float left = (lx < 0 ? 0 : lx) + 0.5,
              right = (rx >= width-0.6 ? width-1.6 : rx) + 0.5;
        float x_rec = 1.0 / float(rx - lx);
        int iy = y + 0.5;
        for(float x = left; x <= right; x++) {
            float alpha = (x - lx) * x_rec;
            float z = 1.0 / (lz * (1 - alpha) + rz * alpha);
            int ix = x + 0.5;
            int i = iy * width + ix;
            if((z < pixBuf.z_buf[i]) && (z >= near)) {
                vertex_texture vt = ltv * (1.0 - alpha) +  rtv * alpha;
                float u = vt.u * z;
                float v = vt.v * z;
                if(u < 0)u = 0;
                if(u > 1)u = 1;
                if(v < 0)v = 0;
                if(v > 1)v = 1;
                pixBuf.buf[i] = calc_color(ix, iy, z, u, v, normal, mtl);
                pixBuf.z_buf[i] = z;
            }
        }

        lx += ldx;
        rx += rdx;
        lz += ldz;
        rz += rdz;

    }
}

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

uint32_t treshold(uint32_t v) {
    return v > 255 ? 255 : v;
}

uint32_t Renderer3D::calc_color(float x, float y, float z, float tu, float tv, const arma::vec& n, const material &mtl) {
//    static int gamma = 8;
    if(world->lights.size()) {
        x -= center_x;
        y -= center_y;

        x *= z / x_mult;
        y *= z / y_mult;

        arma::vec v = arma::normalise(arma::vec{{-x}, {-y}, {-z}});

        rgb illumination = mtl.emission; // wynikowa iluminacja

        for(unsigned i = 0; i < world->lights.size(); i++) {
            light_source& light = world->lights[i];

            arma::vec l = arma::vec{{light.x - x}, {light.y - y}, {light.z - z}};
            float dist2 = l[0] * l[0] + l[1] * l[1] + l[2] * l[2];
            float att = 1.0 / (light.att_a + light.att_b * std::sqrt(dist2) + light.att_c * dist2);

            if(att > 1.0 / 255.0) {
                l = arma::normalise(l);
                arma::vec r = arma::normalise(-l - 2 * arma::dot(-l, n) * n);

                rgb ambient = mtl.ambient * light.ambient;
                rgb diffuse(0, 0, 0);
                rgb specular(0, 0, 0);

                float dot_ln = arma::dot(l, n);
                if(dot_ln > 0) {
                    diffuse = mtl.diffuse * dot_ln * light.diffuse;

                    float dot_rv = arma::dot(r, v);
                    if(dot_rv > 0) {
//                        arma::vec c = arma::cross(r, v);
//                        float lambda = arma::dot(c, c) / 2.0;//1 - dot_rv;
//                        float beta = mtl.ns / gamma;
                        specular = mtl.specular
                                //* std::pow((1.0 - beta * lambda), gamma)
                                * std::pow(dot_rv, mtl.ns)
                                * light.specular;
                    }
                }

                illumination += ambient + (diffuse + specular) * att; // ?
            }
        }

        uint32_t color = mtl.tex.get(tu, tv);
        uint32_t lr = treshold(((color >> 16) & 0xff) * illumination.r);
        uint32_t lg = treshold(((color >> 8) & 0xff) * illumination.g);
        uint32_t lb = treshold(((color) & 0xff) * illumination.b);
        return (color & (0xff << 24)) | lr << 16 | lg << 8 | lb;
    }
    return mtl.tex.get(tu, tv); // tryb bez świateł
    //return 0xff000000;
}


void Renderer3D::draw_ellipse(int x, int y, double r1, double r2, double angle, uint32_t color, int s) {
    double cosa = std::cos(angle);
    double sina = std::sin(angle);

    double px = r1 * cosa + x;
    double py = r1 * sina + y;
    for(int i=1; i<=s; i++)
    {
        double a1 = 2.0 * M_PI / s * i;
        double x1 = r1 * cos(a1);
        double y1 = r2 * sin(a1);

        double cx = x1 * cosa - y1 * sina + x;
        double cy = x1 * sina + y1 * cosa + y;

        bresenham2d(cx, cy, px, py, color);
        px = cx;
        py = cy;
    }
}

void Renderer3D::bresenham3d(render_point a, render_point b, const material& mtl, const arma::vec& normal) {

    bool changed = false;
    int x = a.x;
    int y = a.y;

    int dx = abs(b.x - a.x);
    int dy = abs(b.y - a.y);

    int signx = sgn(b.x - a.x);
    int signy = sgn(b.y - a.y);

    if(dy > dx) {
        std::swap(dx, dy);
        changed = true;
    }

    int e = 2 * dy - dx;
    //---------------------------------
    float dz = (b.v[2] - a.v[2]) / dx;
    float z_rec = a.v[2];

    float alpha = 0;
    float dalpha = 1.0 / dx;
    //---------------------------------
    for(int i = 1; i <= dx; i++) {
        int index = y * width + x;
        float z = 1.0 / z_rec;
        if(y >= 0 && y < height && x >= 0 && x < width) {
            if(pixBuf.z_buf[index] > z) {
                float u = (a.vt.u * (1 - alpha) + b.vt.u * alpha) * z,
                      v = (a.vt.v * (1 - alpha) + b.vt.v * alpha) * z;
                if(u > 1) u = 1;
                if(v > 1) v = 1;
                pixBuf.buf[index] = calc_color(x, y, z, u, v, normal, mtl);//0xff00ff00;//
                pixBuf.z_buf[index] = z;
            }
        }
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
        z_rec += dz;
        alpha += dalpha;
        //---------------------------------
    }
}

void Renderer3D::bresenham2d(int x1, int y1, int x2, int y2, uint32_t color) {
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
       // pixBuf.set_pixel(x, y, color);
        pixBuf.buf[y * width + x] = color;
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

char Renderer3D::cohen_sutherland_clipping(int& x0, int& y0, int& x1, int& y1) {
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

render_point intersection(render_point s, render_point e, vec2d a, vec2d b) {
    float d1 = ((b.y-a.y)*(e.x-s.x)-(e.y-s.y)*(b.x-a.x)),
          d2 = ((e.y-s.y)*(b.x-a.x)-(b.y-a.y)*(e.x-s.x));
    //if((abs(d1) > 0.001) && (abs(d2) > 0.001)) {
    float x=((b.x-a.x)*(e.x*s.y-e.y*s.x)-(e.x-s.x)*(b.x*a.y-b.y*a.x))/d1;
    float y=((e.y-s.y)*(b.x*a.y-b.y*a.x)-(b.y-a.y)*(e.x*s.y-e.y*s.x))/d2;
    float alpha;
    if(abs(e.x - s.x) > 0.0001)
        alpha = (x - s.x) / (e.x - s.x);
    else
        alpha = (y - s.y) / (e.y - s.y);
    return s * (1 - alpha) + e * alpha;
    //}
}

void Renderer3D::clip_poly(std::vector<render_point>& poly) {
    if(poly.size() <= 2) return;
    vec2d clipping_rect[4] = {{0,0}, {float(width-1), 0.0}, {float(width-1), float(height-1)}, {0.0, float(height-1)}};
    for(int i = 0; (i < 4) && !poly.empty(); i++) {
        vec2d a = clipping_rect[i], b = clipping_rect[(i+1)%4];
        arma::vec clipping_edge{{b.x - a.x}, {b.y - a.y}, {0}};
        std::vector<render_point> input = poly;
        poly.clear();
        render_point s = input.back();
        for(render_point& e : input) {
            int e_side = arma::vec(arma::cross(clipping_edge, arma::vec{{e.x - a.x}, {e.y - a. y}, {0}}))[2],
                s_side = arma::vec(arma::cross(clipping_edge, arma::vec{{s.x - a.x}, {s.y - a. y}, {0}}))[2];
            if(e_side > 0) {

                if(s_side < 0) {
                    poly.push_back(intersection(s, e, a, b));
                }
                poly.push_back(e);
            } else if(s_side > 0) {
                poly.push_back(intersection(s, e, a, b));
            }
            s = e;
        }
    }
}


/*
void Renderer3D::texture_tri(render_point v[3], const material& mtl, const arma::vec& normal) {

    render_point top = v[0];
    render_point left = v[1];
    render_point right = v[2];

    int dy = top.y - left.y;
    int sign_y = sgn(dy);
    dy *= sign_y;
    if(!sign_y) return bresenham3d(v[1], v[2], mtl, normal); // jeśli trójkąt jest zdegenerowany to tylko pozioma linia


    float lx = left.x,
          rx = right.x;
    float ldx = (top.x - left.x + sgn(top.x - left.x)) / dy,
          rdx = (top.x - right.x + sgn(top.x - right.x)) / dy;

    float lz = left.v[2],
          rz = right.v[2];
    float ldz = (top.v[2] - left.v[2]) / dy,
          rdz = (top.v[2] - right.v[2]) / dy;

    //
    int middle_bound = left.y + 0.5;
    int top_bound = top.y + 0.5;

    if(((sign_y > 0) && (top_bound < 0 || middle_bound >= height))
    || ((sign_y < 0) && (middle_bound < 0 || top_bound >= height))) return; // poza ekranem

    if(middle_bound < 0) middle_bound = 0;
    else if(middle_bound >= height)middle_bound = height - 1;
    if(top_bound < 0) top_bound = 0;
    else if(top_bound >= height)top_bound = height - 1;
    //

    float temp = middle_bound - left.y;
    lx += (temp) * ldx;
    rx += (temp) * rdx;
    lz += (temp) * ldz;
    rz += (temp) * rdz;

    int y = middle_bound;
    dy++;
    for(int i=0; i < dy; i++) {

        float a = std::abs(y - middle_bound) / float(dy);

        vertex_texture ltv = top.vt * a + left.vt * (1.0 - a);
        vertex_texture rtv = top.vt * a + right.vt * (1.0 - a);

        int left = (lx < 0 ? 0 : lx) + 0.5,
            right = (rx >= width-0.6 ? width-1.6 : rx) + 0.5;
        int dx = rx - lx;
        if(dx) {
            float x_rec = 1.0 / float(dx);
            for(int x = left; x <= right; x++) {
                float alpha = float(x - lx) * x_rec;
                float z = 1.0 / (lz * (1.0 - alpha) + rz * alpha);
                int i = y * width + x;
                if((z < pixBuf.z_buf[i]) && (z >= near)) {
                    vertex_texture vt = (ltv * (1.0 - alpha) +  rtv * alpha) * z;
                    float u = vt.u;
                    float v = vt.v;
                    if(u < 0)u = 0;
                    if(u > 1)u = 1;
                    if(v < 0)v = 0;
                    if(v > 1)v = 1;
                    pixBuf.buf[i] = calc_color(x, y, z, u, v, normal, mtl);
                    pixBuf.z_buf[i] = z;
                }
            }
        } else {
            float z = 1 / lz;
            int x = lx + 0.5;
            int i = y * width + x;
            if((z < pixBuf.z_buf[i]) && (z >= near)) {
                vertex_texture vt = ltv * z;
                float u = vt.u;
                float v = vt.v;
                pixBuf.buf[i] = calc_color(x, y, z, u, v, normal, mtl);
                pixBuf.z_buf[i] = z;
            }
        }

        lx += ldx;
        rx += rdx;
        lz += ldz;
        rz += rdz;

        y += sign_y;
    }
}
*/

/*
void Renderer3D::bresenham_tri_tex(render_point v[3], const material& mtl, const arma::vec& normal) {

    bool lchanged = false;
    bool rchanged = false;

    int tx = v[0].x + 0.5;
    int ty = v[0].y + 0.5;
    int lx = v[1].x + 0.5;
    int ly = v[1].y + 0.5;
    int rx = v[2].x + 0.5;
    int ry = v[2].y + 0.5;

    int ldx = abs(lx - tx);
    int ldy = abs(ly - ty);
    int rdx = abs(rx - tx);
    int rdy = abs(ry - ty);

    int lsignx = sgn(tx - lx);
    int rsignx = sgn(tx - rx);
    int signy = sgn(ty - ly);

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
    float ldz = float(v[1].v[2] - v[0].v[2]) / ldx;
    float lz = v[0].v[2];
    float rdz = float(v[2].v[2] - v[0].v[2]) / rdx;
    float rz = v[0].v[2];

//    float lu = v[1].vt.u,
//          lv = v[1].vt.v,
//          ru = v[2].vt.u,
//          rv = v[2].vt.v;
    //---------------------------------
    int li=0, ri=0;
    //int y = ty;
    int y = ly;
    int w = ly;
//    qDebug("from: %d to: %d", ly, ty);
    while(signy < 0 ? w >= ty : w <= ty) {

        while((y==ly)&&(li<ldx)){
            bresenham_next_pixel(lx, ly, ldx, ldy, lsignx, signy, le, lchanged);
            lz += ldz;
            li++;
        }
        while((y==ry)&&(ri<rdx)) {
            bresenham_next_pixel(rx, ry, rdx, rdy, rsignx, signy, re, rchanged);
            ri++;
            rz += rdz;
        }
        // scanline
        float a = 0,
              da = 1.0 / float(rx - lx + 1);
        for(int x = lx; x <= rx; x++) {
            float z = 1.0/(lz * (1 - a) + rz * a);
            //float u =
            pixBuf.buf[y*width + x] = 0xffffffff;
            a += da;
        }
        y = ly;
        signy < 0 ? w-- : w++;
//        qDebug("y: %d, w: %d", y, w);
    }
}
*/
