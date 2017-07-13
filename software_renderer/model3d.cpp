#include "model3d.h"

#include <QDebug>

#include <fstream>
#include <sstream>

void face::add(std::string &element) {
    unsigned v, vt = 0, vn = 0;
    char c;
    std::stringstream el_ss(element);
    el_ss>>v;
    if(el_ss>>c) {
        unsigned value;
        if(el_ss>>value) {
            vt = value;
        }
        if(el_ss>>c){
            if(el_ss>>value) {
                vn = value;
            }
        }
    }
    //qDebug("add: %u %u %u", v, vt, vn);
    vs_id.push_back(v);
    vts_id.push_back(vt);
    vns_id.push_back(vn);
}

model3d::model3d() {
    vs.push_back({0,0,0});
    vts.push_back({0,0});
    vns.push_back({0,0,0});

    // jeśli w pliku nie zdefiniowano materiału, ani nie użyto żadnego
    mtl_index.insert({"", 0});
    mtls.push_back(material());
}

model3d::model3d(std::stringstream& ss) : model3d() {
    std::string line;
    unsigned curr_mtl = 0;
    while(getline(ss, line)){
        std::stringstream line_ss(line);
        std::string type;
        line_ss >> type;
        if(type == "v") {
            float x=0, y=0, z=0;
            line_ss >> x >> y >> z;
            vs.push_back({x, y, z});

        } else if(type == "vt") {
            float u,v;
            line_ss >> u >> v;
            vts.push_back({u, v});

        } else if(type == "vn") {
            float x, y, z;
            line_ss >> x >> y >> z;
            vns.push_back({x, y, z});

        } else if(type == "f") {
            fs.push_back(face(curr_mtl));
            face& f = fs.back();

            std::string token;
            while(line_ss>>token){
                f.add(token);
            }

        } else if(type == "usemtl") {
            std::string mtl_name;
            line_ss >> mtl_name;
            qDebug("usemtl: %s", mtl_name.c_str());
            curr_mtl = get_mtl_index(mtl_name);
        } else if(type == "mtllib") {
            std::string filename;
            line_ss >> filename;
            std::ifstream file(filename);
            if(file) {
                qDebug("mtl file ok: %s", filename.c_str());
                std::stringstream buffer;
                buffer << file.rdbuf();
                load_mtl(buffer);
                file.close();
            }

        } /*else if(type == "#"){
            // komentarz
        }*/
    }

    scx = scy = scz = 1;
    tx = ty = tz = rx = ry = rz = shx = shy = shz = 0;

    qDebug("model prepared");
}

void model3d::load_mtl(std::stringstream& mtl_ss) {
    std::string line;
    unsigned curr_mtl = 0;
    while(getline(mtl_ss, line)){
        std::stringstream line_ss(line);
        std::string type;
        line_ss >> type;
        if(type == "newmtl") {
            std::string mtl_name;
            line_ss >> mtl_name;
            curr_mtl = get_mtl_index(mtl_name);
            qDebug("newmtl: %s (%u)", mtl_name.c_str(), curr_mtl);
        } else if(type == "map_Kd") {
            std::string tex_path;
            line_ss >> tex_path;
            if(mtls[curr_mtl].tex.load_from_file(tex_path.c_str())) {
                qDebug("tex loaded: %s (%u)", tex_path.c_str(), curr_mtl);
            } else {
                qDebug("tex loading failed: %s (%u)", tex_path.c_str(), curr_mtl);
            }
        } else if(type == "Ka") {
            float r, g, b;
            line_ss >> r >> g >> b;
            mtls[curr_mtl].ambient = rgb{r, g, b};
        } else if(type == "Ke") {
            float r, g, b;
            line_ss >> r >> g >> b;
            mtls[curr_mtl].emission = rgb{r, g, b};
        } else if(type == "Kd") {
            float r, g, b;
            line_ss >> r >> g >> b;
            mtls[curr_mtl].diffuse = rgb{r, g, b};
        } else if(type == "Ks") {
            float r, g, b;
            line_ss >> r >> g >> b;
            mtls[curr_mtl].specular = rgb{r, g, b};
        } else if(type == "Ns") {
            float ns;
            line_ss >> ns;
            mtls[curr_mtl].ns = ns;
        } else if(type == "Tr") {
            float tr;
            line_ss >> tr;
            mtls[curr_mtl].d = 1.0 - tr;
        } else if(type == "d") {
            float d;
            line_ss >> d;
            mtls[curr_mtl].d = d;
        } else if(type == "illum") {
            unsigned illum;
            line_ss >> illum;
            mtls[curr_mtl].illum = illum;
        }
    }
}

unsigned model3d::get_mtl_index(std::string mtl_name) {
    if(mtl_index.find(mtl_name) == mtl_index.end()) {
        mtl_index.insert(std::make_pair(mtl_name, mtls.size()));
        mtls.push_back(material());
    }
    return mtl_index.at(mtl_name);
}
