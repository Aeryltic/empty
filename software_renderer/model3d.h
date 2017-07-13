#ifndef MODEL3D_H
#define MODEL3D_H

#include <vector>
#include <string>
#include <unordered_map>

#include "vertex.h"
#include "texture.h"


struct face {
    unsigned mtl_id;
    std::vector<unsigned> vs_id;
    std::vector<unsigned> vns_id;
    std::vector<unsigned> vts_id;


    face(unsigned mtl) { this->mtl_id = mtl; }

    void add(std::string& element);
    unsigned size() {
        return vs_id.size();
    }
};



struct model3d {
    std::vector<vertex> vs;
    std::vector<vertex_texture> vts;
    std::vector<vertex_normal> vns;
    std::vector<face> fs;

    std::unordered_map<std::string, unsigned> mtl_index;
    std::vector<material> mtls;

    float tx, ty, tz,
          rx, ry, rz,
          shx, shy, shz,
          scx, scy, scz;


    model3d(std::stringstream &ss); // pozbyć się tego i zrobić jakiś load_from_file
    virtual ~model3d() {}

    void load_mtl(std::stringstream& mtl_ss);
    unsigned get_mtl_index(std::string mtl_name);
    material& get_mtl(std::string mtl_name);

    face& operator[](unsigned i) {
        return fs.at(i);
    }
//    vertex& getv(unsigned f, unsigned v) {
//        return vs.at(fs.at(f)[v]);
//    }

private:
    model3d();
};


#endif // MODEL3D_H
