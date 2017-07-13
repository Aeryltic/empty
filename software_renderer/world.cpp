#include "world.h"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>

World::World() {
    cam.px = 0;
    cam.py = 100;
    cam.pz = -250;

    cam.rx = cam.ry = cam.rz = 0;
}

void World::load_model_from_file(const std::string &filename) {
    std::ifstream file(filename);
    if(file) {
        qDebug("file ok: %s", filename.c_str());
        std::stringstream buffer;
        buffer << file.rdbuf();
        models.push_back(model3d(buffer));
        file.close();
    }
}

void World::load_assets()
{
    std::vector<std::string> filename = {
        "Handgun_obj.obj"
        //,"cube.obj"
        //,"oildrum.obj"
    };


    for(unsigned i=0; i< filename.size(); i++) {
        load_model_from_file(filename[i]);
        models[i].scx = 50;
        models[i].scy = 50;
        models[i].scz = 50;
        models[i].ty = 70;
    }

    lights.push_back(light_source(0, 100, -50, rgb(1.0, 1.0, 1.0), rgb(1, 1, 1), rgb(1.0, 1.0, 1.0)));
}
