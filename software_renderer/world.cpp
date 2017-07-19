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
        "table lamp.obj",
        //"Handgun_obj.obj"
        //"Street environment_V01.obj"
        //,"cube.obj"
        "oildrum.obj"

        //"chair.obj",
       // "table.obj",
    };

    std::unordered_map<std::string, float> scale;
    scale["table lamp.obj"] = 1;
    scale["oildrum.obj"] = 100;
    scale["chair.obj"] = 50;
    scale["chair2.obj"] = 1.0/8;
    scale["table.obj"] = 50;
    scale["cube.obj"] = 500;
    scale["desk2.obj"] = 20;
    scale["ModernDeskOBJ.obj"] = 2;
    scale["fruit stand.obj"] = 1.0/50;

    for(unsigned i=0; i< filename.size(); i++) {
        load_model_from_file(filename[i]);
        if(scale.find(filename[i]) != scale.end()) {
            float s = scale[filename[i]];
            models[i].scx = s;
            models[i].scy = s;
            models[i].scz = s;
        }
    }

    lights.push_back(light_source(30, 70, 0, rgb(1, 1, 1), rgb(1, 1 ,1), rgb(1, 1, 1), &models[0]));
    //lights.push_back(light_source(0, -200, 0, rgb(1, 1, 1), rgb(0, 1, 0), rgb(1.0, 1.0, 1.0)));
    //lights.push_back(light_source(0, 200, 0, rgb(1, 1, 1), rgb(0, 0, 1), rgb(1.0, 1.0, 1.0)));

}
