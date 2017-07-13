#ifndef WORLD_H
#define WORLD_H

#include "camera.h"
#include "model3d.h"
#include "light_source.h"

struct World{
    camera cam;
    std::vector<model3d> models;
    std::vector<light_source> lights;


    World();
    World(const World&) = delete;

    void load_assets();
    void load_model_from_file(const std::string& filename);


};

#endif // WORLD_H
