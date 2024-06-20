#pragma once

#include <string>
#include <vector>
#include <openvdb/openvdb.h>


namespace particle_mesher {

struct ConstructMesh {
    std::vector<float> vertices;
    std::vector<unsigned int> triangles;
    std::string info_msg;
};

class OpenVDBMeshConstructor {

public:
    struct Config  {
        float particle_radius;
        float voxel_scale;
        float isovalue;
        float adaptivity;
    };

private:
    float _particle_radius;
    float _voxel_size;
    float _isovalue;
    float _adaptivity;

public:
    OpenVDBMeshConstructor(const Config &config) noexcept :
        _particle_radius{config.particle_radius},
        _voxel_size{config.particle_radius * config.voxel_scale},
        _isovalue{config.isovalue},
        _adaptivity{config.adaptivity}
    { 
        openvdb::initialize();
    }
    [[nodiscard]] ConstructMesh construct(
        const std::vector<float> &positions,
        const std::vector<float> &radii
    );
};

}  // namespace particle_mesher