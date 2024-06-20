#pragma once

#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

#include <openvdb/openvdb.h>

namespace particle_mesher {

struct SurfaceIndices {
    std::vector<bool> is_surface;
    std::string info_msg;
};

class OpenVDBSurfaceSplitter {

public:
    struct Config {
        float particle_radius;
        float voxel_scale;
        float support_scale;
        float half_width;
        size_t surface_neighbor_max;
    };

private:
    float _particle_radius, _support_radius, _voxel_size, _half_width;
    size_t _surface_neighbor_max;

public:
    OpenVDBSurfaceSplitter(const Config &config) noexcept :
        _particle_radius(config.particle_radius),
        _support_radius(config.particle_radius * config.support_scale),
        _voxel_size(config.particle_radius * config.voxel_scale),
        _surface_neighbor_max(config.surface_neighbor_max),
        _half_width(config.half_width)
    {
        openvdb::initialize();
    }
    [[nodiscard]] SurfaceIndices split_surface_count (
        const std::vector<float> &positions
    );
    [[nodiscard]] SurfaceIndices split_surface_sdf (
        const std::vector<float> &positions,
        const std::vector<float> &radii
    );
};

}  // namespace particle_mesher