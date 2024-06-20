#include <tbb/parallel_for.h>

#include "clock.h"
#include "format.h"
#include "surface_splitter.h"
#include "point_grid.h"

namespace particle_mesher {

SurfaceIndices OpenVDBSurfaceSplitter::split_surface_count(
    const std::vector<float> &positions
) {
    if (positions.size() % 3u != 0u)
        throw std::logic_error("Invalid particle count.");

    Clock clock;
    std::vector<std::string> info_pack;
    auto particle_count = positions.size() / 3u;

    // Populate positions (when count > 100000, tbb is faster)
    std::vector<Vec3f> p(particle_count);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, particle_count),
        [&](tbb::blocked_range<size_t> rng) {
            for (size_t i = rng.begin(), N = rng.end(); i < N; ++i) {
                const size_t bi = i * 3u;
                p[i] = Vec3f(positions[bi], positions[bi + 1u], positions[bi + 2u]);
            }
        }
    );

    // std::vector<Vec3f> p;
    // p.reserve(particle_count);
    // for (auto i = 0u; i < particle_count; ++i) {
    //     const auto bi = i * 3u;
    //     p.emplace_back(positions[bi], positions[bi + 1u], positions[bi + 2u]);
    // }

    info_pack.push_back(str_format(
        "Populate positions in %f ms:\n"
        "\tparticles count = %u\n",
        clock.toc(), particle_count
    ));

    // Find point neighbors
    auto point_index_grid_ptr = build_index_grid(p, _voxel_size);
    auto point_neighbor_indices = find_neighbors(p, *point_index_grid_ptr, _support_radius, p);

    SurfaceIndices surface_indices;
    surface_indices.is_surface.resize(particle_count);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, particle_count),
        [&](tbb::blocked_range<size_t> rng) {
            for (size_t i = rng.begin(), N = rng.end(); i < N; ++i) {
                surface_indices.is_surface[i] = point_neighbor_indices[i].size() < _surface_neighbor_max;
            }
        }
    );

    info_pack.push_back(str_format("Copy surface indices in %f ms:\n", clock.toc()));
    combine_message(info_pack, surface_indices.info_msg);

    return surface_indices;
}


SurfaceIndices OpenVDBSurfaceSplitter::split_surface_sdf(
    const std::vector<float> &positions,
    const std::vector<float> &radii
) {
    if (positions.size() % 3u != 0u ||
        radii.size() != 0u && radii.size() * 3 != positions.size())
        throw std::logic_error("Invalid particle count.");

    Clock clock;
    std::vector<std::string> info_pack;
    ParticleList pa(positions, radii, _particle_radius);
    info_pack.push_back(str_format(
        "Add particles in %f ms:\n"
        "\tparticles count = %u, particle radius = %f, voxel size = %f, half width = %f\n",
        clock.toc(), pa.size(), _particle_radius, _voxel_size, _half_width
    ));
    puts(info_pack[info_pack.size() - 1].c_str());

    auto sdf = build_level_set(pa, _voxel_size, _half_width);
    info_pack.push_back(str_format("Particles to SDF in %f ms.\n", clock.toc()));
    puts(info_pack[info_pack.size() - 1].c_str());

    SurfaceIndices surface_indices{find_surface(*sdf, pa), ""};
    info_pack.push_back(str_format("Copy surface indices in %f ms:\n", clock.toc()));
    combine_message(info_pack, surface_indices.info_msg);

    return surface_indices;
}

}  // namespace particle_mesher