#include <exception>
#include <openvdb/tools/VolumeToMesh.h>

#include "clock.h"
#include "format.h"
#include "mesh_construct.h"
#include "point_grid.h"

namespace particle_mesher {

using namespace openvdb::tools;
using openvdb::Vec3s;
using openvdb::Vec3I;
using openvdb::Vec4I;

ConstructMesh OpenVDBMeshConstructor::construct(
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
        "\tparticles count = %u, particle radius = %f, voxel size = %f\n",
        clock.toc(), pa.size(), _particle_radius, _voxel_size
    ));
    
    auto sdf = build_level_set(pa, _voxel_size);
    info_pack.push_back(str_format("Particles to SDF in %f ms.\n", clock.toc()));

    std::vector<Vec3s> points;
    std::vector<Vec3I> tris;
    std::vector<Vec4I> quads;
    volumeToMesh(*sdf, points, tris, quads, _isovalue, _adaptivity, true);

    // puts("CONST 3!!!!!!!!!!!!!");
    info_pack.push_back(str_format("SDF to Mesh in %f ms.\n", clock.toc()));

    ConstructMesh mesh;
    mesh.vertices.resize(points.size() * 3u);
    mesh.triangles.resize((tris.size() + quads.size() * 2u) * 3u);

    for (auto i = 0u; i < points.size(); ++i) {
        auto point = points[i];
        auto bi = i * 3u;
        mesh.vertices[bi    ] = point[0];
        mesh.vertices[bi + 1] = point[1];
        mesh.vertices[bi + 2] = point[2];
    }

    for (auto i = 0u; i < tris.size(); ++i) {
        auto tri = tris[i];
        auto bi = i * 3u;
        mesh.triangles[bi    ] = tri[0];
        mesh.triangles[bi + 1] = tri[1];
        mesh.triangles[bi + 2] = tri[2];
    }
    for (auto i = 0u; i < quads.size(); ++i) {
        auto quad = quads[i];
        auto bi = tris.size() * 3u + i * 6u;
        mesh.triangles[bi    ] = quad[0];
        mesh.triangles[bi + 1] = quad[1];
        mesh.triangles[bi + 2] = quad[2];
        mesh.triangles[bi + 3] = quad[0];
        mesh.triangles[bi + 4] = quad[2];
        mesh.triangles[bi + 5] = quad[3];
    }

    // puts("CONST 4!!!!!!!!!!!!!");
    info_pack.push_back(str_format(
        "Reconstruct mesh surface from particles using OpenVDB in %f ms:\n"
        "\tvertices count = %u, triangles count = %u",
        clock.toc(), mesh.vertices.size(), mesh.triangles.size()
    ));
    combine_message(info_pack, mesh.info_msg);

    return mesh;
}

}  // namespace particle_mesher