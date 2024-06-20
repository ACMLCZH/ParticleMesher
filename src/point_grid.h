#pragma once

#include <openvdb/math/Coord.h>
#include <openvdb/points/PointConversion.h>
#include <openvdb/tree/LeafManager.h>
#include <openvdb/tools/PointIndexGrid.h>
#include <openvdb/tools/ParticlesToLevelSet.h>

namespace particle_mesher {

using namespace openvdb::points;
using namespace openvdb::tree;
using namespace openvdb::tools;
using namespace openvdb::math;
using openvdb::Vec3f;
using openvdb::Vec3R;
using openvdb::Real;
using openvdb::PointIndex32;
using openvdb::FloatTree;
using openvdb::FloatGrid;
using openvdb::createLevelSet;

class ParticleList {
private:
    std::vector<Vec3R> _positions;
    std::vector<Real> _radii;
    float _radius_min, _radius_max;

public:
    typedef Vec3R PosType;
    ParticleList(const std::vector<float> &positions, const std::vector<float> &radii, float radius) noexcept {
        auto particle_count = positions.size() / 3u;
        auto global_radius = radii.size() == 0u;
        _positions.reserve(particle_count);
        _radii.reserve(particle_count);
        _radius_min = 1e9;
        _radius_max = 0;
        for (auto i = 0u; i < particle_count; ++i) {
            const auto bi = i * 3u;
            _positions.emplace_back(positions[bi], positions[bi + 1u], positions[bi + 2u]);
            float r = (global_radius ? radius : radii[i]);
            _radii.emplace_back(r);
            _radius_min = std::min(_radius_min, r);
            _radius_max = std::max(_radius_max, r);
        }
    }
    [[nodiscard]] size_t size() const noexcept {
        return _positions.size();
    }
    void clear() noexcept {
        _positions.clear();
        _radii.clear();
    }
    void getPos(size_t n, Vec3R &pos) const noexcept {   // needed by rasterizeSpheres(pa, r)
        pos = _positions[n];
    }
    void getPosRad(size_t n, Vec3R &pos, Real &radius) const noexcept {    // need by rasterizeSpheres(pa)
        pos = _positions[n];
        radius = _radii[n];
    }
    [[nodiscard]] float radius_min() const noexcept { return _radius_min; }
    [[nodiscard]] float radius_max() const noexcept { return _radius_max; }
};

[[nodiscard]] static FloatGrid::Ptr build_level_set(
    const ParticleList &particle_list, float voxel_size, float half_width = 3.f
) noexcept {
    auto sdf = createLevelSet<FloatGrid>(voxel_size, half_width);
    ParticlesToLevelSet<FloatGrid> rasterizer(*sdf);
    rasterizer.setRmin(particle_list.radius_min() / voxel_size * 0.9);     // Rmin < R(voxel unit) < Rmax
    rasterizer.setRmax(particle_list.radius_max() / voxel_size * 1.1);     // in voxel unit
    rasterizer.rasterizeSpheres(particle_list);
    pruneLevelSet(sdf->tree());
    return std::move(sdf);
}

[[nodiscard]] static std::vector<bool> find_surface(
    const FloatGrid &level_set, const ParticleList &query_list
) noexcept {
    const Transform &transform = level_set.transform();
    auto query_count = query_list.size();
    std::vector<bool> is_surface(query_count);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, query_count),
        [&](tbb::blocked_range<size_t> rng) {
            auto level_set_accessor = level_set.getConstAccessor();
            GridSampler<ValueAccessor<const FloatTree>, BoxSampler> sampler(level_set_accessor, transform);
            for (size_t i = rng.begin(), N = rng.end(); i < N; ++i) {
                Vec3R particle;
                Real radius;
                query_list.getPosRad(i, particle, radius);
                float value = sampler.wsSample(particle);
                // printf("(%f, %f, %f), value = %f, radius = %f\n",
                //         particle(0), particle(1), particle(2), value, radius);
                is_surface[i] = std::abs(value) <= radius;
            }
        }       
    );
    
    return std::move(is_surface);
}

[[nodiscard]] static PointIndexGrid::Ptr build_index_grid(
    const std::vector<Vec3f> &positions, float voxel_size
) noexcept {
    PointAttributeVector<Vec3f> pWrapper(positions);
    auto transform = Transform::createLinearTransform(voxel_size);
    return createPointIndexGrid<PointIndexGrid>(pWrapper, *transform);
}

[[nodiscard]] static std::vector<std::vector<size_t> > find_neighbors(
    const std::vector<Vec3f> &points, const PointIndexGrid &point_index_grid, float support_radius,
    const std::vector<Vec3f> &query_points
) noexcept {
    const Transform &transform = point_index_grid.transform();
    size_t neighbor_search = (size_t)(support_radius / transform.voxelSize()(0));
    auto query_count = query_points.size();
    std::vector<std::vector<size_t> > query_neighbor_indices(query_count);

    // Find point neighbors
    tbb::parallel_for(tbb::blocked_range<size_t>(0, query_count),
        [&](tbb::blocked_range<size_t> rng) {
            auto point_index_accessor = point_index_grid.getConstAccessor();
            for (size_t i = rng.begin(), N = rng.end(); i < N; ++i) {
                const Vec3f &xi = query_points[i];
                const Coord ijk = transform.worldToIndexNodeCentered(xi);
                CoordBBox neighborBox(ijk.offsetBy(-neighbor_search), ijk.offsetBy(neighbor_search));
                const Coord &na(neighborBox.min()), &nb(neighborBox.max());
                
                // Find particle neighbors
                for (auto x = na.x(); x <= nb.x(); ++x)
                    for (auto y = na.y(); y <= nb.y(); ++y)
                        for (auto z = na.z(); z <= nb.z(); ++z) {
                            auto nijk = Coord(x, y, z);
                            auto nleafPtr = point_index_accessor.probeConstLeaf(nijk);
                            if (!nleafPtr) continue;
                            const PointIndex32 *nbegin(nullptr), *nend(nullptr);
                            nleafPtr->getIndices(nijk, nbegin, nend);
                            for (;nbegin < nend; ++nbegin) {
                                const Vec3f &xj = points[*nbegin];
                                Vec3f xixj = xi - xj;
                                if (xixj.length() < support_radius) 
                                    query_neighbor_indices[i].emplace_back(*nbegin);
                            }
                        }
           }
        }
    );
    return std::move(query_neighbor_indices);
}

}  // namespace particle_mesher