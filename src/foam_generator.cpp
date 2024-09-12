#include <tbb/parallel_for.h>
#include <openvdb/math/Coord.h>
#include <openvdb/points/PointConversion.h>
#include <openvdb/tree/LeafManager.h>
#include <openvdb/tools/PointIndexGrid.h>

#include "clock.h"
#include "format.h"
#include "foam_generator.h"
#include "point_grid.h"

namespace particle_mesher {

using namespace openvdb::points;
using namespace openvdb::tree;
using namespace openvdb::tools;
using namespace openvdb::math;
using openvdb::Vec3f;
using openvdb::PointIndex32;

OpenVDBFoamGenerator::OpenVDBFoamGenerator(
    const OpenVDBFoamGenerator::Config &config, const std::string_view object_id
) noexcept :
    _object_id(object_id),
    _kernel(config.particle_radius * config.support_scale),
    _particle_radius(config.particle_radius),
    _particle_mass(6.4f * cube(config.particle_radius) * config.foam_density),
    _support_radius(config.particle_radius * config.support_scale),
    _voxel_size(config.particle_radius * config.voxel_scale),
    _time_step(config.time_step),
    _neighbor_search(size_t(config.support_scale / config.voxel_scale)),
    _surface_neighbor_max(config.surface_neighbor_max),
    _generate_neighbor_min(config.generate_neighbor_min),
    _foam_neighbor_min(config.foam_neighbor_min), _foam_neighbor_max(config.foam_neighbor_max),
    _k_ta(config.k_ta), _k_wc(config.k_wc), _k_bo(config.k_bo), _k_dr(config.k_dr), _k_ad(config.k_ad), _k_foam(config.k_foam),
    _spray_decay{config.spray_decay}, _foam_decay{config.foam_decay}, _bubble_decay{config.bubble_decay},
    _lim_ta(config.lim_ta[0], config.lim_ta[1]),
    _lim_wc(config.lim_wc[0], config.lim_wc[1]),
    _lim_ke(config.lim_ke[0], config.lim_ke[1]),
    _lim_life(config.lim_life[0], config.lim_life[1]),
    _gravity(config.gravity[0], config.gravity[1], config.gravity[2]),
    _lower_bound(config.lower_bound[0], config.lower_bound[1], config.lower_bound[2]),
    _upper_bound(config.upper_bound[0], config.upper_bound[1], config.upper_bound[2])
{
    openvdb::initialize();
}

FoamSpheres OpenVDBFoamGenerator::generate_foams(
    const std::vector<float> &positions,
    const std::vector<float> &velocities
) {
    if (positions.size() % 3u != 0u ||
        velocities.size() % 3u != 0u ||
        positions.size() != velocities.size())
        throw std::logic_error("Invalid particle count.");

    Clock clock;
    std::vector<std::string> info_pack;

    auto particle_count = positions.size() / 3u;
    std::vector<Vec3f> p(particle_count);
    std::vector<Vec3f> v(particle_count);
    std::vector<Vec3f> normals(particle_count);
    std::vector<float> densities(particle_count);
    std::vector<std::vector<Vec3f> > pFoams(particle_count);
    std::vector<std::vector<Vec3f> > vFoams(particle_count);
    std::vector<std::vector<float> > lifeFoams(particle_count);
    std::vector<size_t> foam_offsets(particle_count);

    info_pack.push_back(str_format(
        "Generator initialized with:\n"
        "\ttime step = %f, particle radius = %f, voxel size = %f, support radius = %f, mass = %f, neighbor offset = %u\n"
        "\tlower = (%f, %f, %f), upper = (%f, %f, %f), gravity = (%f, %f, %f)\n"
        "\tta: (k = %f, range = (%f, %f)), wc: (k = %f, range = (%f, %f)), ke: (range = (%f, %f))\n"
        "\tfoam scale = %f, boyancy = %f, drag = %f, air_damp = %f\n"
        "\tspray_decay = %f, foam_decay = %f, bubble_decay = %f\n",
        _time_step, _particle_radius, _voxel_size, _support_radius, _particle_mass, _neighbor_search,
        _lower_bound(0), _lower_bound(1), _lower_bound(2),
        _upper_bound(0), _upper_bound(1), _upper_bound(2), 
        _gravity(0), _gravity(1), _gravity(2),
        _k_ta, _lim_ta(0), _lim_ta(1), _k_wc, _lim_wc(0), _lim_wc(1), _lim_ke(0), _lim_ke(1),
        _k_foam, _k_bo, _k_dr, _k_ad,
        _spray_decay, _foam_decay, _bubble_decay
    ));

    // Populate positions and velocities
    tbb::parallel_for(tbb::blocked_range<size_t>(0, particle_count),
        [&](tbb::blocked_range<size_t> rng) {
            for (size_t i = rng.begin(), N = rng.end(); i < N; ++i) {
                const size_t bi = i * 3u;
                p[i] = Vec3f(positions[bi],  positions[bi + 1u],  positions[bi + 2u]);
                v[i] = Vec3f(velocities[bi], velocities[bi + 1u], velocities[bi + 2u]);
            }
        }
    );

    info_pack.push_back(str_format(
        "Populate positions and velocities in %f ms: \n"
        "\tparticles count = %u\n",
        clock.toc(), particle_count
    ));

    // Find fluid neighbors
    auto point_index_grid_ptr = build_index_grid(p, _voxel_size);
    auto fluid_neighbor_indices = find_neighbors(p, *point_index_grid_ptr, _support_radius, p);
    
    // auto all_count = 0u;
    // for (auto i = 0u; i < particle_count; ++i) {
    //     all_count += fluid_neighbor_indices[i].size();
    // }
    // "\tAverage neighbors count = %f\n",
    info_pack.push_back(str_format("Find neighbors in %f ms:\n", clock.toc()));
    
    // Calculate densities
    // density = mass [r^3] * w [r^-3] -> unchanged with r
    tbb::parallel_for(tbb::blocked_range<size_t>(0, particle_count),
        [&](tbb::blocked_range<size_t> rng) {
            for (size_t i = rng.begin(), N = rng.end(); i < N; ++i) {
                const Vec3f &xi = p[i];
                const std::vector<size_t> &neighbors = fluid_neighbor_indices[i];
                
                float &di = densities[i];
                di = _particle_mass * _kernel.W_zero();
                for (auto j: neighbors) {
                    const Vec3f &xj = p[j];
                    di += _particle_mass * _kernel.W(xi - xj);
                }
            }
        }
    );

    info_pack.push_back(str_format("Calculate densities in %f ms.\n", clock.toc()));

    // Calculate normals
    tbb::parallel_for(tbb::blocked_range<size_t>(0, particle_count),
        [&](tbb::blocked_range<size_t> rng) {
            for (size_t i = rng.begin(), N = rng.end(); i < N; ++i) {
                const Vec3f &xi = p[i];
                const std::vector<size_t> &neighbors = fluid_neighbor_indices[i];
                if (neighbors.size() > _surface_neighbor_max) continue;

                Vec3f &ni = normals[i];
                for (auto j: neighbors) {
                    const Vec3f &xj = p[j];
                    const float &dj = densities[j];
                    ni -= _particle_mass / dj * _kernel.gradW(xi - xj);
                }
                ni.normalize();
            }
        }
    );

    info_pack.push_back(str_format("Calculate normals in %f ms.\n", clock.toc()));

    // Generate foams
    tbb::parallel_for(tbb::blocked_range<size_t>(0, particle_count),
        [&](tbb::blocked_range<size_t> rng) {
            for (size_t i = rng.begin(), N = rng.end(); i < N; ++i) {
                const Vec3f &xi = p[i];
                const Vec3f &vi = v[i];
                const Vec3f &ni = normals[i];
                const std::vector<size_t> &neighbors = fluid_neighbor_indices[i];
                if (neighbors.size() < _generate_neighbor_min) 
                    continue;

                float v_diff = 0.0, curvature = 0.0;
                for (auto j: neighbors) {
                    const Vec3f &xj = p[j];
                    const Vec3f &vj = v[j];
                    const Vec3f &nj = normals[j];
                    const float dj = densities[j];

                    // Trapped air potential 
                    // Compute Eq. 2 in Ihmsen et al., "Unified spray, foam and air bubbles for particle-based fluids", 2012
                    Vec3f vivj = vi - vj;
                    const float vmag = vivj.length();
                    vivj.normalize(1.0e-6);
                    Vec3f xixj = xi - xj;
                    xixj.normalize();
                    v_diff += _particle_mass / dj * vmag * (1.0f - vivj.dot(xixj)) * _kernel.W(xi - xj);

                    // Wave crest curvature
                    // Compute Eq. 4 in Ihmsen et al., "Unified spray, foam and air bubbles for particle-based fluids", 2012
                    if (-xixj.dot(ni) < 0)
                        curvature += _particle_mass / dj * (1.0f - ni.dot(nj)) * _kernel.W(xi - xj);
                }

                // Trapped air potential 
                const float I_ta = clamp_normalize(v_diff, _lim_ta.x(), _lim_ta.y());

                // Wave crest
                // Compute Eq. 7 in Ihmsen et al., "Unified spray, foam and air bubbles for particle-based fluids", 2012
                const float vi_norm = vi.length();
                const Vec3f vi_normed = vi * (1.0f / vi_norm);
                const float delta = vi_normed.dot(ni) >= 0.6 ? 1.0 : 0.0;
                const float I_wc = clamp_normalize(delta * curvature, _lim_wc.x(), _lim_wc.y());

                // Kinetic energy
                const float ke = 0.5f * _particle_mass * vi_norm * vi_norm;
                const float I_ke = clamp_normalize(ke, _lim_ke.x(), _lim_ke.y());

                const unsigned int nd = (unsigned int)(std::max(
                    (_k_ta * I_ta + _k_wc * I_wc) * _k_foam * I_ke * _time_step + 0.5f, 0.0f
                ));

                // if (nd > 0) {
                //     printf("%.12f %.12f %.12f %u\n", v_diff, delta * curvature, ke, nd);
                // }

                // Generate foam particles
                Vec3f e1, e2;
                get_orthogonal(vi_normed, e1, e2);
                e1 = _particle_radius * e1;
                e2 = _particle_radius * e2;
                
                pFoams[i].reserve(nd);
                vFoams[i].reserve(nd);
                lifeFoams[i].reserve(nd);
                
                // Generate a random distribution of the foam particles in a cylinder.
                for (unsigned int j = 0; j < nd; ++j)
                {
                    const float Xr = uniform_random();
                    const float Xtheta = uniform_random();
                    const float Xh = uniform_random();
                    const float Xlt = uniform_random();

                    const float r = _particle_radius * sqrt(Xr);
                    const float theta = Xtheta * (2.0f * M_PI);
                    const float h = (Xh - 0.5f) * _time_step * vi_norm;

                    const Vec3f xd = xi + r * cos(theta) * e1 + r * sin(theta) * e2 + h * vi_normed;
                    const Vec3f vd = r * cos(theta) * e1 + r * sin(theta) * e2 + vi;
                    const float lt = _lim_life.x() + I_ke * Xlt * (_lim_life.y() - _lim_life.x());

                    pFoams[i].emplace_back(xd);
                    vFoams[i].emplace_back(vd);
                    lifeFoams[i].emplace_back(lt);
                }
            }
        }
    );

    info_pack.push_back(str_format("Generate foams in %f ms.\n", clock.toc()));

    // append to foams to total
    size_t foam_count = _tot_pFoams.size(), new_foam_count = 0;
    for (int i = 0; i < particle_count; ++i) {
        new_foam_count += pFoams[i].size();
        foam_offsets[i] = new_foam_count;
    }
    size_t tot_foam_count = foam_count + new_foam_count;
    _tot_pFoams.resize(tot_foam_count);
    _tot_vFoams.resize(tot_foam_count);
    _tot_lifeFoams.resize(tot_foam_count);
    _tot_typeFoams.resize(tot_foam_count);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, particle_count),
        [&](tbb::blocked_range<size_t> rng) {
            for (size_t i = rng.begin(), N = rng.end(); i < N; ++i) {
                size_t start = i == 0 ? 0 : foam_offsets[i - 1];
                size_t size = pFoams[i].size();
                for (size_t j = 0; j < size; ++j) {
                    size_t cur_j = foam_count + start + j;
                    _tot_pFoams[cur_j] = pFoams[i][j];
                    _tot_vFoams[cur_j] = vFoams[i][j];
                    _tot_lifeFoams[cur_j] = lifeFoams[i][j];
                }
            }
        }
    );

    info_pack.push_back(str_format(
        "Append generated foam particles in %f ms:\n"
        "\tfoam count from %u to %u\n",
        clock.toc(), foam_count, tot_foam_count
    ));

    auto foam_neighbor_indices = find_neighbors(
        p, *point_index_grid_ptr, _support_radius, _tot_pFoams
    );
    std::vector<size_t> foams_num_neighbors(tot_foam_count);

    // advect foams
    tbb::parallel_for(tbb::blocked_range<size_t>(0, tot_foam_count),
        [&](tbb::blocked_range<size_t> rng) {
            // auto pointIndexAccessor = pointIndexGridPtr->getAccessor();
            for (size_t i = rng.begin(), N = rng.end(); i < N; ++i) {
                const Vec3f &xi = _tot_pFoams[i];

                // Correct neighbor counts
                foams_num_neighbors[i] = foam_neighbor_indices[i].size();
                if (foams_num_neighbors[i] >= _foam_neighbor_min &&
                    foams_num_neighbors[i] <= _foam_neighbor_max) {
                    auto xmin = xi - _lower_bound, xmax = _upper_bound - xi;
                    for (auto j = 0u; j < 3; ++j) {
                        float dis;
                        dis = std::abs(xmin(j));
                        if (dis < _support_radius) 
                            foams_num_neighbors[i] += (size_t)(std::ceil(
                                foams_num_neighbors[i] * unit_sphere_cap_volume(1 - dis / _support_radius)
                            ));
                        dis = std::abs(xmax(j));
                        if (dis < _support_radius)
                            foams_num_neighbors[i] += (size_t)(std::ceil(
                                foams_num_neighbors[i] * unit_sphere_cap_volume(1 - dis / _support_radius)
                            ));
                    }
                }

                // Handle boundary foams
                for (auto j = 0u; j < 3; ++j) {
                    if (xi(j) < _lower_bound(j) || xi(j) > _upper_bound(j)) {
                        _tot_lifeFoams[i] = -1.0f;
                        break;
                    }
                }

    			// Advect foam particles
                if (foams_num_neighbors[i] < _foam_neighbor_min) {   // spray
                    _tot_typeFoams[i] = FoamType::Spray;
    				_tot_vFoams[i] += _time_step * _gravity;
    				_tot_vFoams[i] *= _k_ad;        // air damping
                    _tot_pFoams[i] += _time_step * _tot_vFoams[i];
    				if (foams_num_neighbors[i] < 2)
    					_tot_lifeFoams[i] = -1.0f;
    				else
    					_tot_lifeFoams[i] -= _spray_decay * _time_step;
                } else {
                    Vec3f sv(0.f, 0.f, 0.f);
                    float sk = 0.0;
                    for (auto nIndex: foam_neighbor_indices[i]) {
                        const Vec3f &xj = p[nIndex];
                        const float kj = _kernel.W(xi - xj);
                        const Vec3f &vj = v[nIndex]; 
                        sv += vj * kj;
                        sk += kj;
                    }
                    sv = sv / sk;

                    if (foams_num_neighbors[i] <= _foam_neighbor_max) {      // foam
                        _tot_typeFoams[i] = FoamType::Foam;
                        _tot_vFoams[i] = sv;
                        _tot_pFoams[i] += sv * _time_step;
                        _tot_lifeFoams[i] -= _foam_decay * _time_step;
                    } else {                                            // bubble
                        _tot_typeFoams[i] = FoamType::Bubble;
                        _tot_vFoams[i] += _time_step * _gravity * -_k_bo + (sv - _tot_vFoams[i]) * _k_dr;
                        _tot_pFoams[i] += _time_step * _tot_vFoams[i];
                        _tot_lifeFoams[i] -= _bubble_decay * _time_step;
                    }
                }
            }
        }
    );

    info_pack.push_back(str_format(
        "Recognize and advect foam particles in %f ms:\n", clock.toc()
    ));

    // remove foam particles
	auto removed_foam_count = 0u;
    std::vector<size_t> render_indices;
	for (auto i = 0u; i < tot_foam_count; ++i)
	{
		if (_tot_lifeFoams[i] <= 0.0)
            ++removed_foam_count;
		else {
            size_t ni = i - removed_foam_count;
			_tot_pFoams[ni] = _tot_pFoams[i];
			_tot_vFoams[ni] = _tot_vFoams[i];
			_tot_lifeFoams[ni] = _tot_lifeFoams[i];
			_tot_typeFoams[ni] = _tot_typeFoams[i];
            foams_num_neighbors[ni] = foams_num_neighbors[i];
            if (foams_num_neighbors[ni] >= 2)
                render_indices.emplace_back(ni);
		}
	}
	if (removed_foam_count > 0u) {
		_tot_pFoams.resize(tot_foam_count - removed_foam_count);
		_tot_vFoams.resize(tot_foam_count - removed_foam_count);
		_tot_lifeFoams.resize(tot_foam_count - removed_foam_count);
		_tot_typeFoams.resize(tot_foam_count - removed_foam_count);
	}

    info_pack.push_back(str_format(
        "Remove foam particles in %f ms:\n"
        "\tremoved foam count = %u\n",
        clock.toc(), removed_foam_count
    ));

    // Copy foam particles into the struct
    auto render_count = render_indices.size();
    FoamSpheres foams;
    foams.positions.resize(render_count * 3);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, render_count),
        [&](tbb::blocked_range<size_t> rng) {
            for (size_t i = rng.begin(), N = rng.end(); i < N; ++i) {
                size_t bi = i * 3u, tot_index = render_indices[i];
                foams.positions[bi     ] = _tot_pFoams[tot_index].x();
                foams.positions[bi + 1u] = _tot_pFoams[tot_index].y();
                foams.positions[bi + 2u] = _tot_pFoams[tot_index].z();
            }
        }
    );

    info_pack.push_back(str_format(
        "Copy foam particles in %f ms:\n"
        "\trender count (neighbor >= 2) = %u\n",
        clock.toc(), render_count
    ));

    combine_message(info_pack, foams.info_msg);

    return foams;
}

}  // namespace particle_mesher