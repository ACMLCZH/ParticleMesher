#pragma once

#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <openvdb/openvdb.h>

namespace particle_mesher {

struct FoamSpheres {
    std::vector<float> positions;
    std::string info_msg;
};

using openvdb::Vec3f;
using openvdb::Vec2f;

[[nodiscard]] inline static float clamp_normalize(
    const float val, const float min_val, const float max_val
) noexcept {
    return (std::min(val, max_val) - std::min(val, min_val)) / (max_val - min_val);
}

[[nodiscard]] inline static float cube(float x) noexcept { return x * x * x; }


[[nodiscard]] inline static float uniform_random() noexcept { 
    return static_cast<float>(std::rand()) / RAND_MAX;
}

[[nodiscard]] inline static float unit_sphere_cap_volume(const float h) noexcept {
    return h * h * (M_PI / 3) * (3 - h);
};
    
inline static void get_orthogonal(
    const Vec3f &vec, Vec3f &x, Vec3f &y
) noexcept {
    Vec3f v(1.0f, 0.0f, 0.0f);
    if (fabs(v.dot(vec)) > 0.999)
        v = Vec3f(0.0f, 1.0f, 0.0f);
    x = vec.cross(v);
    y = vec.cross(x);
    x.normalize();
    y.normalize();
}

class CubicKernel {
protected:
    float m_radius, m_k, m_l, m_W_zero;
public:
    CubicKernel(float radius) noexcept {
        const float h3 = cube(radius);
        m_radius = radius;
        m_k = 8.0f / (M_PI * h3);
        m_l = 48.0f / (M_PI * h3);
        m_W_zero = W(Vec3f(0.0f, 0.0f, 0.0f));
    }
    [[nodiscard]] float W(const float r) noexcept {
        const float q = r / m_radius;
        if (q <= 1.0)
            if (q <= 0.5)
                return m_k * (6.0f * (q - 1.0f) * q * q + 1.0f);
            else 
                return m_k * (2.0f * cube(1.0f - q));
        else
            return 0.0;
    }
    [[nodiscard]] float W(const Vec3f &r) noexcept { return W(r.length()); }
    [[nodiscard]] Vec3f gradW(const Vec3f &r) noexcept {
        Vec3f res(0.0f, 0.0f, 0.0f);
        const float rl = r.length();
        const float q = rl / m_radius;
        if ((rl > 1.0e-9) && (q <= 1.0f)) {
            Vec3f gradq = r / rl;
            gradq /= m_radius;
            if (q <= 0.5f)
                res = m_l * q * (3.0f * q - 2.0f) * gradq;
            else {
                const float p = 1.0f - q;
                res = m_l * (-p * p) * gradq;
            }
        }
        return res;
    }
    [[nodiscard]] float W_zero() noexcept { return m_W_zero; }
};

class OpenVDBFoamGenerator {

public:
    struct Config {
        float particle_radius;          // 0.025 -> 0.01
        float voxel_scale;
        float time_step;
        std::vector<float> lower_bound, upper_bound, gravity;
        // float ta_min, ta_max;               // 5, 20
        // float wc_min, wc_max;               // 2, 8
        // float ke_min, ke_max;               // 5, 50
        // float life_min, life_max;           // 2.0, 5.0
        std::vector<float> lim_ta, lim_wc, lim_ke, lim_life;
        float support_scale;                // 4.0
        int surface_neighbor_max;           // 20
        int generate_neighbor_min;          // 15
        int foam_neighbor_min, foam_neighbor_max;       // default 6, 20
        float k_ta, k_wc, k_bo, k_dr, k_ad, k_foam;     // 4000, 50000, 2.0(Buoyancy), 0.8(Drag), 0.99(Air damp), 1000.0(Foam scale)
        float spray_decay, foam_decay, bubble_decay;    // 2, 1, 0
        float foam_density;                             // 1000.0
    };

private:
    enum FoamType { Spray, Foam, Bubble, NUM_PARTICLE_TYPES };

private:
    std::string _object_id;     // Just a reminder that you cannot apply one generator to different objects.
    CubicKernel _kernel;
    float _particle_radius, _particle_mass, _support_radius, _voxel_size, _time_step;
    size_t _neighbor_search;
    size_t _surface_neighbor_max;
    size_t _generate_neighbor_min;
    size_t _foam_neighbor_min, _foam_neighbor_max;
    float _k_ta, _k_wc, _k_bo, _k_dr, _k_ad, _k_foam;
    float _spray_decay, _foam_decay, _bubble_decay;
    Vec2f _lim_ta, _lim_wc, _lim_ke, _lim_life;
    Vec3f _gravity, _lower_bound, _upper_bound;
    std::vector<Vec3f> _tot_pFoams, _tot_vFoams;
    std::vector<float> _tot_lifeFoams;
    std::vector<FoamType> _tot_typeFoams;

public:
    OpenVDBFoamGenerator(const Config &config, const std::string_view object_id) noexcept;
    [[nodiscard]] FoamSpheres generate_foams(
        const std::vector<float> &positions,
        const std::vector<float> &velocities
    );
};

}  // namespace particle_mesher