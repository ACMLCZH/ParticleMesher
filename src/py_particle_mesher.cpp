#include <vector>
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "foam_generator.h"
#include "mesh_construct.h"
#include "surface_splitter.h"

using namespace particle_mesher;
namespace py = pybind11;
using PyFloatArr = py::array_t<float>;
using PyUIntArr = py::array_t<unsigned int>;
using PyBoolArr = py::array_t<bool>;
using StdFloatArr = std::vector<float>;

struct PyConstructMesh {
    PyFloatArr vertices;
    PyUIntArr triangles;
    std::string info_msg;
};

struct PyFoamSpheres {
    PyFloatArr positions;
    std::string info_msg;
};

struct PySurfaceIndices {
    PyBoolArr is_surface;
    std::string info_msg;
};

class PyMeshConstructor : public OpenVDBMeshConstructor {

public:
    PyMeshConstructor(const OpenVDBMeshConstructor::Config &config) noexcept :
        OpenVDBMeshConstructor(config) { }
    PyConstructMesh construct(const PyFloatArr &positions, const PyFloatArr &radii) {
        std::vector<float> positions_vec(positions.data(), positions.data() + positions.size());
        std::vector<float> radii_vec(radii.data(), radii.data() + radii.size());
        auto construct_mesh = OpenVDBMeshConstructor::construct(positions_vec, radii_vec);
        return {
            PyFloatArr(construct_mesh.vertices.size(), construct_mesh.vertices.data()),
            PyUIntArr(construct_mesh.triangles.size(), construct_mesh.triangles.data()),
            construct_mesh.info_msg
        };
    }
};

class PyFoamGenerator : public OpenVDBFoamGenerator {
    
public:
    PyFoamGenerator(const OpenVDBFoamGenerator::Config &config, const std::string_view object_id) noexcept : 
        OpenVDBFoamGenerator(config, object_id) { }
    PyFoamSpheres generate_foams(
        const PyFloatArr &positions, const PyFloatArr &velocities
    ) {
        std::vector<float> positions_vec(positions.data(), positions.data() + positions.size());
        std::vector<float> velocities_vec(velocities.data(), velocities.data() + velocities.size());
        auto foam_spheres = OpenVDBFoamGenerator::generate_foams(positions_vec, velocities_vec);

        return {
            PyFloatArr(foam_spheres.positions.size(), foam_spheres.positions.data()),
            foam_spheres.info_msg
        };
    }
};

class PySurfaceSplitter : public OpenVDBSurfaceSplitter {

public:
    PySurfaceSplitter(const OpenVDBSurfaceSplitter::Config &config) noexcept :
        OpenVDBSurfaceSplitter(config) { }
    PySurfaceIndices split_surface_sdf(const PyFloatArr &positions, const PyFloatArr &radii) {
        std::vector<float> positions_vec(positions.data(), positions.data() + positions.size());
        std::vector<float> radii_vec(radii.data(), radii.data() + radii.size());
        auto surface_indices = OpenVDBSurfaceSplitter::split_surface_sdf(positions_vec, radii_vec);
        std::vector<uint8_t> is_surface_tmp(surface_indices.is_surface.begin(), surface_indices.is_surface.end());
        return {
            PyBoolArr(is_surface_tmp.size(), (const bool*)is_surface_tmp.data()),
            surface_indices.info_msg
        };
    }
    PySurfaceIndices split_surface_count(const PyFloatArr &positions) {
        std::vector<float> positions_vec(positions.data(), positions.data() + positions.size());
        auto surface_indices = OpenVDBSurfaceSplitter::split_surface_count(positions_vec);
        std::vector<uint8_t> is_surface_tmp(surface_indices.is_surface.begin(), surface_indices.is_surface.end());
        return {
            PyBoolArr(is_surface_tmp.size(), (const bool*)is_surface_tmp.data()),
            surface_indices.info_msg
        };
    }
};


PYBIND11_MODULE(ParticleMesherPy, m) {
    m.doc() = "Python binding of ParticleMesher";
    py::class_<PyMeshConstructor>(m, "MeshConstructor")
        .def(py::init<OpenVDBMeshConstructor::Config>(),
            py::arg("config")
        )
        .def("construct", &PyMeshConstructor::construct,
            py::arg("positions"),
            py::arg("radii")
        );
    py::class_<OpenVDBMeshConstructor::Config>(m, "MeshConstructorConfig")
        .def(py::init<float, float, float, float>(),
            py::arg("particle_radius"),
            py::arg("voxel_scale"),
            py::arg("isovalue") = 0.f,
            py::arg("adaptivity") = 0.f
        );
    py::class_<PyConstructMesh>(m, "ConstructMesh")
        .def_readonly("vertices", &PyConstructMesh::vertices)
        .def_readonly("triangles", &PyConstructMesh::triangles)
        .def_readonly("info_msg", &PyConstructMesh::info_msg);
    py::class_<PyFoamGenerator>(m, "FoamGenerator")
        .def(py::init<OpenVDBFoamGenerator::Config, std::string_view>(),
            py::arg("config"),
            py::arg("object_id")        // just a reminder
        )
        .def("generate_foams", &PyFoamGenerator::generate_foams,
            py::arg("positions"),
            py::arg("velocities")
        );
    py::class_<OpenVDBFoamGenerator::Config>(m, "FoamGeneratorConfig")
        .def(py::init<float, float, float, StdFloatArr, StdFloatArr, StdFloatArr,
                      StdFloatArr, StdFloatArr, StdFloatArr, StdFloatArr,
                      float, int, int, int, int,
                      float, float, float, float, float, float,
                      float, float, float, float>(),
            py::arg("particle_radius"),
            py::arg("voxel_scale"),
            py::arg("time_step"),
            py::arg("lower_bound"),
            py::arg("upper_bound"),
            py::arg("gravity"),
            py::arg("lim_ta") = py::make_tuple(0.32f, 3.2f),
            py::arg("lim_wc") = py::make_tuple(0.128f, 1.28f),
            py::arg("lim_ke") = py::make_tuple(0.5f, 5.0f),
            py::arg("lim_life") = py::make_tuple(2.f, 5.f),
            // py::arg("ta_min") = 5.f, py::arg("ta_max") = 20.f,
            // py::arg("wc_min") = 2.f, py::arg("wc_max") = 8.f,
            // py::arg("ke_min") = 5.f, py::arg("ke_max") = 50.f,
            // py::arg("life_min") = 2.f, py::arg("life_max") = 5.f,
            py::arg("support_scale") = 4.f,
            py::arg("surface_neighbor_max") = 20,
            py::arg("generate_neighbor_min") = 15, 
            py::arg("foam_neighbor_min") = 6,
            py::arg("foam_neighbor_max") = 20,
            py::arg("k_ta") = 10.f,
            py::arg("k_wc") = 10.f,
            // py::arg("k_ta") = 4000.f,
            // py::arg("k_wc") = 50000.f,
            py::arg("k_bo") = 2.f,
            py::arg("k_dr") = 0.8f,
            py::arg("k_ad") = 0.99f,
            py::arg("k_foam") = 1000.f,
            py::arg("spray_decay") = 2.f,
            py::arg("foam_decay") = 1.f,
            py::arg("bubble_decay") = 0.f,
            py::arg("foam_density") = 1000.f
        );
    py::class_<PyFoamSpheres>(m, "FoamSpheres")
        .def_readonly("positions", &PyFoamSpheres::positions)
        .def_readonly("info_msg", &PyFoamSpheres::info_msg);
    py::class_<PySurfaceSplitter>(m, "SurfaceSplitter")
        .def(py::init<OpenVDBSurfaceSplitter::Config>(),
            py::arg("config")
        )
        .def("split_surface_sdf", &PySurfaceSplitter::split_surface_sdf,
            py::arg("positions"),
            py::arg("radii")
        )
        .def("split_surface_count", &PySurfaceSplitter::split_surface_count,
            py::arg("positions")
        );
    py::class_<OpenVDBSurfaceSplitter::Config>(m, "SurfaceSplitterConfig")
        .def(py::init<float, float, float, float, int>(),
            py::arg("particle_radius"),
            py::arg("voxel_scale"),
            py::arg("support_scale") = 4.f,
            py::arg("half_width") = 3.f,
            py::arg("surface_neighbor_max") = 20
        );
    py::class_<PySurfaceIndices>(m, "SurfaceIndices")
        .def_readonly("is_surface", &PySurfaceIndices::is_surface)
        .def_readonly("info_msg", &PySurfaceIndices::info_msg);
}