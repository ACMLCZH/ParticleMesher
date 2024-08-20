add_library(particle_mesher SHARED
    mesh_construct.cpp
    foam_generator.cpp
    surface_splitter.cpp
)
target_link_libraries(particle_mesher PRIVATE OpenVDB::openvdb)

if(DEFINED PYTHON_VERSIONS)
    foreach(version ${PYTHON_VERSIONS})
        set(PYBIND11_PYTHON_VERSION ${version})
        find_package(Python ${version} EXACT REQUIRED COMPONENTS Interpreter Development)
        pybind11_add_module(ParticleMesherPy py_particle_mesher.cpp)
        target_link_libraries(ParticleMesherPy PRIVATE particle_mesher)
    endforeach()
else()
    find_package(Python REQUIRED COMPONENTS Interpreter Development)
    message("Default Python ${Python_VERSION} in environment was found")
    pybind11_add_module(ParticleMesherPy py_particle_mesher.cpp)
    target_link_libraries(ParticleMesherPy PRIVATE particle_mesher)
    # set_target_properties(ParticleMesherPy PROPERTIES POSITION_INDEPENDENT_CODE ON)
    # set_target_properties(ParticleMesherPy PROPERTIES LINK_FLAGS "-static")
endif()
