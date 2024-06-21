from setuptools import setup, Extension

ext_module = Extension(
    name='ParticleMesherPy',  # Name of the Python module
    sources=[],  # No source files needed
    library_dirs=['./ParticleMesherPy'],
    libraries=[          # Names of the prebuilt libraries
        'ParticleMesherPy',
        'boost_iostreams',
        'blosc',
        'openvdb',
        'tbb',
        'particle_mesher',
    ],  
)

setup(
    name='ParticleMesherPy',
    version='0.1',
    python_requires='>=3.9',
    ext_modules=[ext_module],
)

# setup(
#     name='ParticleMesherPy',
#     version='0.1',
#     package_data={'': ['*.so']},
#     install_requires=[],  # Add any dependencies here
# )