from setuptools import setup, find_packages

setup(
    name='ParticleMesherPy',
    version='0.1',
    python_requires='>=3.9',
    packages=find_packages(),
    package_data={'': ['*.so']},
    install_requires=[],  # Add any dependencies here
)