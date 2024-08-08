# ParticleMesher
 
### Dependencies
- OpenVDB
    Follow https://www.openvdb.org/documentation/doxygen/build.html to install OpenVDB (Building Standalone).
    You should first clone OpenVDB repository a place you wish to install:
    ```
    git clone https://github.com/AcademySoftwareFoundation/openvdb.git
    ```
    Then install dependencies that building OpenVDB requires:
    ```
    sudo apt-get install cmake                   # CMake
    sudo apt-get install libtbb-dev              # TBB
    sudo apt-get install zlibc                   # zlib
    sudo apt-get install libboost-iostreams-dev  # Boost::iostream
    sudo apt-get install libblosc-dev            # Blosc
    ```
    Then build and install OpenVDB:
    ```
    cd openvdb
    mkdir build
    cd build
    cmake ../
    cmake --build .
    cmake --build . --target install
    ```
    If there's complaint about boost version: https://stackoverflow.com/a/75893242

    After installation you will find your OpenVDB module path. It is often `/usr/local/lib/cmake/OpenVDB`

    ```
    libboost_iostreams.so.1.73.0
    libblosc.so.1
    libopenvdb.so.11.0
    libtbb.so.12
    libparticle_mesher.so
    ```

### Compile
- Build ParticleMesher and its python binding:
    ```
    cd ./ext/ParticleMesher
    cmake -S . -B build -D CMAKE_BUILD_TYPE=Release -D PYTHON_VERSIONS=3.9 -D OPENVDB_MODULE_PATH=[OpenVDB module path]
    cmake --build build -j $(nproc)
    ```

### Add to Path
- add to your `.bashrc`
    ```
    echo "export PMP_PATH=${PWD}/build/bin"
    source ~/.bashrc
    ```