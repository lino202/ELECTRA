# README #

# Description

ELECTRA is a command-line software and library for solving PDEs for cardiac electrophysiology simulation using Finite Element and meshless methods. 

## Features
- Cardiac Electrophysiology Simulator
- Finite Element Method
- Meshless Methods
- Non-distributed memory
- Ensight binary results
- Cardiac conduction system
- Monodomain solver

# Installation 

You can use ELECTRA by building from source (Linux) or by using the Docker image (Windows). 

## Build form source 

We assume you are working on a Unix based system.

You will need to install cmake higher than 3.14.0, and other dependencies like libboost. You can perfectly follow the Dockerfile for the dependencies and others posterior steps. 

After you have minimal dependencies, you need to download this repository and go to its location in your file system. So we suppose you are in [/any/path]/ELECTRA. 

Now, You need to build and install ELECTRA dependencies which are all under the deps folder (some are collect by cmake on the fly), so first make the bash-local environment variable

```
export ELECTRA_DEPS_DIR=[path-to-electra]/ELECTRA/deps/
```

Then we can install dependencies as follow. First we locate us in the ELECTRA source code folder
```
cd [/any/path]/ELECTRA
```

Install eigen
```
cd ${ELECTRA_DEPS_DIR}/eigen \
    && cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=RELEASE \
    -DCMAKE_INSTALL_PREFIX="${ELECTRA_DEPS_DIR}/eigen/install" \
    -DEIGEN_TEST_CXX11=ON -DEIGEN_TEST_OPENMP=ON -DEIGEN_TEST_SSE2=ON -DEIGEN_TEST_SSE3=ON -DEIGEN_TEST_SSE4_1=ON -DEIGEN_TEST_SSE4_2=ON \
    && cd build && make -j4 && make install
```
Install armadillo
```
cd ${ELECTRA_DEPS_DIR}/armadillo \
    && cmake -S. -Bbuild -DCMAKE_INSTALL_PREFIX="${ELECTRA_DEPS_DIR}/armadillo/install" -DCMAKE_BUILD_TYPE=RELEASE \
    && cd build && make -j4 && make install
```

Install CGAL
```
cd ${ELECTRA_DEPS_DIR}/cgal \
    && cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=RELEASE -DCGAL_HEADER_ONLY=ON \
    -DWITH_Eigen3=ON -DWITH_GMP=ON -DWITH_MPFR=ON -DCMAKE_INSTALL_PREFIX="${ELECTRA_DEPS_DIR}/cgal/install" \
    && cd build && make -j4 && make install
```

Install IMP
```
cd ${ELECTRA_DEPS_DIR}/IMP \
    && cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_APPS=OFF -DBUILD_DOC=OFF -DBUILD_TESTS=OFF \
    -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="${ELECTRA_DEPS_DIR}/IMP/install" \
    && cd build && make -j4 && make install
```

Install CLOUDEA
```
cd ${ELECTRA_DEPS_DIR}/CLOUDEA \
    && cmake -S. -Bbuild -DCMAKE_INSTALL_PREFIX="${ELECTRA_DEPS_DIR}/CLOUDEA/install" \
    -DCMAKE_BUILD_TYPE=RELEASE -DCLOUDEA_USE_CGAL=ON -DBUILD_APPS=OFF -DBUILD_DOC=OFF -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_PREFIX_PATH="${ELECTRA_DEPS_DIR}/IMP/install;${ELECTRA_DEPS_DIR}/armadillo/install;${ELECTRA_DEPS_DIR}/cgal/install" \
    && cd build && make -j4 && make install
```

Then you need to build ELECTRA
```
cd /home/ELECTRA \
    && cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=RELEASE \
    -DBUILD_DOC=OFF -DELECTRA_WITH_CUDA=OFF -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_PREFIX_PATH="$ELECTRA_DEPS_DIR/IMP/install;$ELECTRA_DEPS_DIR/CLOUDEA/install;$ELECTRA_DEPS_DIR/armadillo/install;$ELECTRA_DEPS_DIR/cgal/install" \
    && cd build && make -j4 && cd ../
```

Under [/any/path]/ELECTRA/build/bin you can find the ElectraSim application

## Docker 

You can download the latest docker image released from the release section [Last Release](https://github.com/lino202/ELECTRA/releases/latest).
Once you have the image and Docker is installed in your machine you should:

```
docker load -i [path-to-dokcer-image]/electra-docker_[version-tag].tar.gz
```

Then, create the container (without any mounted volumen) 
```
docker run -it --name [desired-nname] electra-docker:[version-tag]
```

Or with a mounted container for sharing files and folders between the host and your container
```
docker run -it --name [desired-nname] -v [host-folder-path]:[container-path] electra-docker:[version-tag]
```

Then you will be inside the container. Here you can see ElectraSim app under /home/ELECTRA/build/bin, the entire source code should be under /home/ELECTRA and you should be able to use the command ElectraSim to call the binary under /home/ELECTRA/build/bin.
If it is a new day and you can do the following for opening your ELECTRA Docker container:
```
docker start [container-name]
```

After starting, you can enter the container command line iteratively doing:
```
docker exec -it [container-name] /bin/bash
```

In docker you can use the application and moreover you have all the dependencies and source code for building, debugging and contributing to ELECTRA. The build can be done as in the build section for unix based system

```
cd /home/ELECTRA \
    && cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=RELEASE \
    -DBUILD_DOC=OFF -DELECTRA_WITH_CUDA=OFF -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_PREFIX_PATH="$ELECTRA_DEPS_DIR/IMP/install;$ELECTRA_DEPS_DIR/CLOUDEA/install;$ELECTRA_DEPS_DIR/armadillo/install;$ELECTRA_DEPS_DIR/cgal/install" \
    && cd build && make -j4
```

If you want to debug change in the above command RELEASE for DEBUG


### Docker image generation

If you are on your local machine or in the docker container you can generate a docker image file to share your ELECTRA version with others.
Execute 

```
cd [path-to-electra]/ELECTRA
./dockerGenImage.sh
```

This sript dockerGenImage would ask you to input the version tag for your docker image which should be equal to the one in ./CMakeLists.txt VERSION field. Afterwards, the script execute the DockerFile instructions for creating the container, copying the folder files and building ELECTRA with its dependencies. If the script finish without errors, you should see a new folder images where your electra-docker_[version].tar.gz image should be.

Now you can install and use your image as described in [Docker section](#docker)

# Usage

The application of ElectraSim needs a .json as an input for setting up the simulation. Several fields are needed. See [Examples](https://github.com/lino202/ELECTRA/tree/main/examples) for getting started. 

An explanation and options of the .json files will be here anytime soon.


# Contributors

This software was mainly developed by [@KMountris](https://github.com/KMountris)