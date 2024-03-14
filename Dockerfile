FROM ubuntu:22.04

# Copy src files from host to image.
COPY ./ /home/ELECTRA/

# Install basic libs
RUN apt-get update && apt-get upgrade -y \
    && apt-get install -y apt-utils \
    git \
    g++ \
    make \
    wget \
    nano \
    xz-utils \
    mlocate \
    libssl-dev \
    libncurses5 \
    libncurses5-dev \
    libomp-dev \
    libgmp-dev \
    libmpfr-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libboost-thread-dev \
    libboost-program-options-dev \
    libopenblas-base \
    libopenblas-dev \
    liblapack-dev \
    libatlas-base-dev \
    libsuitesparse-dev \
    libsuperlu-dev

# Download and install cmake
RUN cd /home \
    && wget https://github.com/Kitware/CMake/releases/download/v3.16.5/cmake-3.16.5.tar.gz \
    && tar xvf cmake-3.16.5.tar.gz \
    && rm cmake-3.16.5.tar.gz \
    && mv cmake-3.16.5 cmake \
    && mkdir /home/cmake/build \
    && cd /home/cmake/build \
    && ../bootstrap --parallel=4 \
    && make -j4 && make install \
    && cd ../.. \
    && rm -rf cmake

# Add Deps env variable
ENV ELECTRA_DEPS_DIR="/home/ELECTRA/deps"
RUN cd /home/ELECTRA


# Build the dependencies which are under ELECTRA_DEPS_DIR already for simplicity ----------------------------------------------------

#Eigen
RUN cd ${ELECTRA_DEPS_DIR}/eigen \
    && cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=RELEASE \
    -DCMAKE_INSTALL_PREFIX="${ELECTRA_DEPS_DIR}/eigen/install" \
    -DEIGEN_TEST_CXX11=ON -DEIGEN_TEST_OPENMP=ON -DEIGEN_TEST_SSE2=ON -DEIGEN_TEST_SSE3=ON -DEIGEN_TEST_SSE4_1=ON -DEIGEN_TEST_SSE4_2=ON \
    && cd build && make -j4 && make install

#Armadillo
RUN cd ${ELECTRA_DEPS_DIR}/armadillo \
    && cmake -S. -Bbuild -DCMAKE_INSTALL_PREFIX="${ELECTRA_DEPS_DIR}/armadillo/install" -DCMAKE_BUILD_TYPE=RELEASE \
    && cd build && make -j4 && make install

#CGAL
RUN cd ${ELECTRA_DEPS_DIR}/cgal \
    && cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=RELEASE -DCGAL_HEADER_ONLY=ON \
    -DWITH_Eigen3=ON -DWITH_GMP=ON -DWITH_MPFR=ON -DCMAKE_INSTALL_PREFIX="${ELECTRA_DEPS_DIR}/cgal/install" \
    && cd build && make -j4 && make install


#IMP
RUN cd ${ELECTRA_DEPS_DIR}/IMP \
    && cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_APPS=OFF -DBUILD_DOC=OFF -DBUILD_TESTS=OFF \
    -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="${ELECTRA_DEPS_DIR}/IMP/install" \
    && cd build && make -j4 && make install

#CLOUDEA
RUN cd ${ELECTRA_DEPS_DIR}/CLOUDEA \
    && cmake -S. -Bbuild -DCMAKE_INSTALL_PREFIX="${ELECTRA_DEPS_DIR}/CLOUDEA/install" \
    -DCMAKE_BUILD_TYPE=RELEASE -DCLOUDEA_USE_CGAL=ON -DBUILD_APPS=OFF -DBUILD_DOC=OFF -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_PREFIX_PATH="${ELECTRA_DEPS_DIR}/IMP/install;${ELECTRA_DEPS_DIR}/armadillo/install;${ELECTRA_DEPS_DIR}/cgal/install" \
    && cd build && make -j4 && make install

# BUILD ELECTRA -----------------------------------------------------------------------------------------------
RUN cd /home/ELECTRA \
    && cmake -S. -Bbuild -DCMAKE_BUILD_TYPE=RELEASE \
    -DBUILD_DOC=OFF -DELECTRA_WITH_CUDA=OFF -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_PREFIX_PATH="$ELECTRA_DEPS_DIR/IMP/install;$ELECTRA_DEPS_DIR/CLOUDEA/install;$ELECTRA_DEPS_DIR/armadillo/install;$ELECTRA_DEPS_DIR/cgal/install" \
    && cd build && make -j4 && cd ../../..

# Add shortcut to ELECTRA applications
RUN echo "alias ElectraSim='/home/ELECTRA/build/bin/ElectraSim'" >> ~/.bashrc
RUN echo "alias ElectraSim='/home/ELECTRA/build/bin/ElectraCell'" >> ~/.bashrc

# Finished successfully
# constructing ELECTRA image.
CMD ["bash"]
