FROM quay.io/pypa/manylinux2014_x86_64

ENV BOOST_MAJOR=1
ENV BOOST_MINOR=73
ENV BOOST_PATCH=0
ENV GMP_VERSION=6.3.0
ENV MPFR_VERSION = 4.2.1
ENV METIS_VERSION=5.1.1
ENV LAPACK_VERSION=3.12.0
ENV MUMPS_VERSION=3.0.5
ENV IPOPT_VERSION=3.14.14
ENV TBB_VERSION=2021.10.0
ENV CRITERION_VERSION=2.4.2
ENV PAPILO_REVISION=v2.2.0
ENV SOPLEX_REVISION=release-700

RUN yum update -y && \
    yum install -y \
    gcc \
    gcc-c++ \
    libgfortran \
    make \
    cmake \
    git \
    wget \
    glibc-static && \
    yum clean all

# lapack and blas
WORKDIR /home/builds

RUN wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v$LAPACK_VERSION.tar.gz && \
    tar xfz v$LAPACK_VERSION.tar.gz && \
    cd /home/builds/lapack-$LAPACK_VERSION && \
    mkdir build && \
    cd /home/builds/lapack-$LAPACK_VERSION/build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DCMAKE_POSITION_INDEPENDENT_CODE=ON && \
    make -j && \
    mv lib/libblas.a /usr/local/lib && \
    mv lib/liblapack.a /usr/local/lib && \
    rm -rf /home/builds/*

# gmp
RUN wget https://gmplib.org/download/gmp/gmp-$GMP_VERSION.tar.xz && \
    tar xf gmp-$GMP_VERSION.tar.xz && \
    cd /home/builds/gmp-$GMP_VERSION && \
    ./configure --with-pic --enable-cxx --disable-shared && \
    make install -j && \
    rm -rf /home/builds/*

# gklib and metis
RUN wget https://github.com/KarypisLab/GKlib/archive/refs/tags/METIS-v$METIS_VERSION-DistDGL-0.5.tar.gz && \
    tar xfz METIS-v$METIS_VERSION-DistDGL-0.5.tar.gz && \
    cd /home/builds/GKlib-METIS-v$METIS_VERSION-DistDGL-0.5 && \
    make config prefix=/home/builds/GKlib-METIS-v$METIS_VERSION-DistDGL-0.5 && \
    make -j && \
    make install && \
    cd /home/builds && \
    wget https://github.com/KarypisLab/METIS/archive/refs/tags/v$METIS_VERSION-DistDGL-v0.5.tar.gz && \
    tar xfz v$METIS_VERSION-DistDGL-v0.5.tar.gz && \
    cd /home/builds/METIS-$METIS_VERSION-DistDGL-v0.5 && \
    make config prefix=/usr/local gklib_path=/home/builds/GKlib-METIS-v$METIS_VERSION-DistDGL-0.5 && \
    make && \
    make install && \
    rm -rf /home/builds/*

# mumps
RUN wget https://github.com/coin-or-tools/ThirdParty-Mumps/archive/refs/tags/releases/$MUMPS_VERSION.zip && \
    unzip $MUMPS_VERSION.zip && \
    cd /home/builds/ThirdParty-Mumps-releases-$MUMPS_VERSION && \
    ./get.Mumps && ./configure --disable-shared --enable-static --with-lapack-lflags="-llapack -lblas -lgfortran -lquadmath -lm" && \
    make -j && \
    make install && \
    rm -rf /home/builds/*

# ipopt
RUN wget https://github.com/coin-or/Ipopt/archive/releases/$IPOPT_VERSION.zip && \
    unzip $IPOPT_VERSION.zip && \
    cd /home/builds/Ipopt-releases-$IPOPT_VERSION && \
    ./configure --enable-static --disable-shared --with-lapack="-llapack -lblas -lgfortran -lquadmath -lm" && \
    make -j && make install && \
    rm -rf /home/builds/*

# boost
RUN echo "Download boost from https://boostorg.jfrog.io/artifactory/main/release/${BOOST_MAJOR}.${BOOST_MINOR}.${BOOST_PATCH}/source/boost_${BOOST_MAJOR}_${BOOST_MINOR}_${BOOST_PATCH}.tar.gz" && \
    wget https://boostorg.jfrog.io/artifactory/main/release/${BOOST_MAJOR}.${BOOST_MINOR}.${BOOST_PATCH}/source/boost_${BOOST_MAJOR}_${BOOST_MINOR}_${BOOST_PATCH}.tar.gz && \
    tar xfz boost_${BOOST_MAJOR}_${BOOST_MINOR}_${BOOST_PATCH}.tar.gz && \
    cd /home/builds/boost_${BOOST_MAJOR}_${BOOST_MINOR}_${BOOST_PATCH} && \
    ./bootstrap.sh --with-libraries=program_options,serialization && \
    ./b2 -j3 install && \
    rm -rf /home/builds/*

# soplex
RUN git clone https://github.com/scipopt/soplex.git && \
    cd /home/builds/soplex && \
    git checkout $SOPLEX_REVISION && \
    mkdir /home/builds/soplex/build && \
    cd /home/builds/soplex/build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make -j && \
    make install && \
    rm -rf /home/builds/*

# tbb
RUN wget https://github.com/oneapi-src/oneTBB/archive/refs/tags/v$TBB_VERSION.tar.gz && \
    tar xvf v$TBB_VERSION.tar.gz && \
    mkdir /home/builds/oneTBB-$TBB_VERSION/build && \
    cd /home/builds/oneTBB-$TBB_VERSION/build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release -DTBB_TEST=OFF -DTBB_EXAMPLES=OFF -DTBB4PY_BUILD=OFF -DBUILD_SHARED_LIBS=off && \
    make -j3 && make install && \
    rm -rf /home/builds/*

# papilo
RUN git clone https://github.com/scipopt/papilo.git && \
    cd /home/builds/papilo && \
    git checkout $PAPILO_REVISION && \
    mkdir /home/builds/papilo/build && \
    cd /home/builds/papilo/build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release -DQUADMATH=off -DBLA_STATIC=ON -DBLA_VENDOR=OpenBLAS && \
    make -j5 && \
    make install && \
    rm -rf /home/builds/*

# criterion
RUN wget https://github.com/Snaipe/Criterion/releases/download/v2.4.2/criterion-$CRITERION_VERSION-linux-x86_64.tar.xz && \
    tar xf criterion-$CRITERION_VERSION-linux-x86_64.tar.xz && \
    cp -r criterion-$CRITERION_VERSION/* /usr/local && \
    rm -rf /home/builds/*
