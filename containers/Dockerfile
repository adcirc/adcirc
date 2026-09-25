# syntax=docker/dockerfile:1

# ADCIRC is compiled in the devel base image and shipped in the runtime base
# image, which has the same libraries but no compilers. Both come from
# https://github.com/adcirc/adcirc-ci-docker and must be the same version.
ARG BASE_VERSION=2026.1.0
ARG BUILD_IMAGE=adcircorg/adcirc-base:${BASE_VERSION}-devel
ARG RUNTIME_IMAGE=adcircorg/adcirc-base:${BASE_VERSION}-runtime

FROM ${BUILD_IMAGE} AS builder

ARG TARGETARCH

SHELL [ "/bin/bash", "-o", "pipefail", "-c" ]

# Root is needed to write to /opt, and this stage is never published
# hadolint ignore=DL3002
USER 0

# The source tree lives at the same path in the builder and final images so
# that file paths embedded in debug info and tracebacks resolve in the container
COPY src /opt/adcirc-source/src
COPY prep /opt/adcirc-source/prep
COPY scripts /opt/adcirc-source/scripts
COPY util /opt/adcirc-source/util
COPY wind /opt/adcirc-source/wind
COPY cmake /opt/adcirc-source/cmake
COPY thirdparty /opt/adcirc-source/thirdparty
COPY CMakeLists.txt /opt/adcirc-source/CMakeLists.txt
COPY .git /opt/adcirc-source/.git

WORKDIR /opt/adcirc-source

RUN <<EOF
set -eo pipefail
source /opt/spack-environment/activate.sh

# x86 images use the Intel compilers, ARM images use GCC
case "${TARGETARCH}" in
  amd64)
    compiler_args=(
      -DCMAKE_C_COMPILER=icx
      -DCMAKE_CXX_COMPILER=icpx
      -DCMAKE_Fortran_COMPILER=ifx
      -DCMAKE_Fortran_FLAGS_RELEASE="-O2 -g -traceback -fp-model=precise"
      -DCMAKE_C_FLAGS_RELEASE="-O2 -g -DNDEBUG -fp-model=precise"
      -DCMAKE_CXX_FLAGS_RELEASE="-O2 -g -DNDEBUG -fp-model=precise" ) ;;
  arm64)
    compiler_args=(
      -DCMAKE_C_COMPILER=gcc
      -DCMAKE_CXX_COMPILER=g++
      -DCMAKE_Fortran_COMPILER=gfortran
      -DCMAKE_Fortran_FLAGS_RELEASE="-O2 -g -fbacktrace"
      -DCMAKE_C_FLAGS_RELEASE="-O2 -g -DNDEBUG"
      -DCMAKE_CXX_FLAGS_RELEASE="-O2 -g -DNDEBUG" ) ;;
  *)
    echo "Unsupported architecture: ${TARGETARCH}" >&2
    exit 1 ;;
esac

cmake -S . -B build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/opt/adcirc \
    "${compiler_args[@]}" \
    -DBUILD_ADCIRC=ON \
    -DBUILD_PADCIRC=ON \
    -DBUILD_ADCSWAN=ON \
    -DBUILD_PADCSWAN=ON \
    -DBUILD_ADCPREP=ON \
    -DBUILD_ASWIP=ON \
    -DBUILD_UTILITIES=ON \
    -DENABLE_OUTPUT_NETCDF=ON \
    -DENABLE_OUTPUT_XDMF=ON \
    -DENABLE_DATETIME=ON \
    -DENABLE_GRIB2=ON \
    -DNETCDF_F90_ROOT="${NETCDF_FORTRAN_HOME}" \
    -DNETCDF_DIR="${NETCDFHOME}" \
    -DXDMFHOME="${XDMFHOME}"
cmake --build build --parallel
cmake --install build
EOF

FROM ${RUNTIME_IMAGE} AS adcirc

SHELL [ "/bin/bash", "-o", "pipefail", "-c" ]

USER 0

COPY --from=builder /opt/adcirc /opt/adcirc
COPY --from=builder /opt/adcirc-source/src /opt/adcirc-source/src
COPY --from=builder /opt/adcirc-source/prep /opt/adcirc-source/prep
COPY --from=builder /opt/adcirc-source/wind /opt/adcirc-source/wind
COPY --from=builder /opt/adcirc-source/util /opt/adcirc-source/util
COPY --from=builder /opt/adcirc-source/thirdparty/swan /opt/adcirc-source/thirdparty/swan

# activate.sh sets PATH outright, so /opt/adcirc/bin is added there for shells
# started through the entrypoint as well as through ENV for everything else
# (e.g. singularity exec, docker run --entrypoint)
RUN echo "export PATH=/opt/adcirc/bin:\$PATH" >> /opt/spack-environment/activate.sh
ENV PATH=/opt/adcirc/bin:${PATH}

# The base image's non-root user (see adcirc-ci-docker for how it supports
# arbitrary UIDs)
USER 1001
WORKDIR /home/adcirc

# Fail the build rather than publish an image whose executables cannot start
RUN <<EOF
set -e
if ldd /opt/adcirc/bin/* | grep "not found"; then
  exit 1
fi
adcirc -v
padcirc -v
adcprep --help > /dev/null
aswip --help > /dev/null
EOF


LABEL org.opencontainers.image.title="ADCIRC"
LABEL org.opencontainers.image.description="ADCIRC coastal circulation and storm surge model"
LABEL org.opencontainers.image.authors="Zach Cobell <zcobell@gmail.com>"
LABEL org.opencontainers.image.source="https://github.com/adcirc/adcirc"
LABEL org.opencontainers.image.licenses="LicenseRef-ADCIRC"
LABEL io.k8s.description="ADCIRC Docker Container"
LABEL io.openshift.expose-services="None"

ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash" ]
