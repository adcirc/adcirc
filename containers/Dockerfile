FROM adcircorg/adcirc-base:2025.0.1 AS builder

ARG TARGETARCH
ARG TARGETPLATFORM

SHELL [ "/bin/bash", "-o", "pipefail", "-c" ]

USER root

RUN mkdir -p /home/adcirc/adcirc-build/build
COPY src /home/adcirc/adcirc-build/src
COPY prep /home/adcirc/adcirc-build/prep
COPY scripts /home/adcirc/adcirc-build/scripts
COPY util /home/adcirc/adcirc-build/util
COPY wind /home/adcirc/adcirc-build/wind
COPY cmake /home/adcirc/adcirc-build/cmake
COPY thirdparty /home/adcirc/adcirc-build/thirdparty
COPY CMakeLists.txt /home/adcirc/adcirc-build/CMakeLists.txt 
COPY .git /home/adcirc/adcirc-build/.git

RUN cd /home/adcirc/adcirc-build/build && \
    source /opt/spack-environment/activate.sh && \
    if [[ "x${TARGETARCH}" == "xamd64" ]]; then \
      cmake .. \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX=/opt/adcirc \
          -DCMAKE_C_COMPILER=icx \
          -DCMAKE_CXX_COMPILER=icpx \
          -DCMAKE_Fortran_COMPILER=ifx \
          -DBUILD_ADCIRC=ON \
          -DBUILD_PADCIRC=ON \
          -DBUILD_ADCSWAN=ON \
          -DBUILD_PADCSWAN=ON \
          -DBUILD_ADCPREP=ON \
          -DBUILD_ASWIP=ON \
          -DENABLE_OUTPUT_NETCDF=ON \
          -DENABLE_OUTPUT_XDMF=ON \
          -DENABLE_GRIB2=ON \
          -DCMAKE_Fortran_FLAGS_RELEASE="-O2 -fp-model=precise" \
          -DCMAKE_C_FLAGS_RELEASE="-O2 -DNDEBUG -fp-model=precise" \
          -DCMAKE_CXX_FLAGS_RELEASE="-O2 -DNDEBUG -fp-model=precise" \
          -DNETCDF_F90_ROOT=$NETCDF_FORTRAN_HOME \
          -DNETCDFHOME=$NETCDFHOME \
          -DXDMFHOME=$XDMFHOME \
          -G Ninja; \ 
    elif [[ "x${TARGETARCH}" == "xarm64" ]]; then \
      cmake .. \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_INSTALL_PREFIX=/opt/adcirc \
          -DCMAKE_C_COMPILER=gcc \
          -DCMAKE_CXX_COMPILER=g++ \
          -DCMAKE_Fortran_COMPILER=gfortran \
          -DBUILD_ADCIRC=ON \
          -DBUILD_PADCIRC=ON \
          -DBUILD_ADCSWAN=ON \
          -DBUILD_PADCSWAN=ON \
          -DBUILD_ADCPREP=ON \
          -DBUILD_ASWIP=ON \
          -DENABLE_OUTPUT_NETCDF=ON \
          -DENABLE_OUTPUT_XDMF=ON \
          -DENABLE_GRIB2=ON \
          -DCMAKE_Fortran_FLAGS_RELEASE="-O2" \
          -DCMAKE_C_FLAGS_RELEASE="-O2 -DNDEBUG" \
          -DCMAKE_CXX_FLAGS_RELEASE="-O2 -DNDEBUG" \
          -DNETCDF_F90_ROOT=$NETCDF_FORTRAN_HOME \
          -DNETCDFHOME=$NETCDFHOME \
          -DXDMFHOME=$XDMFHOME \
          -G Ninja; \
    fi; \
    ninja && \
    ninja install 

FROM adcircorg/adcirc-base:2025.0.1 AS adcirc

USER root 
COPY --from=builder /opt/adcirc /opt/adcirc 
RUN echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/opt/views/view/lib64" >> /opt/spack-environment/activate.sh
RUN echo "export PATH=\$PATH:/opt/adcirc/bin" >> /opt/spack-environment/activate.sh

USER adcirc
WORKDIR /home/adcirc

LABEL "maintainer"="Zach Cobell <zcobell@gmail.com>"
LABEL "io.k8s.description"="ADCIRC Docker Container"
LABEL "io.openshift.expose-services"="None"

ENTRYPOINT [ "/entrypoint.sh" ]
CMD [ "/bin/bash" ]

