FROM zcobell/adcirc-ci-container:v56 AS builder

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

RUN cd /home/adcirc/adcirc-build/build && \
    source /opt/spack-environment/activate.sh && \
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=/opt/adcirc \
        -DCMAKE_C_COMPILER=icx \
        -DCMAKE_CXX_COMPILER=icx \
        -DCMAKE_Fortran_COMPILER=ifx \
        -DBUILD_ADCIRC=ON \
        -DBUILD_PADCIRC=ON \
        -DBUILD_ADCSWAN=ON \
        -DBUILD_PADCSWAN=ON \
        -DBUILD_ADCPREP=ON \
        -DBUILD_ASWIP=ON \
        -DENABLE_OUTPUT_NETCDF=ON \
        -DENABLE_OUTPUT_XDMF=ON \
        -DENABLE_DATETIME=ON \
        -DENABLE_GRIB2=ON \
        -DCMAKE_Fortran_FLAGS_RELEASE="-O2 -fp-model=precise" \
        -DCMAKE_C_FLAGS_RELEASE="-O2 -DNDEBUG -fp-model=precise" \
        -DCMAKE_CXX_FLAGS_RELEASE="-O2 -DNDEBUG -fp-model=precise" \
        -DNETCDF_F90_ROOT=$NETCDF_FORTRAN_HOME \
        -DNETCDFHOME=$NETCDFHOME \
        -DXDMFHOME=$XDMFHOME \
        -G Ninja && \
        ninja && \
        ninja install 

FROM zcobell/adcirc-ci-container:v56 AS adcirc

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

