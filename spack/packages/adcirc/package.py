# ----------------------------------------------------------------------------#
#
#     spack install adcirc
#
# You can edit this file again by typing:
#
#     spack edit adcirc
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------#
#
# ----------------------------------------------------------------------------#
#                                                                             #
#                              ADCIRC                                         #
#                                                                             #
#    A (PARALLEL) ADVANCED CIRCULATION MODEL FOR SHELVES, COASTAL SEAS        #
#                         AND ESTUARIES                                       #
#                                                                             #
#                                                                             #
#                          DEVELOPED BY:                                      #
#                                                                             #
#                      DR. R.A. LUETTICH, JR                                  #
#                                                                             #
#             UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL                     #
#                   INSTITUTE OF MARINE SCIENCES                              #
#                                                                             #
#                        DR. J.J. WESTERINK                                   #
#                                                                             #
#          DEPARTMENT OF CIVIL ENGINEERING AND GEOLOGICAL SCIENCES            #
#                     UNIVERSITY OF NOTRE DAME                                #
#                       NOTRE DAME, IN 46556                                  #
#                                                                             #
#                                                                             #
#        MAJOR FUNDING FOR THE DEVELOPMENT OF ADCIRC WAS PROVIDED BY          #
#                                                                             #
#                       DEPARTMENT OF THE ARMY                                #
#                    WATERWAYS EXPERIMENT STATION                             #
#                 COASTAL ENGINEERING RESEARCH CENTER                         #
#                        3909 HALLS FERRY RD                                  #
#                      VICKSBURG, MI 39180-6199                               #
#                                                                             #
# ----------------------------------------------------------------------------#
#                                                                             #
#          THE ADCIRC SOURCE CODE IS COPYRIGHTED, 1994-2022 BY:               #
#                                                                             #
#                 R.A. LUETTICH, JR AND J.J. WESTERINK                        #
#                                                                             #
# ----------------------------------------------------------------------------#

from spack import *


class Adcirc(CMakePackage):
    """ADCIRC: The ADvanced CIRCulation model for simulation of time dependent, free surface circulation and transport in two and three dimensions"""

    # ...Metadata
    homepage = "https://www.adcirc.org"
    maintainers = ["zcobell"]

    # ...Package location and sample archive name
    git = "https://github.com/adcirc/adcirc.git"

    # ...ADCIRC versions
    version("55.01", commit="58ddcb5f63a390fa7b53ce791c9f756f5ea0bc4b", preferred=True)
    version("55.00", commit="b804a461db68ec73b669d1d82469f5a7a32d66a5", deprecated=True)
    version("54.02", commit="29bb38ca684647eef3bfa1152785b1e762137f37", deprecated=True)
    version("54.01", commit="87e9cfc07fc04d5ff2b29ffe36cefadbae4024ac", deprecated=True)
    version("54.00", commit="fcac69501a050f26d8129bdc65bf802709c406e8", deprecated=True)

    version("master", branch="master")
    version("develop", branch="development")

    # ...Build variants
    variant(
        "netcdf", default=True, description="Build with with netCDF4 format enabled"
    )
    variant(
        "grib",
        default=False,
        description="Builds the model with GRIB format enabled",
        when=("@55:+netcdf"),
    )
    variant("mpi", default=True, description="Builds the parallel executables")
    variant(
        "swan", default=False, description="Builds the tightly coupled SWAN wave model"
    )
    variant(
        "libshared", default=False, description="Builds libadcirc as a shared library"
    )
    variant(
        "libstatic", default=False, description="Builds libadcirc as a static library"
    )
    variant("aswip", default=False, description="Builds aswip")
    variant(
        "utilities", default=False, description="Builds the adcirc utilities package"
    )

    # ...Dependencies
    depends_on("cmake", type="build")
    depends_on("perl", type="build", when="+swan")
    depends_on("hdf5~threadsafe~mpi", when="+netcdf", type=("build", "link"))
    depends_on("netcdf-c~mpi", when="+netcdf", type=("build", "link"))
    depends_on("netcdf-fortran", when="+netcdf", type=("build", "link"))
    depends_on("mpi", when="+mpi", type=("build", "link"))


    def cmake_args(self):
        args = []

        if "+netcdf" in self.spec:
            args.append(self.define("ENABLE_OUTPUT_NETCDF", True))
            args.append(
                self.define("NETCDF_F90_ROOT", self.spec["netcdf-fortran"].prefix)
            )

        if "+grib" in self.spec:
            args.append(self.define("ENABLE_GRIB2", True))
            args.append(self.define("ENABLE_DATETIME", True))

        if "+swan" in self.spec:
            args.append(self.define("BUILD_ADCSWAN", True))

        if "+mpi" in self.spec:
            args.append(self.define("BUILD_PADCIRC", True))
            args.append(self.define("BUILD_ADCPREP", True))
            if "+swan" in self.spec:
                args.append(self.define("BUILD_PADCSWAN", True))

        if "+libshared" in self.spec:
            args.append(self.define("BUILD_LIBADCIRC_SHARED", True))

        if "+libstatic" in self.spec:
            args.append(self.define("BUILD_LIBADCIRC_STATIC", True))

        if "+aswip" in self.spec:
            args.append(self.define("BUILD_ASWIP", True))

        if "+utilities" in self.spec:
            args.append(self.define("BUILD_UTILITIES", True))

        # ...The gcc10+ fix. ADCIRC is generally ok without this but certain
        #   mpi flavors will cause issues
        if self.spec.satisfies("%gcc@10:"):
            args.append(self.define("CMAKE_Fortran_FLAGS", "-fallow-argument-mismatch"))

        args.append(self.define("BUILD_ADCIRC", True))
        return args
