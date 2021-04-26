#!/usr/bin/env python3

class Configuration:
    def __init__(self,name,cc,cxx,fc,netcdfhome,modules,copt,cxxopt,fcopt,cdbg,cxxdbg,fcdbg):
        self.__name = name
        self.__cc = cc
        self.__cxx = cxx
        self.__fc = fc
        self.__netcdfhome = netcdfhome
        self.__modules = modules
        self.__copt = copt
        self.__cxxopt = cxxopt
        self.__fcopt = fcopt
        self.__cdbg = cdbg
        self.__cxxdbg = cxxdbg
        self.__fcdbg = fcdbg

    def name(self):
        return self.__name
    def cc(self): 
        return self.__cc
    def cxx(self):
        return self.__cxx
    def fc(self):
        return self.__fc
    def netcdfhome(self):
        return self.__netcdfhome
    def modules(self):
        return self.__modules
    def copt(self):
        return self.__copt
    def cxxopt(self):
        return self.__cxxopt
    def fcopt(self):
        return self.__fcopt
    def cdbg(self):
        return self.__cdbg
    def cxxdbg(self):
        return self.__cxxdbg
    def fcdbg(self):
        return self.__fcdbg

system_configurations=[ 
    Configuration("generic",
        "gcc","g++","gfortran",
        "none","none",
        "-DNDEBUG -O3","-DNDEBUG -O3","-O3",
        "-O0 -g","-O0 -g","-O0 -g"),
    Configuration("psc-bridges2",
        "icc","icpc","ifort",
        "none","intel/20.4;openmpi/4.0.2-intel20.4",
        "-DNDEBUG -O3 -mavx2","-DNDEBUG -O3 -mavx2","-O3 -mavx2",
        "-O0 -g","-O0 -g","-O0 -g -traceback"),
    Configuration("und-aegaeon",
        "icc","icpc","ifort",
        "none","intel/19.0;mvapich2/2.3.1/intel/19.0",
        "-DNDEBUG -O3 -xCORE-AVX2","-DNDEBUG -O3 -xCORE-AVX2","-O3 -xCORE-AVX2",
        "-O0 -g","-O0 -g","-O0 -g -traceback")
]

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def available_systems():
    s = []
    for S in system_configurations:
        s.append(S.name())
    return s

def generate_executable_defs(build):
    defs=""
    if "adcirc" in build:
        defs+="-DBUILD_ADCIRC=ON "
    if "padcirc" in build:
        defs+="-DBUILD_PADCIRC=ON "
    if "adcprep" in build:
        defs+="-DBUILD_ADCPREP=ON "
    if "adcswan" in build:
        defs+="-DBUILD_ADCSWAN=ON "
    if "padcswan" in build:
        defs+="-DBUILD_PADCSWAN=ON "
    if "aswip" in build:
        defs+="-DBUILD_ASWIP=ON "
    return defs

def get_configuration(system_name):
    for S in system_configurations:
        if S.name() == system_name:
            return S
    print("System configuration for "+system_name+" could not be found.")
    sys.exit(1)

def main():
    import argparse
    import sys
    import os
    
    parser = argparse.ArgumentParser(description="ADCIRC CMake shortcuts for HPC systems")
    parser.add_argument('--system',metavar="<name>",type=str,help="Name of the system to configure CMake for")
    parser.add_argument('--list',help="List the names of available systems",action="store_true")
    parser.add_argument('--debug',help="Configure build in debug mode (i.e. -O0)",action="store_true")
    parser.add_argument('--build',metavar="executable",help="Name of the executable to build. Repeat this argument for multiple executables",action="append")
    parser.add_argument('--netcdf',metavar="path",help="Allows use of custom netcdf library directory")
    parser.add_argument('--show',help="Show the generated command line arguments for cmake",action="store_true")
    
    args = parser.parse_args()
    
    if args.list:
        system_list = available_systems()
        for S in system_list:
            print(S)
        sys.exit(0)
    
    if not args.build:
        print("ERROR: No executables selected for building.")
        sys.exit(1)
    executable_defs = generate_executable_defs(args.build)

    if not args.system:
        print("ERROR: No system defined")
        sys.exit(1)

    config = get_configuration(args.system)
    if args.debug:
        build_type = "Debug"
    else:
        build_type = "Release"

    cmake_options = "-DCMAKE_C_COMPILER="+config.cc()+" -DCMAKE_CXX_COMPILER="+config.cxx()+ \
        " -DCMAKE_Fortran_COMPILER="+config.fc()+" "+executable_defs+ \
        "-DCMAKE_C_FLAGS_DEBUG=\""+config.cdbg()+"\" -DCMAKE_CXX_FLAGS_DEBUG=\""+config.cxxdbg()+ \
        "\" -DCMAKE_Fortran_FLAGS_DEBUG=\""+config.fcdbg()+"\" -DCMAKE_C_FLAGS_RELEASE=\""+config.copt()+ \
        "\" -DCMAKE_CXX_FLAGS_RELEASE=\""+config.cxxopt()+"\" -DCMAKE_Fortran_FLAGS_RELEASE=\""+config.fcopt() + \
        "\" -DCMAKE_BUILD_TYPE="+build_type+" "

    if args.netcdf:
        cmake_options += "-DENABLE_OUTPUT_NETCDF=ON -DNETCDFHOME="+args.netcdf
    elif not config.netcdfhome() == "none":
        cmake_options += "-DENABLE_OUTPUT_NETCDF=ON -DNETCDFHOME="+config.netcdfhome()
    
    cmake_exe = which("cmake")
    if not cmake_exe:
        print("ERROR: No cmake executable found")
        sys.exit(1)

    if not config.modules() == "none":
        module_cmd = ""
        for M in config.modules().split(";"):
            module_cmd += "module load "+M+"; "
        cmake_exe = module_cmd + cmake_exe


    if not os.path.exists("build"):
        os.mkdir("build")

    os.chdir("build")
    os.system(cmake_exe+" --build . --target clean") 
    full_cmd = cmake_exe+" .. "+cmake_options
    if args.show:
        print(full_cmd)
    ierr = os.system(full_cmd)
    if ierr == 0:
        print("CMake ran successfully. Now run:")
        print("  "+module_cmd+"cd build; make")
    else:
        print("ERROR: CMake exited with errors")
        sys.exit(1)
        

if __name__ == "__main__":
    main()
