! metis_interface.f90 -- Fortran METIS Interface
!
! Copyright 2019 Ivan Pribec <ivan.pribec@gmail.com>
! 
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
! 
!     http://www.apache.org/licenses/LICENSE-2.0
! 
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!*****************************************************************************************
!> author: Ivan Pribec
!  date: 7/2019
!
! A Fortran interface to the METIS graph partitioning library.
!
! http://glaros.dtc.umn.edu/gkhome/node/877
!
module metis_interface

    use iso_c_binding, only: c_int32_t, c_int64_t, c_float, c_double, c_ptr

    implicit none
    private

    !
    ! Width of elementary data types (should match those used in metis.h)
    !
#ifdef INT64
    integer, parameter, public :: idx_t = c_int64_t ! <--- modify integer size here (c_int32_t or c_int64_t)
#else
    integer, parameter, public :: idx_t = c_int32_t ! <--- modify integer size here (c_int32_t or c_int64_t)
#endif

#ifdef REAL64
    integer, parameter, public :: real_t = c_double  ! <--- modify real size here (c_float or c_double)
#else
    integer, parameter, public :: real_t = c_float
#endif

    !
    ! Number of METIS options
    !
    integer(kind=idx_t), parameter, public :: METIS_NOPTIONS = 40       !! Number of METIS OPTIONS

    !
    !
    ! Enum type definitions
    !

    ! Return codes
    integer(kind=idx_t), parameter, public :: METIS_OK = 1              !! Returned normally.
    integer(kind=idx_t), parameter, public :: METIS_ERROR_INPUT = -2    !! Returned due to erroneous inputs and/or options.
    integer(kind=idx_t), parameter, public :: METIS_ERROR_MEMORY = -3   !! Returned due to insufficient memory.
    integer(kind=idx_t), parameter, public :: METIS_ERROR = -4          !! Return status: some other type of error.

    ! Operation type codes
    ! integer(kind=idx_t), parameter, public :: METIS_OP_PMETIS = 0       
    ! integer(kind=idx_t), parameter, public :: METIS_OP_KMETIS = 1
    ! integer(kind=idx_t), parameter, public :: METIS_OP_OMETIS = 2

    ! Options codes (i.e., `options`)
    integer(kind=idx_t), parameter, public :: METIS_OPTION_PTYPE     = 0  !! Specifies the partitioning method.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_OBJTYPE   = 1  !! Specifies the type of objective.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_CTYPE     = 2  !! Specifies the matching scheme to be used during coarsening.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_IPTYPE    = 3  !! Determines the algorithm used during initial partitioning.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_RTYPE     = 4  !! Determines the algorithm used for refinement.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_DBGLVL    = 5  !! Specifies the amount of progress/debugging information will be printed.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_NITER     = 6  !! Specifies the number of iterations for the refinement algorithm.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_NCUTS     = 7  !! Specifies the number of different partitionings that it will compute.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_SEED      = 8  !! Specifies the seed for the random number generator.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_NO2HOP    = 9  !! Specifies that the coarsening will not perform any 2–hop matchings when the standard matching approach fails to sufficiently coarsen the graph.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_MINCONN   = 10 !! Specifies that the partitioning routines should try to minimize the maximum degree of the subdomain graph.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_CONTIG    = 11 !! Specifies that the partitioning routines should try to produce partitions that are contigous.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_COMPRESS  = 12 !! Specifies that the graph should be compressed by combining together vertices that have identical adjacency lists.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_CCORDER   = 13 !! Specifies if the connected components of the graph should first be identifies and ordered separately.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_PFACTOR   = 14 !! Specifies the minimum degree of the vertices that will be ordered last.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_NSEPS     = 15 !! Specifies the number of different separators that it will compute at each level of nested dissection.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_UFACTOR   = 16 !! Specifies the maximum allowed load imbalance among the partitions.
    integer(kind=idx_t), parameter, public :: METIS_OPTION_NUMBERING = 17 !! Used to indicate which numbering scheme is used for the adjacency structure of a graph or the element-node structure of a mesh

    ! Partitioning Schemes
    integer(kind=idx_t), parameter, public :: METIS_PTYPE_RB   = 0 !! Multilevel recursive bisectioning.
    integer(kind=idx_t), parameter, public :: METIS_PTYPE_KWAY = 1 !! Multilevel \(k\)-way partitioning.

    ! Graph types for meshes
    ! integer(kind=idx_t), parameter, public :: METIS_GTYPE_DUAL  = 0
    ! integer(kind=idx_t), parameter, public :: METIS_GTYPE_NODAL = 1               

    ! Coarsening Schemes
    integer(kind=idx_t), parameter, public :: METIS_CTYPE_RM   = 0 !! Random matching.
    integer(kind=idx_t), parameter, public :: METIS_CTYPE_SHEM = 1 !! Sorted heavy-edge matching.

    ! Initial partitioning schemes
    integer(kind=idx_t), parameter, public :: METIS_IPTYPE_GROW    = 0 !! Grows a bisection using a greedy strategy.
    integer(kind=idx_t), parameter, public :: METIS_IPTYPE_RANDOM  = 1 !! Computes a bisection at random followed by a refinement.
    integer(kind=idx_t), parameter, public :: METIS_IPTYPE_EDGE    = 2 !! Derives a separator form an edge cut.
    integer(kind=idx_t), parameter, public :: METIS_IPTYPE_NODE    = 3 !! Grows a bisection using a greedy node-based strategy.
    integer(kind=idx_t), parameter, public :: METIS_IPTYPE_METISRB = 4

    ! Refinement schemes
    integer(kind=idx_t), parameter, public :: METIS_RTYPE_FM        = 0 !! FM-based cut refinement.
    integer(kind=idx_t), parameter, public :: METIS_RTYPE_GREEDY    = 1 !! Greedy-based cut and volume refinement.
    integer(kind=idx_t), parameter, public :: METIS_RTYPE_SEP2SIDED = 2 !! Two-sided node FM refinement.
    integer(kind=idx_t), parameter, public :: METIS_RTYPE_SEP1SIDED = 3 !! One-sided node FM refinement.

    ! Debug Levels
    integer(kind=idx_t), parameter, public :: METIS_DBG_INFO       = 1       !! Shows various diagnostic messages.
    integer(kind=idx_t), parameter, public :: METIS_DBG_TIME       = 2       !! Perform timing analysis.
    integer(kind=idx_t), parameter, public :: METIS_DBG_COARSEN    = 4       !! Show the coarsening progress.
    integer(kind=idx_t), parameter, public :: METIS_DBG_REFINE     = 8       !! Show the refinement progress.
    integer(kind=idx_t), parameter, public :: METIS_DBG_IPART      = 16      !! Show info on initial partitioning.
    integer(kind=idx_t), parameter, public :: METIS_DBG_MOVEINFO   = 32      !! Show info on vertex moves during refinement.
    integer(kind=idx_t), parameter, public :: METIS_DBG_SEPINFO    = 64      !! Show info on vertex moves during sep refinement.
    integer(kind=idx_t), parameter, public :: METIS_DBG_CONNINFO   = 128     !! Show info on minimization of subdomain connectivity.
    integer(kind=idx_t), parameter, public :: METIS_DBG_CONTIGINFO = 256     !! Show info on elimination of connected components. 
    integer(kind=idx_t), parameter, public :: METIS_DBG_MEMORY     = 2048    !! Show info related to wspace allocation.
    
    ! Types of objectives
    integer(kind=idx_t), parameter, public :: METIS_OBJTYPE_CUT  = 0 !! Edge-cut minimization.
    integer(kind=idx_t), parameter, public :: METIS_OBJTYPE_VOL  = 1 !! Total communication volume minimization.
    integer(kind=idx_t), parameter, public :: METIS_OBJTYPE_NODE = 2

    !
    ! METIS API
    !

    ! Graph partitioning routines
    public :: METIS_PartGraphRecursive
    public :: METIS_PartGraphKway

    ! Mesh partitioning routines
    public :: METIS_PartMeshDual
    public :: METIS_PartMeshNodal

    ! Sparse matrix reordering routines
    public :: METIS_NodeND

    ! Mesh-to-graph conversion routines
    public :: METIS_MeshToDual
    public :: METIS_MeshToNodal

    ! Utility routines
    public :: METIS_SetDefaultOptions
    public :: METIS_Free


    interface

    !
    ! Graph partitioning routines
    !

!*****************************************************************************************
!> This function is used to partition a graph into `nparts` parts using recursive bisection.
! 
!  If `tpwgt` is present, the *target partition weight* for the `i`-th partition and `j`-th constraint should
!  be specified at `tpwgts(i*ncon+j)` (the numbering for both partitions and constraints starts from 0).
!  For each constraint, the sum of the `tpwgts`entries must be 1.0.
!
!  The following options are valid: <br />
!  `METIS_OPTION_CTYPE`, `METIS_OPTION_IPTYPE`, `METIS_OPTION_RTYPE`, 
!  `METIS_OPTION_NO2HOP`, `METIS_OPTION_NCUTS`, `METIS_OPTION_NITER`, 
!  `METIS_OPTION_SEED`, `METIS_OPTION_UFACTOR`, `METIS_OPTION_NUMBERING`,
!  `METIS_OPTION_DBGLVL`
!
    function METIS_PartGraphRecursive(nvtxs,ncon,xadj,adjncy,&
        vwgt,vsize,adjwgt,nparts,tpwgts,ubvec,options,objval,part) result(ierr) bind(C,name="METIS_PartGraphRecursive")
        import idx_t, real_t, METIS_NOPTIONS
        
        ! Parameters
        integer(kind=idx_t), intent(in) :: nvtxs
            !! The number of vertices in the graph.
        integer(kind=idx_t), intent(in) :: ncon
            !! The number of balancing constraints. It should be atleast 1.
        integer(kind=idx_t), intent(in), dimension(*) :: xadj, adjncy
            !! The adjacency structure of the graph as described in section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(kind=idx_t), intent(in), dimension(*), optional :: vwgt ! NULL
            !! The weights of the vertices as described in Section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(kind=idx_t), intent(in), dimension(*), optional :: vsize ! NULL
            !! The size of the vertices for computing the total communication volume as described in section 5.7 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(kind=idx_t), intent(in), dimension(*), optional :: adjwgt ! NULL
            !! The weights of the edges as describe in Section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(kind=idx_t) , intent(in) :: nparts
            !! The number of parts to partition the graph.
        real(kind=real_t), intent(in), optional :: tpwgts(nparts*ncon)
            !! An array of size `nparts*ncon` that specifies the desired weight for each partition and constraint.
            !! If not present, the graph is divided equally among the partitions. More in the description.
        real(kind=real_t), intent(in), optional :: ubvec(ncon)
            !! An array of size `ncon` that specifies the allowed load imbalance for each constraint. 
            !! For the `i`-th partition and `j`-th constraint the allowed weight is the `ubvec(j)*tpwgts(i*ncon+j)`
            !! fraction of the `j`-th's constraint total weight. If not present, the load imbalance
            !! tolerance is 1.001 (for `ncon = 1`) or 1.01 (for `ncon > 1`).
        integer(kind=idx_t), intent(in), optional :: options(METIS_NOPTIONS)
            !! An array of options as described in Section 5.4 of the METIS manual. See description for valid options.
        integer(kind=idx_t), intent(out) :: objval
            !! Upon successful completion, this variable stores the edge-cut or the total communication volume of the partitioning
            !! solution. The value returned depends on the partitioning's objective function.
        integer(kind=idx_t), intent(out) :: part(nvtxs)
            !! This is a vector of size `nvtxs` that upon successful completion stores the partition vector of the graph.
            !! The numbering of this vector starts from either 0 or 1, depending on the value of `options(METIS_OPTION_NUMBERING)`.

        ! Returns
        integer(kind=idx_t) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.<br /> 
            !! `METIS_ERROR_INPUT` - Indicates an input error.<br /> 
            !! `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.<br /> 
            !! `METIS_ERROR` - Indicates some other type of error.
    end function
!*****************************************************************************************

!*****************************************************************************************
!> This function is used to partition a graph into `nparts` parts using multilevel \(k\)-way partitioning.
!
!  If `tpwgt` is present, the *target partition weight* for the `i`-th partition and `j`-th constraint should
!  be specified at `tpwgts(i*ncon+j)` (the numbering for both partitions and constraints starts from 0).
!  For each constraint, the sum of the `tpwgts`entries must be 1.0.
!
!  The following options are valid: <br />
!  `METIS_OPTION_OBJTYPE`, `METIS_OPTION_CTYPE`, `METIS_OPTION_IPTYPE`,
!  `METIS_OPTION_RTYPE`, `METIS_OPTION_NO2HOP`, `METIS_OPTION_NCUTS`,
!  `METIS_OPTION_NITER`, `METIS_OPTION_UFACTOR`, `METIS_OPTION_MINCONN`,
!  `METIS_OPTION_CONTIG`, `METIS_OPTION_SEED`, `METIS_OPTION_NUMBERING`,
!  `METIS_OPTION_DBGLVL`
!
    function METIS_PartGraphKway(nvtxs,ncon,xadj,adjncy,&
        vwgt,vsize,adjwgt,nparts,tpwgts,ubvec,options,objval,part) result(ierr) bind(C,name="METIS_PartGraphKway")
        import idx_t, real_t, METIS_NOPTIONS
        
        ! Parameters
        integer(kind=idx_t), intent(in) :: nvtxs
            !! The number of vertices in the graph.
        integer(kind=idx_t), intent(in) :: ncon
            !! The number of balancing constraints. It should be at least 1.
        integer(kind=idx_t), intent(in), dimension(*) :: xadj, adjncy
            !! The adjacency structure of the graph as described in section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(kind=idx_t), intent(in), dimension(*), optional :: vwgt
            !! The weights of the vertices as described in Section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(kind=idx_t), intent(in), dimension(*), optional :: vsize
            !! The size of the vertices for computing the total communication volume as described in section 5.7 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(kind=idx_t), intent(in), dimension(*), optional :: adjwgt
            !! The weights of the edges as describe in Section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(kind=idx_t), intent(in) :: nparts
            !! The number of parts to partition the graph.
        real(kind=real_t), intent(in), optional :: tpwgts(nparts*ncon)
            !! An array of size `nparts*ncon` that specifies the desired weight for each partition and constraint.
            !! If not present, the graph is divided equally among the partitions. More in the description.
        real(kind=real_t), intent(in), optional :: ubvec(ncon)
            !! An array of size `ncon` that specifiew the allowed load imbalance for each constraint. 
            !! For the `i`-th partition and `j`-th constraint the allowed weight is the `ubvec(j)*tpwgts(i*ncon+j)`
            !! fraction of the `j`-th's constraint total weight. If not present, the load imbalance
            !! tolerance is 1.001 (for `ncon == 1`) or 1.01 (for `ncon > 1`).
        integer(kind=idx_t), intent(in), optional :: options(METIS_NOPTIONS)
            !! An array of options as described in Section 5.4 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf). See description for valid options.
        integer(kind=idx_t), intent(out) :: objval
            !! Upon successful completion, this variable stores the edge-cut or the total communication volume of the partitioning
            !! solution. The value returned depends on the partitioning's objective function.
        integer(kind=idx_t), intent(out) :: part(nvtxs)
            !! This is a vector of size `nvtxs` that upon successful completion stores the partition vector of the graph.
            !! The numbering of this vector starts from either 0 or 1, depending on the value of `options(METIS_OPTION_NUMBERING)`.
        
        ! Returns
        integer(kind=idx_t) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.<br /> 
            !! `METIS_ERROR_INPUT` - Indicates an input error.<br /> 
            !! `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.<br /> 
            !! `METIS_ERROR` - Indicates some other type of error.
    end function
!*****************************************************************************************

    !
    ! Mesh partitioning routines
    !

!*****************************************************************************************
!> This function is used to partition a mesh into `nparts` parts based on a partitioning of the mesh's dual graph.
!
!  The following options are valid: <br />
!  `METIS_OPTION_PTYPE`, `METIS_OPTION_OBJTYPE`, `METIS_OPTION_CTYPE`,
!  `METIS_OPTION_IPTYPE`, `METIS_OPTION_RTYPE`, `METIS_OPTION_NCUTS`,
!  `METIS_OPTION_NITER`, `METIS_OPTION_SEED`, `METIS_OPTION_UFACTOR`,
!  `METIS_OPTION_NUMBERING`, `METIS_OPTION_DBGLVL`
!
    function METIS_PartMeshDual(ne,nn,eptr,eind,vwgt,vsize,ncommon, &
        nparts,tpwgts,options,objval,epart,npart) result(ierr) bind(C,name="METIS_PartMeshDual")
        import idx_t, METIS_NOPTIONS
        integer(kind=idx_t), intent(in) :: ne
            !! The number of elements in the mesh.
        integer(kind=idx_t), intent(in) :: nn
            !! The number of nodes in the mesh.
        integer(kind=idx_t), intent(in), dimension(*) :: eptr,eind
            !! The pair of arrays storing the mesh as described in Section 5.6 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(kind=idx_t), intent(in), optional :: vwgt(ne)
            !! An array of size `ne` specifying the weights of the elements. If not present,
            !! all elements have an equal weight.
        integer(kind=idx_t), intent(in), optional :: vsize(ne)
            !! An array of size `ne` specifying the size of the elements that is used
            !! for computing the total comunication volume as described in Section 5.7 of the  [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
            !! If not present, the objective is cut or all elements have an equal size.
        integer(kind=idx_t), intent(in) :: ncommon
        integer(kind=idx_t), intent(in) :: nparts
            !! The number of parts to partition the mesh.
        integer(kind=idx_t), intent(in), optional :: tpwgts(nparts)
            !! An array of size `nparts` that specifies the desired weight for each partition. The *target
            !! partition weight* for the `i`-th partition is specified at `tpwgts(i)` (the numbering for the 
            !! partitions starts from 0). The sum of the `tpwgts` entries must be 1.0. <br /> If not present, the graph
            !! is divided equally among the partitions.
        integer(kind=idx_t), intent(in), optional :: options(METIS_NOPTIONS)
            !! An array of options as described in Section 5.4 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf). See description for valid options.
        integer(kind=idx_t), intent(out) :: objval
            !! Upon successful completion, this variable stores either the edgecut or the total communication
            !! volume of the dual graph's partitioning.
        integer(kind=idx_t), intent(out) :: epart(ne)
            !! A vector of size `ne` that upon successful completion stores the partition vector for the elements
            !! of the mesh. The numbering of this vector starts from either 0 or 1, depending on the value of
            !! `options(METIS_OPTION_NUMBERING)`.
        integer(kind=idx_t), intent(out) :: npart(nn)
            !! A vector of size `nn` that upon successful completion stores the partition vector for the nodes 
            !! of the mesh. The numbering of this vector starts from either 0 or 1, depending on the value of
            !! `options(METIS_OPTION_NUMBERING)`.

        ! Returns
        integer(kind=idx_t) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.<br /> 
            !! `METIS_ERROR_INPUT` - Indicates an input error.<br /> 
            !! `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.<br /> 
            !! `METIS_ERROR` - Indicates some other type of error.
    end function
!*****************************************************************************************

!*****************************************************************************************
!> This function us used to partition a mesh into `nparts` parts based on a 
!  partitioning of the mesh's nodal graph.
!
!  The following options are valid: <br />
!  `METIS_OPTION_PTYPE`, `METIS_OPTION_OBJTYPE`, `METIS_OPTION_CTYPE`,
!  `METIS_OPTION_IPTYPE`, `METIS_OPTION_RTYPE`, `METIS_OPTION_NCUTS`,
!  `METIS_OPTION_NITER`, `METIS_OPTION_SEED`, `METIS_OPTION_UFACTOR`,
!  `METIS_OPTION_NUMBERING`, `METIS_OPTION_DBGLVL`
!
    function METIS_PartMeshNodal(ne,nn,eptr,eind,vwgt,vsize, &
        nparts,tpwgts,options,objval,epart,npart) result(ierr) bind(C,name="METIS_PartMeshNodal")
        import idx_t, METIS_NOPTIONS
        integer(kind=idx_t), intent(in) :: ne
            !! The number of elements in the mesh.
        integer(kind=idx_t), intent(in) :: nn
            !! The number of nodes in the mesh.
        integer(kind=idx_t), intent(in), dimension(*) :: eptr,eind
            !! The pair of arrays storing the mesh as described in Section 5.6 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(kind=idx_t), intent(in), optional :: vwgt(nn)
            !! An array of size `nn` specifying weights of the nodes. If not passed, all nodes have an equal weight.
        integer(kind=idx_t), intent(in), optional :: vsize(nn)
            !! An array of size `nn` specifying the size of the nodes that is used for computing the
            !! total comunication volume as described in Section 5.7 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf). If not passed,
            !! the objective is cut or all nodes have an equal size.
        integer(kind=idx_t), intent(in) :: nparts
            !! The number of parts to partition the mesh.
        integer(kind=idx_t), intent(in), optional :: tpwgts(nparts)
            !! An array of size `nparts` that specifies the desired weight for each partition. The *target
            !! partition weight* for the `i`-th partition is specified at `tpwgts(i)` (the numbering for the 
            !! partitions starts from 0). The sum of the `tpwgts` entries must be 1.0. If not passed, the graph
            !! is divided equally among the partitions.
        integer(kind=idx_t), intent(in), optional :: options(METIS_NOPTIONS)
            !! An array of options as described in Section 5.4 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf). See description for valid options.
        integer(kind=idx_t), intent(out) :: objval
            !! Upon successful completion, this variable stores either the edgecut or the total communication
            !! volume of the nodal graph's partitioning.
        integer(kind=idx_t), intent(out) :: epart(ne)
            !! A vector of size `ne` that upon successful completion stores the partition vector for the elements
            !! of the mesh. The numbering of this vector starts from either 0 or 1, depending on the value of
            !! `options(METIS_OPTION_NUMBERING)`.
        integer(kind=idx_t), intent(out) :: npart(nn)
            !! A vector of size `nn` that upon successful completion stores the partition vector for the nodes 
            !! of the mesh. The numbering of this vector starts from either 0 or 1, depending on the value of
            !! `options(METIS_OPTION_NUMBERING)`.

        ! Returns
        integer(kind=idx_t) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.<br /> 
            !! `METIS_ERROR_INPUT` - Indicates an input error.<br /> 
            !! `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.<br /> 
            !! `METIS_ERROR` - Indicates some other type of error.
    end function
!*****************************************************************************************

    !
    ! Sparse Matrix Reordering Routines
    !

!*****************************************************************************************
!> This function computes fill reducing orderings of sparse matrices using the
!  multilevel nested dissection algorithm.
!
!  Let \(A\) be the original matrix and \(A'\) be the permuted matrix. The arrays
!  `perm` and `iperm` are defined as follows. Row (column) `i` of \(A'\) is the
!  `perm(i)` row (column) of \(A\), and row (column) `i` of \(A\) is the `iperm(i)`
!  row (column) of \(A'\). the numbering of this vector starts from either 0 or 1,
!  depending on the value of `options(METIS_OPTION_NUMBERING)`.
!
!  If the graph is weighted, meaning `vgwt` was provided, the nested dissection ordering computes
!  vertex separators that minimize the sum of the weights of the vertices on the separators.
!
!  The following options are valid: <br />
!  `METIS_OPTION_CTYPE`, `METIS_OPTION_RTYPE`, `METIS_OPTION_NO2HOP`,
!  `METIS_OPTION_NSEPS`, `METIS_OPTION_NITER`, `METIS_OPTION_UFACTOR`,
!  `METIS_OPTION_COMPRESS`, `METIS_OPTION_CCORDER`, `METIS_OPTION_SEED`,
!  `METIS_OPTION_PFACTOR`, `METIS_OPTION_NUMBERING`, `METIS_OPTION_DBGLVL`
!
!#Example
! The code below generates the nodal graph of the following mesh:
!```Fortran
! use metis_interface, only: idx_t, METIS_SetDefaultOptions, METIS_NodeND
!integer(kind=idx_t), parameter :: n = 15, m = 22
!integer(kind=idx_t) :: xadj(n+1), adjncy(2*m)
!integer(kind=idx_t) :: perm(n), iperm(n)
!integer(kind=idx_t) :: options(0:39), ierr
!
!xadj = [1,3,6,9,12,14,17,21,25,29,32,34,37,40,43,45]
!adjncy = [2,6,1,3,7,2,4,8,3,5,9,4,10,1,7,11,2,6, &
!                  8,12,3,7,9,13,4,8,10,14,5,9,15,6,12,7,11,13, &
!                  8,12,14,9,13,15,10,14]
!
!ierr = METIS_SetDefaultOptions(options)
!options(18) = 1
!
!ierr = METIS_NodeND(n,xadj,adjncy,options=options,perm=perm,iperm=iperm)
!end
!```
    function METIS_NodeND(nvtxs,xadj,adjncy,vwgt,options,perm,iperm) result(ierr) bind(C,name="METIS_NodeND")
        import idx_t, METIS_NOPTIONS

        ! Parameters
        integer(kind=idx_t), intent(in) :: nvtxs
            !! The number of vertices in the graph.
        integer(kind=idx_t), intent(in), dimension(*) :: xadj, adjncy
            !! The adjacency structure of the graph as described in Section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(kind=idx_t), intent(in), optional :: vwgt(nvtxs)
            !! An array of size `nvtxs` specifying the weights of the vertices.
        integer(kind=idx_t), intent(in), optional :: options(METIS_NOPTIONS)
            !! This is the array of options as described in Section 5.4 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf). See description for valid options.
        integer(kind=idx_t), intent(out) :: perm(nvtxs), iperm(nvtxs)
            !! Vectors of size `nvtxs`. Upon successful completion, they store the fill-reducing
            !! permutation and inverse-permutation. More in the description.

        ! Returns
        integer(kind=idx_t) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.<br /> 
            !! `METIS_ERROR_INPUT` - Indicates an input error.<br /> 
            !! `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.<br /> 
            !! `METIS_ERROR` - Indicates some other type of error.
    end function
!*****************************************************************************************

    !
    ! Mesh-to-graph conversion routines
    !

!*****************************************************************************************
!> This function is used to generate the dual graph of a mesh.
!
!@note
! To use the returned arrays `xadj` and `adjncy`, these must be first converted from
! a C pointer to a Fortran pointer using the subroutine `c_f_pointer(cptr,fptr,shape)` 
! that assigns the target of the C pointer `cptr` to the Fortran pointer `fptr` and
! specifies its shape. The `shape` is an integer rank-one array, storing the size `ne+1`
! in case of the dual graph. The size of the new `adjncy` array is stored in the 
! last element of `xadj` when using C-style numbering. An example is shown below.
!@endnote
!
!@warning
! Memory for the returned arrays `xadj` and `adjncy` is allocated by METIS' API in C
! using the standard `malloc` function. It is the responsibility of the application to free
! this memory by calling `free`. Therefore, METIS provides the [[METIS_Free]] function that is a wrapper to
! C's `free`function.
!@endwarning
!
!# Example
! The code below generates the nodal graph of the following mesh:
! __image__
!```Fortran
!use iso_c_binding, only : c_ptr, c_f_pointer
!use metis_interface, only: idx_t, METIS_MeshToDual, METIS_Free
!integer(kind=idx_t), parameter :: ne = 3, nn = 8, npel = 4
!integer(kind=idx_t) :: ierr, eptr(ne+1), eind(ne*npel), numflag, ncommon
!type(c_ptr) :: xadj, adjncy
!integer(kind=idx_t), dimension(:), pointer :: fxadj => null(), fadjncy => null()
!
!numflag = 0 ! C-style numbering
!ncommon = 2 ! 2 common nodes per edge
!eptr = [0,4,8,12]
!eind = [0,1,2,3,1,4,5,2,4,6,7,5] ! note four nodes per element
!
!ierr = METIS_MeshToDual(ne,nn,eptr,eind,ncommon,numflag,xadj,adjncy)
!call c_f_pointer(xadj,fxadj,shape=[ne+1]) ! xadj is of the size ne+1
!call c_f_pointer(adjncy,fadjncy,shape=[fxadj(ne+1)]) ! last value in xadj is the size of adjncy
!
! !... use values in fxadj and fadjncy ...
!
!call METIS_Free(xadj)
!call METIS_Free(adjncy)
!end
!```
    function METIS_MeshToDual(ne,nn,eptr,eind,ncommon,numflag,xadj,adjncy) result(ierr) bind(C,name="METIS_MeshToDual")
        import idx_t, c_ptr
        
        ! Parameters
        integer(kind=idx_t), intent(in) :: ne
            !! The number of elements in the mesh.
        integer(kind=idx_t), intent(in) :: nn
            !! The number of nodes in the mesh.
        integer(kind=idx_t), intent(in), dimension(*) :: eptr, eind
            !! The pair of arrays storing the mesh as described in Section 5.6 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(kind=idx_t), intent(in) :: ncommon
            !! The number of common nodes that two elements must have in order to put
            !! an edge between them in the dual graph.
        integer(kind=idx_t), intent(in) :: numflag
            !! Used to indicate which numbering scheme is used for `eptr` and `eind`. 
            !! The possible values are: <br />
            !! 0 - C-style numbering is assumed that starts from 0 <br />
            !! 1 - Fortran-style numbering is assumed that starts from 1
        type(c_ptr), intent(out) :: xadj, adjncy
            !! These arrays store the adjacency structure of the generated dual graph. 
            !! The format of the adjacency structure is described in Section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).

        ! Returns
        integer(kind=idx_t) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.<br /> 
            !! `METIS_ERROR_INPUT` - Indicates an input error.<br /> 
            !! `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.<br /> 
            !! `METIS_ERROR` - Indicates some other type of error.
    end function
!*****************************************************************************************

!*****************************************************************************************
!> This function is used to generate the nodal graph of a mesh.
!
!@note
! To use the returned arrays `xadj` and `adjncy`, these must be first converted from
! a C pointer to a Fortran pointer using the subroutine `c_f_pointer(cptr,fptr,shape)` 
! that assigns the target of the C pointer `cptr` to the Fortran pointer `fptr` and
! specifies its shape. The `shape` is an integer rank-one array, storing the size `nn+1`
! in case of the nodal graph. The size of the new `adjncy` array is stored in the 
! last element of `xadj` when using C-style numbering. An example is shown below.
!@endnote
!
!@warning
! Memory for the returned arrays `xadj` and `adjncy` is allocated by METIS' API in C
! using the standard `malloc` function. It is the responsibility of the application to free
! this memory by calling `free`. Therefore, METIS provides the [[METIS_Free]] function that is a wrapper to
! C's `free`function.
!@endwarning
!
!# Example
! The code below generates the nodal graph of the following mesh:
! __image__
!```Fortran
!use iso_c_binding, only : c_ptr, c_f_pointer
!use metis_interface, only: idx_t, METIS_MeshToNodal, METIS_Free
!integer(kind=idx_t), parameter :: ne = 3, nn = 8, npel = 4
!integer(kind=idx_t) :: ierr, eptr(ne+1), eind(ne*npel), numflag
!type(c_ptr) :: xadj, adjncy
!integer(kind=idx_t), dimension(:), pointer :: fxadj => null(), fadjncy => null()
!
!numflag = 0 ! C-style numbering
!eptr = [0,4,8,12]
!eind = [0,1,2,3,1,4,5,2,4,6,7,5] ! note four nodes per element
!
!ierr = METIS_MeshToNodal(ne,nn,eptr,eind,numflag,xadj,adjncy)
!call c_f_pointer(xadj,fxadj,shape=[nn+1]) ! xadj is of the size nn+1
!call c_f_pointer(adjncy,fadjncy,shape=[fxadj(nn+1)]) ! last value in xadj is the size of adjncy
!
! !... use values in fxadj and fadjncy ...
!
!call METIS_Free(xadj)
!call METIS_Free(adjncy)
!end
!```
    function METIS_MeshToNodal(ne,nn,eptr,eind,numflag,xadj,adjncy) result(ierr) bind(C,name="METIS_MeshToNodal")
        import idx_t, c_ptr
        
        ! Parameters
        integer(kind=idx_t), intent(in) :: ne
            !! The number of elements in the mesh.
        integer(kind=idx_t), intent(in) :: nn
            !! The number of nodes in the mesh.
        integer(kind=idx_t), intent(in), dimension(*) :: eptr, eind
            !! The pair of arrays storing the mesh as described in Section 5.6 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).
        integer(kind=idx_t), intent(in) :: numflag
            !! Used to indicate which numbering scheme is used for `eptr` and `eind`.
            !! The possible values are: <br />
            !! 0 - C-style numbering is assumed that starts from 0 <br />
            !! 1 - Fortran-style numbering is assumed that starts from 1
        type(c_ptr), intent(out) :: xadj, adjncy
            !! These arrays store the adjacency structure of the generated dual graph. 
            !! The format of the adjacency structure is described in Section 5.5 of the [manual](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf).

        ! Returns
        integer(kind=idx_t) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.<br /> 
            !! `METIS_ERROR_INPUT` - Indicates an input error.<br /> 
            !! `METIS_ERROR_MEMORY` - Indicates that it could not allocate the required memory.<br /> 
            !! `METIS_ERROR` - Indicates some other type of error.
    end function
!*****************************************************************************************

    !
    ! Utility routines
    !

!*****************************************************************************************
!> Initializes the options array into its default values.
!
!@note
! The passed array `options` must have the size `METIS_NOPTIONS` (40).
! To be able to use the option parameters specified in the [[metis_interface]] module 
! it is recommended to use zero-based indexing for the options array: 
!```Fortran
!integer(kind=idx_t) :: opts(0:39)
!```
!@endnote
!
!# Examples
! To set Fortran style index-numbering use:
!```Fortran
!use metis_enum, only: METIS_OPTION_NUMBERING
!use metis_interface, only: idx_t, METIS_SetDefaultOptions
!integer(kind=idx_t) :: opts(0:39)
! 
!call METIS_SetDefaultOptions(opts)
!opts(METIS_OPTION_NUMBERING) = 1 ! Fortran-style index numbering
!```
!
! Other options can also be changed using parameters specified in the [[metis_interface]] module.
!```Fortran
!use metis_enum, only : METIS_OPTION_DBGLVL, METIS_DBG_INFO
!use metis_interface, only : idx_t, METIS_SetDefaultOptions
!integer(kind=idx_t) :: opts(0:39)
! 
!call METIS_SetDefaultOptions(opts)
!opts(METIS_OPTION_DBGLVL) = METIS_DBG_INFO ! Show various diagnostic messages
!end
!```
    function METIS_SetDefaultOptions(options) result(ierr) bind(C,name="METIS_SetDefaultOptions")
        import idx_t, METIS_NOPTIONS
        
        ! Parameters
        integer(kind=idx_t), intent(out) :: options(METIS_NOPTIONS)
            !! The array of options that will be initialized.

        ! Returns
        integer(kind=idx_t) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.
    end function
!*****************************************************************************************

!*****************************************************************************************
!> Frees the memory that was allocated by either the [[METIS_MeshToDual]] or the
!  [[METIS_MeshToNodal]] routines for returning the dual or nodal graph of a mesh.
!
!@warning Memory deallocation should always happen on the same side it was allocated!
! Also check the descriptions of the above-mentioned routines.
!
!# Example
!
!```Fortran
! type(c_ptr) :: xadj(:),adjncy(:)
! 
! call METIS_MeshToNodal(...,xadj,adjncy)
! 
! ! xadj and adjncy should be deallocated on the C side! ;)
! call METIS_Free(xadj)
! call METIS_Free(adjncy)
!```   
    function METIS_Free(ptr) result(ierr) bind(C,name="METIS_Free")
        import idx_t, c_ptr
        
        ! Parameters
        type(c_ptr), value :: ptr
            !! The pointer to be freed. This pointer should be one of the `xadj` or `adjncy`
            !! arrays returned by METIS' API routines.

        ! Returns
        integer(kind=idx_t) :: ierr
            !! `METIS_OK` - Indicates that the function returned normally.
    end function
!*****************************************************************************************

    end interface

end module
!*****************************************************************************************





