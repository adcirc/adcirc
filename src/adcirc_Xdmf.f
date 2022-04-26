      integer xdmfaddattribute
      integer xdmfaddinformation
      integer xdmfsetgeometry
      integer xdmfsettopology
      integer xdmfsetdimensions
      integer xdmfsetorigin
      integer xdmfsetbrick
      integer xdmfstoremap
      integer xdmfaddcoordinate
      integer xdmfaddset

      ! Array Type
      integer XDMF_ARRAY_TYPE_INT8
      integer XDMF_ARRAY_TYPE_INT16
      integer XDMF_ARRAY_TYPE_INT32
      integer XDMF_ARRAY_TYPE_INT64
      integer XDMF_ARRAY_TYPE_UINT8
      integer XDMF_ARRAY_TYPE_UINT16
      integer XDMF_ARRAY_TYPE_UINT32
      integer XDMF_ARRAY_TYPE_FLOAT32
      integer XDMF_ARRAY_TYPE_FLOAT64

      ! Attribute Center
      integer XDMF_ATTRIBUTE_CENTER_GRID
      integer XDMF_ATTRIBUTE_CENTER_CELL
      integer XDMF_ATTRIBUTE_CENTER_FACE
      integer XDMF_ATTRIBUTE_CENTER_EDGE
      integer XDMF_ATTRIBUTE_CENTER_NODE

      ! Attribute Type
      integer XDMF_ATTRIBUTE_TYPE_SCALAR
      integer XDMF_ATTRIBUTE_TYPE_VECTOR
      integer XDMF_ATTRIBUTE_TYPE_TENSOR
      integer XDMF_ATTRIBUTE_TYPE_MATRIX
      integer XDMF_ATTRIBUTE_TYPE_TENSOR6
      integer XDMF_ATTRIBUTE_TYPE_GLOBALID
      integer XDMF_ATTRIBUTE_TYPE_NOTYPE

      ! Geometry Type
      integer XDMF_GEOMETRY_TYPE_XYZ
      integer XDMF_GEOMETRY_TYPE_XY

      ! Grid Collection Type
      integer XDMF_GRID_COLLECTION_TYPE_SPATIAL
      integer XDMF_GRID_COLLECTION_TYPE_TEMPORAL

      ! Topology Type
      integer XDMF_TOPOLOGY_TYPE_POLYVERTEX
      integer XDMF_TOPOLOGY_TYPE_POLYLINE
      integer XDMF_TOPOLOGY_TYPE_POLYGON
      integer XDMF_TOPOLOGY_TYPE_TRIANGLE
      integer XDMF_TOPOLOGY_TYPE_QUADRILATERAL
      integer XDMF_TOPOLOGY_TYPE_TETRAHEDRON
      integer XDMF_TOPOLOGY_TYPE_PYRAMID
      integer XDMF_TOPOLOGY_TYPE_WEDGE
      integer XDMF_TOPOLOGY_TYPE_HEXAHEDRON
      integer XDMF_TOPOLOGY_TYPE_EDGE_3
      integer XDMF_TOPOLOGY_TYPE_TRIANGLE_6
      integer XDMF_TOPOLOGY_TYPE_QUADRILATERAL_8
      integer XDMF_TOPOLOGY_TYPE_QUADRILATERAL_9
      integer XDMF_TOPOLOGY_TYPE_TETRAHEDRON_10
      integer XDMF_TOPOLOGY_TYPE_PYRAMID_13
      integer XDMF_TOPOLOGY_TYPE_WEDGE_15
      integer XDMF_TOPOLOGY_TYPE_WEDGE_18
      integer XDMF_TOPOLOGY_TYPE_HEXAHEDRON_20
      integer XDMF_TOPOLOGY_TYPE_HEXAHEDRON_24
      integer XDMF_TOPOLOGY_TYPE_HEXAHEDRON_27
      integer XDMF_TOPOLOGY_TYPE_HEXAHEDRON_64
      integer XDMF_TOPOLOGY_TYPE_HEXAHEDRON_125
      integer XDMF_TOPOLOGY_TYPE_HEXAHEDRON_216
      integer XDMF_TOPOLOGY_TYPE_HEXAHEDRON_343
      integer XDMF_TOPOLOGY_TYPE_HEXAHEDRON_512
      integer XDMF_TOPOLOGY_TYPE_HEXAHEDRON_729
      integer XDMF_TOPOLOGY_TYPE_HEXAHEDRON_1000
      integer XDMF_TOPOLOGY_TYPE_HEXAHEDRON_1331
      integer XDMF_TOPOLOGY_TYPE_MIXED

      ! Set Type
      integer XDMF_SET_TYPE_NODE
      integer XDMF_SET_TYPE_CELL
      integer XDMF_SET_TYPE_FACE
      integer XDMF_SET_TYPE_EDGE

      ! Grid Type
      integer XDMF_GRID_TYPE_CURVILINEAR
      integer XDMF_GRID_TYPE_RECTILINEAR
      integer XDMF_GRID_TYPE_REGULAR
      integer XDMF_GRID_TYPE_UNSTRUCTURED

      !------------------------------------------------------

      parameter (XDMF_ARRAY_TYPE_INT8    = 0)
      parameter (XDMF_ARRAY_TYPE_INT16   = 1)
      parameter (XDMF_ARRAY_TYPE_INT32   = 2)
      parameter (XDMF_ARRAY_TYPE_INT64   = 3)
      parameter (XDMF_ARRAY_TYPE_UINT8   = 4)
      parameter (XDMF_ARRAY_TYPE_UINT16  = 5)
      parameter (XDMF_ARRAY_TYPE_UINT32  = 6)
      parameter (XDMF_ARRAY_TYPE_FLOAT32 = 7)
      parameter (XDMF_ARRAY_TYPE_FLOAT64 = 8)

      parameter (XDMF_ATTRIBUTE_CENTER_GRID = 100)
      parameter (XDMF_ATTRIBUTE_CENTER_CELL = 101)
      parameter (XDMF_ATTRIBUTE_CENTER_FACE = 102)
      parameter (XDMF_ATTRIBUTE_CENTER_EDGE = 103)
      parameter (XDMF_ATTRIBUTE_CENTER_NODE = 104)

      parameter (XDMF_ATTRIBUTE_TYPE_SCALAR   = 200)
      parameter (XDMF_ATTRIBUTE_TYPE_VECTOR   = 201)
      parameter (XDMF_ATTRIBUTE_TYPE_TENSOR   = 202)
      parameter (XDMF_ATTRIBUTE_TYPE_MATRIX   = 203)
      parameter (XDMF_ATTRIBUTE_TYPE_TENSOR6  = 204)
      parameter (XDMF_ATTRIBUTE_TYPE_GLOBALID = 205)
      parameter (XDMF_ATTRIBUTE_TYPE_NOTYPE   = 206)

      parameter (XDMF_GEOMETRY_TYPE_XYZ  = 301)
      parameter (XDMF_GEOMETRY_TYPE_XY   = 302)

      parameter (XDMF_GRID_COLLECTION_TYPE_SPATIAL  = 400)
      parameter (XDMF_GRID_COLLECTION_TYPE_TEMPORAL = 401)

      parameter (XDMF_TOPOLOGY_TYPE_POLYVERTEX       = 500)
      parameter (XDMF_TOPOLOGY_TYPE_POLYLINE         = 501)
      parameter (XDMF_TOPOLOGY_TYPE_POLYGON          = 502)
      parameter (XDMF_TOPOLOGY_TYPE_TRIANGLE         = 503)
      parameter (XDMF_TOPOLOGY_TYPE_QUADRILATERAL    = 504)
      parameter (XDMF_TOPOLOGY_TYPE_TETRAHEDRON      = 505)
      parameter (XDMF_TOPOLOGY_TYPE_PYRAMID          = 506)
      parameter (XDMF_TOPOLOGY_TYPE_WEDGE            = 507)
      parameter (XDMF_TOPOLOGY_TYPE_HEXAHEDRON       = 508)
      parameter (XDMF_TOPOLOGY_TYPE_EDGE_3           = 509)
      parameter (XDMF_TOPOLOGY_TYPE_TRIANGLE_6       = 510)
      parameter (XDMF_TOPOLOGY_TYPE_QUADRILATERAL_8  = 511)
      parameter (XDMF_TOPOLOGY_TYPE_QUADRILATERAL_9  = 512)
      parameter (XDMF_TOPOLOGY_TYPE_TETRAHEDRON_10   = 513)
      parameter (XDMF_TOPOLOGY_TYPE_PYRAMID_13       = 514)
      parameter (XDMF_TOPOLOGY_TYPE_WEDGE_15         = 515)
      parameter (XDMF_TOPOLOGY_TYPE_WEDGE_18         = 516)
      parameter (XDMF_TOPOLOGY_TYPE_HEXAHEDRON_20    = 517)
      parameter (XDMF_TOPOLOGY_TYPE_HEXAHEDRON_24    = 518)
      parameter (XDMF_TOPOLOGY_TYPE_HEXAHEDRON_27    = 519)
      parameter (XDMF_TOPOLOGY_TYPE_HEXAHEDRON_64    = 520)
      parameter (XDMF_TOPOLOGY_TYPE_HEXAHEDRON_125   = 521)
      parameter (XDMF_TOPOLOGY_TYPE_HEXAHEDRON_216   = 522)
      parameter (XDMF_TOPOLOGY_TYPE_HEXAHEDRON_343   = 523)
      parameter (XDMF_TOPOLOGY_TYPE_HEXAHEDRON_512   = 524)
      parameter (XDMF_TOPOLOGY_TYPE_HEXAHEDRON_729   = 525)
      parameter (XDMF_TOPOLOGY_TYPE_HEXAHEDRON_1000  = 526)
      parameter (XDMF_TOPOLOGY_TYPE_HEXAHEDRON_1331  = 527)
      parameter (XDMF_TOPOLOGY_TYPE_MIXED            = 528)

      parameter (XDMF_SET_TYPE_NODE = 601)
      parameter (XDMF_SET_TYPE_CELL = 602)
      parameter (XDMF_SET_TYPE_FACE = 603)
      parameter (XDMF_SET_TYPE_EDGE = 604)

      parameter (XDMF_GRID_TYPE_CURVILINEAR	= 701)
      parameter (XDMF_GRID_TYPE_RECTILINEAR	= 702)
      parameter (XDMF_GRID_TYPE_REGULAR	= 703)
      parameter (XDMF_GRID_TYPE_UNSTRUCTURED	= 704)

