============
Installation
============

Prerequisites
=============

Before installing FortranDate, ensure you have the following prerequisites:

* A modern C++ compiler that supports at least C++17 (C++20 recommended)
* A Fortran compiler compatible with your C++ compiler
* CMake version 3.20 or newer
* `Howard Hinnant's date library <https://github.com/HowardHinnant/date>`_ (included)

Supported Compilers
===================

The library has been tested with the following compiler combinations:

* GCC 9.0+ with GFortran 9.0+
* Intel oneAPI 2021+ (icx and ifx)
* NVHPC 21.5+ (nvc++ and nvfortran)

Building from Source
====================

Clone the Repository
--------------------

.. code-block:: bash

   git clone https://github.com/zcobell/fdate.git
   cd fdate

Configure with CMake
--------------------

FortranDate uses CMake to manage the build process. The basic build configuration is:

.. code-block:: bash

   mkdir build
   cd build
   cmake ..

Advanced CMake Options
----------------------

There are several CMake options that control the build:

+----------------------+------------------+-------------------------------------------+
| Option               | Default          | Description                               |
+======================+==================+===========================================+
| FDATE_CXX_STANDARD   | 20               | C++ standard to use (17 or 20)            |
+----------------------+------------------+-------------------------------------------+
| FDATE_BUILD_SHARED   | OFF              | Build as shared library                   |
+----------------------+------------------+-------------------------------------------+
| FDATE_ENABLE_TESTING | OFF              | Build and run tests                       |
+----------------------+------------------+-------------------------------------------+

For example, to build a shared library with C++17 and enable testing:

.. code-block:: bash

   cmake -DFDATE_CXX_STANDARD=17 -DFDATE_BUILD_SHARED=ON -DFDATE_ENABLE_TESTING=ON ..

Build the Library
-----------------

.. code-block:: bash

   cmake --build .

Install (Optional)
------------------

You can install the library system-wide with:

.. code-block:: bash

   cmake --install .

By default, this will install to `/usr/local`. To specify a different installation prefix:

.. code-block:: bash

   cmake --install . --prefix /path/to/install

Testing (Optional)
------------------

If you enabled testing with ``FDATE_ENABLE_TESTING=ON``, you can run the tests with:

.. code-block:: bash

   ctest

Linking with Your Project
=========================

CMake Integration
-----------------

If you're using CMake for your project, the simplest way to use FortranDate is:

.. code-block:: cmake

   # In your CMakeLists.txt
   find_package(FortranDate REQUIRED)
   
   add_executable(your_program main.f90)
   target_link_libraries(your_program PRIVATE FortranDate::fortrandate)

Manual Linking
--------------

If you're not using CMake, you'll need to:

1. Add the include directory to your compiler flags
2. Link against the libraries

.. code-block:: bash

   # Compilation
   gfortran -c main.f90 -I/path/to/fortrandate/include
   
   # Linking
   gfortran main.o -o program -L/path/to/fortrandate/lib -lfortrandate

Troubleshooting
===============

Common Issues
-------------

**Linking errors with C++ standard library:**

Ensure your Fortran compiler is compatible with the C++ compiler used to build the library. For GCC/GFortran, they should be from the same version series.

**Runtime errors about missing symbols:**

If using a shared library build, ensure the library is in your ``LD_LIBRARY_PATH`` (Linux/Unix) or ``PATH`` (Windows).
