============
Contributing
============

Authors
"""""""
* Joannes Westerink
* Rick Luettich

Maintainer
""""""""""
* Zach Cobell, <zcobell@thewaterinstitute.org>

Development Strategy
====================

The ADCIRC development strategy consists of the maintenance of two separate versions at any given time:
the stable version and the development version. The stable version only receives bug fixes, while the
development version receives bug fixes as well as new features.

At some point, the development version becomes stable. At that time, a stable branch is created from it,
and then new features are added to the previously stable code, creating a new development version.
A tag will be applied to the development version of the code in the form ``v49-dev`` in order to act
as a marker for the start of development of a new version. Stable versions of the code will have a tag
applied using semantic versioning.

Branching Strategy
------------------

Apart from work branches, the branching strategy is as follows:

* ``main`` - the main branch of the repository. This branch should always be in a working state.
  It considered the bleeding edge of development
* ``vXXRelease`` - the stable branch for version XX of the code. This branch should only receive bug fixes.

We ask that contributors name their branches in the form:
``<your_username>/<feature_name>``, e.g. ``zcobell/bugfix-1234``

Considerations
--------------

Please consider the following when conducting development on the ADCIRC code base:

1. You should regularly keep your code up to date with the main branch. To do this, you should rebase from the
   main branch. This will also ensure a linear history. Use the following command to rebase from the main branch:

   .. code-block:: bash

      git fetch origin  # Assuming origin is the remote name for https://github.com/adcirc/adcirc
      git rebase origin/main

   Code must have a linear history in order to be merged into the main branch. If your branch becomes extremely out
   of date with the main branch, you will likely encounter merge conflicts here, and these conflicts are best resolved
   by the developer submitting the pull request than the maintainer of the main branch.

2. When you are ready to merge your branch into the main branch, you should use a pull request. This will allow
   the code to be reviewed by other developers before it is merged into the main branch. This will also allow the
   code to be tested by the continuous integration system.

3. If you are submitting a bug fix, you should also submit a pull request to the stable branch (i.e. ``vXXRelease``)
   once it is accepted into the main branch. This will allow the bug fix to be incorporated into the next minor release.

4. Pull requests should be submitted against the ``main`` branch. Prior ADCIRC development has maintained a ``main`` and
   ``development`` branch, however, in practical use, it seems that only the ``main`` branch in addition to the release branches
   are necessary.

5. Pull requests for new features will not be accepted without an acceptance test. This test can be included in the pull
   request and should be briefly described in the pull request description, including what the expected result should be. The
   acceptance test should run in a few minutes on a single core of a modern desktop computer. If the test takes longer than
   this, it will not be viable to run within the continuous integration system.

Coding Standards
----------------

ADCIRC has had many individual contributors and has received code accretions over many years. As such,
the code base has accrued many different styles. The goal of the standards below should be applied to newly
contributed code such that over time we can move in the direction of a preferred style.

File Naming
^^^^^^^^^^^
All new files should be named using a ".F90" extension where historically we have used a "*.F" extension.
There are multiple reasons for this. First, the code style in the older portions of the code are a mix of fixed and free form Fortran and don't adhere to
any specific standard. To move toward a more modern Fortran style and clearly communicate the
expectations, contributors should use the ".F90" extension. Second, this allows us to use the extension
as a marker for new code which can be subject to additional linting and compiler checks.

Fortran Standards
^^^^^^^^^^^^^^^^^
The ADCIRC code base has historically been written in a mix of Fortran 77 and Fortran 90, requiring
compiler flags to allow for the use of both. The goal of the ADCIRC development team is to move toward
more modern Fortran standards, and we specifically name Fortran 2008 as our target. However, the legacy
portions of the code are large, and it will take time to move toward this goal. As such, we will enforce
this guideline only for new code.

Documentation
^^^^^^^^^^^^^
All new code should be documented using Doxygen style comments. The benefits to this are the inline
documentation can be converted to a manual for developers in time and there is a defined style for
documentation that is easy to follow.

Pre-Commit
^^^^^^^^^^
A series of linting and formatting tools are run using the ``pre-commit`` package.
The rules are defined in the ``.pre-commit-config.yaml`` file in the root of the repository.

You can configure ``pre-commit`` to run the linter on every commit automatically before every commit by installing
the ``pre-commit`` hooks using:

.. code-block:: bash

   pre-commit install

Alternatively, you can run the linter manually using:

.. code-block:: bash

   pre-commit run --all

When a pull request is submitted, all pre-commit hooks will run and the pull request will be blocked
until all hooks pass.

Fortitude
"""""""""
We have developed a set of Fortran linting rules using the ``fortitude`` package. Fortitude is a Fortran linter
based on the Python Ruff linter. The rules can be enforced and some can be automatically fixed. Note that because
there are large amounts of legacy code, the linter will only run for "*.F90" files.

Some of the basic guidelines that Fortitude enforces are:

* Required intent for function arguments
* Required ``IMPLICIT NONE`` in all subroutines
* Required module default access specifier. Note that our own guideline is that ALL modules should be ``private`` by default.

Additional pre-commit hooks
"""""""""""""""""""""""""""

We also run additional hook on the code base such as:

* ``end-of-file-fixer`` - ensures that the file ends with a newline
* ``mixed-line-ending`` - ensures that the file uses the same line ending style throughout
* ``cmake-format`` - formats CMake files
* ``fprettify`` - formats Fortran files in a consistent style
* ``ruff-format`` - formats Python files in a consistent style

Compiler Warnings and Developer Mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Compilers can be configured to provide warnings about code that is not standard compliant,
does not follow best practices, may not do what you expect, or does not follow certain good
programming practices. The CMake build system has been configured to use have a "developer mode"
where the GNU compiler suite will compile with maximum warnings, all of which are treated as
errors. Note that we do not have the same level of strictness with third party code that is built
within ADCIRC since we do not control that code, however, the ADCIRC code should compile warning
free with the maximum set of warnings. These checks will be run as part of the continuous integration
system and will block the pull request until all warnings are resolved.

The warnings that are enabled in developer mode are:

.. code-block:: text

   -Werror
   -Wall
   -Wextra
   -pedantic
   -fimplicit-none
   -Wconversion
   -Wconversion-extra
   -Wuninitialized
   -Wsurprising
   -Wuse-without-only
   -Wimplicit-procedure
   -Winteger-division

You can enable developer mode when working with the GNU compilers and the CMake build system by
adding the following to your CMake command line:

.. code-block:: bash

   -DADCIRC_DEVELOPER_MODE=ON

Note that because there are many legacy portions of the code, these warnings are only enforced for
new "*.F90" files. The goal is that over time, we will slowly have more of the codebase passing these
strict warning checks.

Maintainability
^^^^^^^^^^^^^^^

It is very attractive to simply add new code either within a subroutine or function or
at the end or middle of an existing module. It is important to think about why you are
adding the code you are adding, and we encourage you to think about the modularity of any
change. In the general case, if you are adding a new behavior to the code for an existing
feature, this is likely a good candidate for a new subroutine or function. In the case that
you are adding a new feature, this is likely a good candidate for a new module.

In all cases, adding additional variables to the ``GLOBAL`` module is highly discouraged. By
adding variables to the ``GLOBAL`` module, you are more likely to create spooky action at a
distance, i.e., the values held in those variables can change due to things you can not
reason about in your local view of the code, and thus will be more challenging to understand.
Instead, you should  consider passing variables to subroutines and  functions as arguments.
Compilers are quite good at optimizing this and the performance impact is negligible.
Additionally, just because  the variables are not within the ``GLOBAL`` module, but in some other
module instead, does not mean that they are not global in their scope. We recommend that
when you need to maintain a state variable, you should consider using a data structure or
derived type rather than pointing at a variable in some global scope or other module.

Lastly, we ask developers to adhere to the DRY principle, i.e., Don't Repeat Yourself. If
you find yourself feeling the urge to copy and paste code, or copy and paste code even if
it requires a small to medium-sized change, that's typically a sign that you should be
creating a new subroutine or function. This will help to keep the code base
clean and maintainable. If you find yourself needing similar functionality in multiple modules,
consider that this is a good candidate for the start of a new module which might encompass
other pieces of functionality as well.

Please follow the good boy scouts rule and leave the code cleaner than you found it.

Testing
-------

ADCIRC has a long history of development and was developed long before modern coding practices of unit testing.
As a result, the code was not necessarily developed to be unit-testable. Therefore, ADCIRC uses a test method
approximately called "approval testing." ADCIRC runs a suite of tests located in the testing
`repository <https://github.com/adcirc/adcirc-testsuite>`__. These tests are run for any pull request to the code
and verify that the solution has not changed outside of minor differences due to compiler differences. The test
suite uses a docker container available on Docker Hub (``adcircorg/adcirc-ci:2025.0.2``).