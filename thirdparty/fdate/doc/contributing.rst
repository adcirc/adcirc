============
Contributing
============

Thank you for your interest in contributing to FortranDate! This guide will help you get started with contributing to the project.

Getting Started
===============

Setting Up Your Development Environment
---------------------------------------

1. **Fork the Repository**

   Start by forking the FortranDate repository on GitHub.

2. **Clone Your Fork**

   .. code-block:: bash

      git clone --recursive https://github.com/your-username/fortrandate.git
      cd fortrandate

3. **Set Up the Build Environment**

   .. code-block:: bash

      mkdir build
      cd build
      cmake -DFDATE_ENABLE_TESTING=ON ..
      cmake --build .

4. **Run the Tests**

   .. code-block:: bash

      ctest

   All tests should pass before you start development.

Development Workflow
====================

Branching Strategy
------------------

* `main` - The production branch
* `develop` - The development branch where features are integrated
* Feature branches - Named as `feature/descriptive-name`
* Bugfix branches - Named as `bugfix/issue-number-description`

Creating a New Feature
----------------------

1. Create a new branch from `develop`:

   .. code-block:: bash

      git checkout develop
      git pull origin develop
      git checkout -b feature/your-feature-name

2. Implement your changes, including tests.

3. Run the tests locally to ensure they pass.

4. Commit your changes with a descriptive commit message:

   .. code-block:: bash

      git commit -m "Add feature: Description of the feature"

5. Push to your fork:

   .. code-block:: bash

      git push origin feature/your-feature-name

6. Create a pull request from your branch to the `develop` branch of the main repository.

Fixing a Bug
------------

1. Create a new branch from `develop`:

   .. code-block:: bash

      git checkout develop
      git pull origin develop
      git checkout -b bugfix/issue-number-description

2. Fix the bug and add tests to ensure it doesn't recur.

3. Run the tests locally.

4. Commit your changes:

   .. code-block:: bash

      git commit -m "Fix #123: Description of the bug fix"

5. Push to your fork:

   .. code-block:: bash

      git push origin bugfix/issue-number-description

6. Create a pull request from your branch to the `develop` branch.

Coding Standards
================

C++ Standards
-------------

* Follow the C++17/20 standard features
* Use modern C++ idioms (RAII, auto, lambdas, etc. where appropriate)
* Use `auto` when it improves readability
* Prefer compile-time constants and computations
* Use `constexpr` where possible
* Avoid raw pointers; use references, smart pointers, or optional
* Write exception-safe code

.. code-block:: cpp

   // Good example
   constexpr auto add(const int a, const int b) noexcept -> int {
       return a + b;
   }

   // Not preferred
   int add(int a, int b) {
       return a + b;
   }

Fortran Standards
-----------------

* Follow the Fortran 2008 standard
* Use `implicit none` in all modules and procedures
* Use modules instead of common blocks
* Use derived types for encapsulation
* Use parameter attributes for constants
* Use intent attributes for all dummy arguments
* Use pure or elemental procedures where appropriate

.. code-block:: fortran

   ! Good example
   pure function add(a, b) result(res)
      implicit none
      integer, intent(in) :: a, b
      integer :: res
      
      res = a + b
   end function add
   
   ! Not preferred
   function add(a, b)
      integer :: a, b, add
      
      add = a + b
   end function add

Documentation Standards
=======================

* Use Doxygen-style comments in C++ code
* Use FORD-compatible comments in Fortran code
* Document all public functions, types, and methods
* Include examples in documentation
* Update documentation when changing APIs

Testing Standards
=================

* All new features must have tests
* All bug fixes must have tests that would have caught the bug
* Tests should be both unit and integration level
* Aim for high code coverage

C++ Tests
---------

* Use Google Test framework for C++ tests
* Test both normal and edge cases
* Test error conditions where applicable

Fortran Tests
-------------

* Use pFUnit for Fortran tests
* Test all public procedures and interfaces
* Test with various input combinations

Pull Request Process
====================

1. Ensure all tests pass
2. Update the documentation if needed
3. Add a description of your changes
4. Reference any related issues
5. Request a review from one of the maintainers

Your pull request will be reviewed by a maintainer who may request changes or provide feedback. Once approved, your changes will be merged.

Reporting Bugs
==============

When reporting bugs, please include:

* Description of the issue
* Steps to reproduce
* Expected behavior
* Actual behavior
* Your environment:
  * Operating system
  * Compiler versions
  * CMake version
* Any relevant logs or output

Feature Requests
================

When requesting new features, please include:

* A description of the feature
* Why it would be valuable
* Any relevant use cases
* If possible, a proposed API design

Community
=========

Join the community to discuss development, ask questions, or provide feedback:

* GitHub Discussions
* Mailing List
* Chat Channel

Code of Conduct
===============

We follow a code of conduct to ensure a positive and inclusive environment for all contributors. Please be respectful and professional in all interactions.

Thank you for contributing to FortranDate!
