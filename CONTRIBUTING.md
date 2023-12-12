# The ADCIRC Developers Guide

#### Authors
  * Joannes Westerink
  * Rick Luettich

#### Maintainer
  * Zach Cobell, <zcobell@thewaterinstitute.org>

## Development Strategy

The ADCIRC development strategy consists of the maintenance of two separate versions at any given time: 
the stable version and the development version. The stable version only receives bug fixes, while the 
development version receives bug fixes as well as new features.

At some point, the development version becomes stable. At that time, a stable branch is created from it, 
and then new features are added to the previously stable code, creating a new development version. 
A tag will be applied to the development version of the code in the form `v49-dev` in order to act
as a marker for the start of development of a new version. Stable versions of the code will have a tag 
applied using semantic versioning.


### Branching Strategy

Apart from work branches, the branching strategy is as follows:

* `main` - the main branch of the repository. This branch should always be in a working state. 
It considered the bleeding edge of development
* `vXXRelease` - the stable branch for version XX of the code. This branch should only receive bug fixes.

### Considerations

Please consider the following when conducting development on the ADCIRC code base:

1. You should regularly keep your code up to date with the main branch. To do this, you should rebase from the 
main branch. This will also ensure a linear history. Use the following command to rebase from the main branch:

```bash
git fetch origin  # Assuming origin is the remote name for https://github.com/adcirc/adcirc
git rebase origin/main
```
Code must have a linear history in order to be merged into the main branch. If your branch becomes extremely out 
of date with the main branch, you will likely encounter merge conflicts here, and these conflicts are best resolved
by the developer submitting the pull request than the maintainer of the main branch.
 
2. When you are ready to merge your branch into the main branch, you should use a pull request. This will allow
the code to be reviewed by other developers before it is merged into the main branch. This will also allow the
code to be tested by the continuous integration system.
3. If you are submitting a bug fix, you should also submit a pull request to the stable branch (i.e. `vXXRelease`) 
once it is accepted into the main branch. This will allow the bug fix to be incorporated into the next minor release.
3. Pull requests should be submitted against the `main` branch. Prior ADCIRC development has maintained a `main` and 
`development` branch, however, in practical use, it seems that only the `main` branch in addition to the release branches
are necessary.
4. Pull requests for new features will not be accepted without an acceptance test. This test can be included in the pull
request and should be briefly described in the pull request description, including what the expected result should be. The
acceptance test should run in a few minutes on a single core of a modern desktop computer. If the test takes longer than
this, it will not be viable to run within the continuous integration system.


### Coding Standards

ADCIRC has had many individual contributors and has received code accretions over many years. A set of 
uniform coding standards has not been defined, and as a result, ADCIRC contains many different styles 
of Fortran. This section provides a set of style guidelines for contributing code to ADCIRC.

#### The Basics

* `IMPLICIT NONE` at the beginning of each subroutine.
* When adding code to an existing subroutine, make the new code
match the style of the surrounding code, even if you'd prefer
another style.
* New source files should use F90-style free-form source so that 
the code can slowly be upgraded to modern F90 style.

#### Maintainability

When adding code that will be used in slightly different ways in different contexts, make it a subroutine, 
rather than cutting and pasting several times and making small changes to each cut-and-pasted section. 
Although it is faster to write code with cut-and-paste, the resulting code is harder to maintain, 
since each cut-and-pasted section will have to be modified individually later. Also, it is easier 
to make mistakes when working with cut-and-pasted code.

Rick has expressed a desire for greater modularity ... particularly with the number of variables 
in the global module. When adding a major new feature, please consider the modularity of the 
data it requires. In other words, if new variables can be made private to a particular module 
rather than global, please do so. Avoid adding variables to the global module and if possible
remove them from it. 

Please follow the good scouts rule and leave the code cleaner than you found it.

### Release Process

The release process is informal, but generally includes the following steps:

* Set and stick to a list of new features as the goal for the new release.
* Have a consensus among ADCIRC developers that the goal has been achieved.
* Run tests that cover the new features as well as making sure old features still work.
* Fix issues revealed by failed tests.
* Distribute the release candidate code to collaborative users to make sure their runs perform as expected on the new code.
* Fix issues revealed by user dismay.
* Update the docs to reflect the changes. If input files or file formats have changed, produce release notes to highlight the changes.
* Publicly announce and release a shiny new version of ADCIRC.

### Testing

ADCIRC has a long history of development and was developed long before modern coding practices of unit testing. 
As a result, the code was not necessarily developed to be unit-testable. Therefore, ADCIRC uses a test method
approximately called "approval testing." ADCIRC runs a suite of tests located in the testing 
[repository](https://github.com/adcirc/adcirc-testsuite). These tests are run for any pull request to the code
and verify that the solution has not changed outside of minor differences due to compiler differences. The test
suite uses a docker container available on Docker Hub (`zcobell/adcirc-ci-container:latest`) and can be build using
the Dockerfile in the CI docker repository [here](https://github.com/adcirc/adcirc-ci-docker).
