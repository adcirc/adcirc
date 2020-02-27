# ADCIRC NUOPC/ESMF fork
ADCIRC NUOPC cap to couple into NOAA NEMS applications

### To Cite:
# Development of a Flexible Coupling Interface for ADCIRC Model for Coastal Inundation Studies

## Authors
Saeed Moghimi, Sergey Vinogradov, Edward P Myers, Yuji Funakoshi, Andre J Van der Westhuysen, Ali Abdolali, Zaizhong Ma, Fei Liu
## Publication date
2019
## Description
To enable flexible model coupling in storm surge studies, a coupling cap for ADvanced CIRCulation model (ADCIRC) was developed. The cap is essentially a wrap-around ADCIRC model which enables the model to communicate seamlessly with other model components, eg, surface wave and numerical weather prediction models. All the model components advertise their imported and exported fields at the runtime and connect to each other for exchanging data based on the availability of the advertised fields. Models can operate on structured or unstructured grids and the regridding capability will be provided by Earth System Modelling Framework (ESMF) and National Unified Operational Prediction Capability (NUOPC) infrastructures. We implemented the coupled application including ADCIRC cap as well as NUOPC compliant caps to read WaveWatchIII and Hurricane Weather Research and Forecasting Model (HWRF) generated forcing fields. We validated the coupled application for hurricane Ike on very high resolution mesh that covers the entire US Atlantic coastal water. We also showed that inclusion of the surface waves improves the model performance of both total water level and coastal inundation. Also shown how the maximum wave set-up and maximum surge regions may happen at various time and locations depending on the storm track and its landfalling region.

https://repository.library.noaa.gov/view/noaa/20609/noaa_20609_DS1.pdf






# The ADCIRC Developers Guide

[![CircleCI](https://circleci.com/gh/adcirc/adcirc-cg.svg?style=shield&circle-token=b93f1b49c76ccd6cbe1f016b0f4c9c0d86a44734)](https://circleci.com/gh/adcirc/adcirc-cg)

#### Authors
  * Jason Fleming, <jason.fleming@seahorsecoastal.com>
  * Zach Cobell, <zachary.cobell@arcadis.com>

## Development Strategy

The ADCIRC development strategy consists of the maintenance of two separate versions at any given time: the stable version and the development version. The stable version only receives bug fixes, while the development version receives bug fixes as well as new features.

At some point, the development version becomes stable. At that time, a stable branch is created from it, and then new features are added to the previously stable code, creating a new development version. A tag will be applied to the development version of the code in the form ```v49-dev```. Stable versions of the code will have a tag applied in the form ```v49.xx``` where ```xx``` represents the minor update revisions applied to fix bugs in the release versions.

In prior development of the code, ```version.F``` was used to track versioning. However, this is no longer necessary. ```version.F``` will reflect the current development state of the code in the development branch and the current stable version number of the code in the stable branch. There should be no editing necessary.

Similarly, ```header.F``` is no longer needed as these Git will automatically track the version number and log message. The ```git log``` command will provide you with this information in a more understandable way. *Please see the section on providing descriptive commit messages.* Release packages of the code will contain this information generated from Git internally so that users can view the commit log.

### Coding Standards

ADCIRC has had many individual contributors and has received code accretions over many years. A set of uniform coding standards has not been defined, and as a result, ADCIRC contains many different styles of Fortran. This section provides a set of style guidelines for contributing code to ADCIRC.

#### The Basics

* ```IMPLICIT NONE``` at the beginning of each subroutine.
* Fixed form, never exceeding the 72 column limit, even for comment
lines.
* When adding code to an existing subroutine, make the new code
match the style of the surrounding code, even if you'd prefer
another style.

#### Maintainability


When adding code that will be used in slightly different ways in different contexts, make it a subroutine, rather than cutting and pasting several times and making small changes to each cut-and-pasted section. Although it is faster to write code with cut-and-paste, the resulting code is harder to maintain, since each cut-and-pasted section will have to be modified individually later. Also, it is easier to make mistakes when working with cut-and-pasted code.

Rick has expressed a desire for greater modularity ... particularly with the number of variables in the global module. When adding a major new feature, please consider the modularity of the data it requires. In other words, if new variables can be made private to a particular module rather than global, please do so.

### Release Process

The release process informal, but generally includes the following steps:

* Set and stick to a list of new features as the goal for the new release.
* Have a consensus among ADCIRC developers that the goal has been achieved.
* Run tests that cover the new features as well as making sure old features still work.
* Fix issues revealed by failed tests.
* Distribute the release candidate code to collaborative users to make sure their runs perform as expected on the new code.
* Fix issues revealed by user dismay.
* Update the docs to reflect the changes. If input files or file formats have changed, produce release notes to highlight the changes.
* Publicly announce, release, and brag about a shiny new version of ADCIRC.

These tests are often referred to as "regression tests," i.e., a means to detect a regression in correctness or functionality due to changes in the code. A regression test suite can be built up pretty easily.  There are a number of small test cases on adcirc.org, so the initial set of tests should consist of those. We may also want to differentiate between tests for correctness and tests for functionality. The former is making sure you're getting the expected model results and the latter is making sure that nothing is broken (i.e., adcprep, i/o, message passing, etc).

The test suite gets built up based on different purposes.  One may want to exercise more options and more cases, but it is also important to test for known bugs that might have been reintroduced. In otherwords, once a bug is discovered and fixed, a test should be created specifically for that bug to make sure it doesn't appear again.

And lastly, a small test suite should be included with each distribution and a more comprehensive should be created that runs during development and before a release. A step in this direction has been made using ant in the autotest directory.

## Git/Github

### TL;DR

One time configuration steps:
```
[create your own Fork of the repository]
git clone https://github.com/YOUR_USERNAME/adcirc-cg.git
git remote add upstream https://github.com/adcirc/adcirc-cg.git
git config --global core.editor "vim"
git config commit.template .commit-template
```

General Work Flow:
```
git checkout myBranch

---Make Edits to Code---

git add [edited files]
git commit
git push
```

Note: In some cases you may need to be more specific with your checkout:
```
git checkout --track origin/myBranch
```

### Introduction

Git is a version control system. That is, it provides a means for many software developers to work on a  project while maintaining a historical record of changes to source code.

Github is a website built to enable collaboration between developers and provides many utilities for creating easy to use code repositories.

A Git repository is a storehouse of code for a particular project. Git repositories consist of a number of branches. The ```master``` branch is the mainline code that everyone shares and is considered the development version of the code. Another separate branch is maintained to house the current release version of the code.

Unlike other version control tools, Git creates commit numbers based upon a SHA1 hash. The SHA1 hash is a unique identifier that is completely unique to a specific commit as well as its history.

Git is also a fully distributed verrsion control system. When a repository is cloned, not only is the latest checked out commit downloaded, but also the entire history of the repository. This allows users to use the repository fully without an internet connection.

A commit ID looks like the following:

```
04b3c1f526b352f98203d7cf8ffe73b41657a4fe
```

These IDs are so unique that generally you only need the first 7 characters to fully describe a commit:
```
04b3c1f
```
Git calls these references "commitish".

### Commits

Each change set that is made to the code is called a "commit". Commits should be made regularly and should help tell a story of your changes. This is useful for other developers to understand what was done to the code and to track down the changes which may have caused a potential problem later. The greater the change in each commit, the more difficult it becomes to fix issues later.

#### Commit Messages
Your commit messages should be descriptive and use the template provided by this repository. Commit messages consist of two parts. A short description and a long description. The template for commits that is strongly suggested is as follows:

```
# ==[ Subject: One line only short meaningful description for logs ]===|


# ==[ Details: Describe what changed and explain why it changed ]======|


# ==[ Please wrap at 72 characters ]===================================|
```
This template is included in the repository and can be automatically used each time after running the following command:

```
git config commit.template .commit-template
```

Commits happen in two stages. First, the files that should be committed are added to the commit using ```git add```. For example:
```
git add timestep.F adcirc.F
```
This places the files in a staging area. Note that any further changes made to the files will mean they need to be added to the staging area again. You can check the current status of the repository using ```git status```.

To create the commit and add it to your current branch, use the ```git commit``` command. Note that git provides the ```-m``` option to include your message on the command line, however, this does not allow you to include the detail section of the commit message.

You can choose which text editor is used by git to craft your commit message using:
```
git config --global core.editor "vim"
```
This command would use the VIM editor, however, you can set it to any editor installed on your system.

### Branches

Branches are core to the way Git functions and should be used whenever possible. For the purposes of this repository, no work should be done in the ```master``` branch. A small number of branches will be officially supported in the main repository, however most work should be conducted within Forks, which are explained later in this document.

You can retrieve a list of branches currently configured in your repository:
```
git branch --list
```
If you'd like to retrieve a branch not currently set up in your system, you can list them using:
```
git branch --list -a
```
which will display something like:
```
remotes/origin/GAHM_jie
remotes/origin/GPGPU-CUDA
remotes/origin/HEAD -> origin/master
remotes/origin/SSD_jie
```
To checkout the GPGPU-CUDA branch, simply type:

```
git checkout GPGPU-CUDA
```

Note: In some cases you may need to be more specific with your checkout:
```
git checkout --track origin/myBranch
```

Branches can be created using the command:
```
git checkout -b myNewFeature
```

Finally, new branches can be pushed to the server using:
```
git push -u origin myNewFeature
```
Note that the ```-u origin``` command only needs to be done the first time a new branch is pushed.


### Forks

Forks are fully fledged repositories that are available for development work. Forks allow developers to create as many branches as they would like without the approval of repository administrators to conduct experiments and develop ideas. A forked repository is registered to a particular user but maintains an association with the master ADCIRC repository without cluttering it.

To create your own Fork, use the Fork button from the ```adcirc-cg``` repository page on Github. When you do this, a new repository will be created under your account. This repository is your personal playground.

Note that the ADCIRC repository is private. Your local repository will be private as well.

You can obtain a copy of your repository using the following command:

```
git clone https://github.com/YOUR_USERNAME/adcirc-cg.git
```

Next, you can configure your personal repository to maintain an association with the master ADCIRC repository (known as the "Upstream") using:

```
git remote add upstream https://github.com/adcirc/adcirc-cg.git
```
[More Information about Forks](https://help.github.com/articles/configuring-a-remote-for-a-fork/)

You can retrieve changes from the upstream repository using:
```
git fetch upstream
```
To incorporate these changes into your repository, you can use:

```
git merge upstream/master
```
It is critically important that you make sure you are on the correct branch that you want to merge changes into when running the above command.
[More Information about Syncing Forks](https://help.github.com/articles/syncing-a-fork/)

### Submitting Changes to the Upstream

Changes are accepted into the ADCIRC repository via a Pull Request. From the Github web page from your Forked repository, click the "Compare and Pull Request" button. You will be asked to select a ```base``` and a ```compare```. The base is the target for your changes, i.e. the master branch in the upstream repository. The compare is your local branch.

Warning: Your local branch must contain all changes from the upstream master to be accepted!

After clicking the pull request button, you will be given a section to describe your changes. Please describe what has been done throughly. The repository administrators will be given the option to merge your pull request into the upstream repository.

### Continuous Integration (CI)

Continuous Integration, or CI, involves testing each change made to the model against a known result. The ensures that as soon as a change (accidental or intentional) to the solution occurs, the developers are aware of it.

The ADCIRC repository uses the CircleCI.com service to conduct these tests. When a pull request is submitted on Github, the CI server will do the following:

1. Build the code without netCDF enabled
2. Build the code with netCDF enabled
3. Run the test suite found [here](https://github.com/adcirc/adcirc-cg-testsuite). The test suite consists of the following derived from examples found on ADCIRC.org:
    1. Quarter Annular 2D example (serial and parallel)
    2. Shinnecock Inlet example (serial and parallel)
    3. Idealized Inlet example (serial and parallel)
    4. APES Pamlico Sound example (serial and parallel)

The CI server is controlled via the ```circle.yml``` file found in the root directory. It contains the instructions to prepare the build server. The actual test criteria is maintained directly in the test suite repository. Keeping the repositories separate is important so that the code repository does not become bloated.

### General Policies

* Communicate. Jason Fleming is responsible for making sure ADCIRC
development is smooth, pain-free and productive---when in doubt,
email him (jason.fleming@seahorsecoastal.com).
* Also use the adcirc-dev mailing list to keep everyone informed of
what you are doing.
* If you are not on the adcirc-dev mailing list and would like to be,
email Jason Fleming.
* Make a fork to develop new features. No commits are allowed to the repository without a pull request.

### Important Git Tools
The following are some tools available in the Git suite that are extremely useful.

#### Git Bisect
Git Bisect allows users to determine when a particular issue occurred. In its most basic form, it can track the origin of a bug very quickly using a binary search. To use it, do the following:

1. Start the tool:

    ```
    git bisect start
    ```
2. Inform the tool of the last time the code worked:

    ```
    git bisect good [commit-id]
    ```
3. Inform the code of a version of the code that does not work:

    ```
    git bisect bad [commit-id]
    ```
4. Git will then check out revisions of the code and ask you to let it know if it works or does not work. When you've determined if a particular version works or does not work, type ```git bisect good``` or ```git bisect bad```. When Git has determined when the issue started, it will give you the details of that commit.

[More Git Bisect Documentation](https://git-scm.com/docs/git-bisect)

#### Git Diff
Git diff is used to difference commits or files within the repository. You can either difference individual files or entire commits.

```
git diff 4f5c3a HEAD
```

This will difference the specified revision with the current ```HEAD``` revision. ```HEAD``` is a reference to the currently checked out commit. You can also use short hand such as ```HEAD~3``` which means "3 commits before HEAD"

#### Git Describe
Git describe is used to reference a commit. It allows a more user friendly representation of the commit id using tags. For example, a branch with two commits since the ```v53-dev``` tag was applied would look something like:
```
% git describe --tags --always
v53-dev-2-g7fae56e
```
This is a unique reference to this commit.
