# ADCIRC

ADCIRC is a system of computer programs for solving time dependent, free surface circulation and transport problems in
two and three dimensions. These programs utilize the finite element method in space allowing the use of highly flexible,
unstructured grids. Typical ADCIRC applications have included:

* Prediction of storm surge and flooding
* Modeling tides and wind driven circulation
* Larval transport studies
* Near shore marine operations
* Dredging feasibility and material disposal studies

ADCIRC has been used around the world for various studies, including those conducted by the United States Army Corps of
Engineers (USACE), Federal Emergency Management Agency (FEMA), National Oceanographic and Atmospheric Administration 
(NOAA), and many others.

## Authors
* Joannes Westerink, University of Notre Dame
* Rick Luettich, University of North Carolina at Chapel Hill

### Development Group
* Brian Blanton - RENCI
* Zachary Cobell - The Water Institute of the Gulf
* Clint Dawson - University of Texas at Austin
* Casey Dietrich - North Carolina State University
* Randall Kolar - University of Oklahoma at Norman
* Chris Massey - US Army Corps of Engineers Research and Development Center, Coastal and Hydraulics Laboratory

## CI Status
The current main branch is tested on CircleCI

[![CircleCI](https://dl.circleci.com/status-badge/img/gh/adcirc/adcirc/tree/main.svg?style=svg&circle-token=468312e3a9341f3a519bbdfb4df0cda07c98bd91)](https://dl.circleci.com/status-badge/redirect/gh/adcirc/adcirc/tree/main)

# Gallery

Louisiana ADCIRC model simulating Hurricane Katrina storm surge and waves developed by The Water Institute of the Gulf.

<img src="https://i0.wp.com/www.psc.edu/wp-content/uploads/2021/07/katrina_aws_z0_0250-scaled.jpg?resize=1080%2C448&ssl=1" alt="Hurricane Katrina Simulation" width="800"/>

ADCIRC mesh in the Chesapeake Bay area used for the FEMA Coastal Storm Surge Study

<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/8/89/FEMA_Region_III_Coastal_Storm_Surge_Study_%28page_7_crop%29.jpg/1280px-FEMA_Region_III_Coastal_Storm_Surge_Study_%28page_7_crop%29.jpg" alt="ADCIRC mesh in the Chesapeake Bay area used for the FEMA Coastal Storm Surge Study" width="400"/>

## Usage

Conditions for obtaining the code are the following

1. That it is used for research or education purposes,
2. That ADCIRC is identified by name in publications, and
3. That individuals do not redistribute ADCIRC but rather have interested parties email Crystal Fulcher to obtain the
   code

# Versions

Code versions are published based approximately on semantic versioning. Using the major version number (i.e. XX in the
version XX.YY) will maintain solution consistency across minor version numbers except when there are critical bug fixes.
Changes to the major version do not guarantee solution consistency, which may be due to improvements in the algorithm or
other fixes. The general opinion is that the solution is improved in greater major version numbers and if exact
consistency is required, it's recommended that the same major version is used. Note that other factors, including
compiler versions and optimization may also impact solution consistency and the user should understand their compiler.

## Documentation

Documentation is presently undergoing upgrades, however, the main documentation locations for users are:

1. ADCIRC website - https://adcirc.org
2. ADCIRC Wiki - https://wiki.adcirc.org/Main_Page

## Examples

The ADCIRC testing repository (http://github.com/adcirc/adcirc-testsuite) doubles as a set of examples which can be used
for new users to become acquainted with the model. Since version 55, the branches are annotated with the expected
version numbers that would allow the tests to run successfully.
