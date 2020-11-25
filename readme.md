# Dispersion toolkit

## Purpose

a) compute dispersion of a given point set
b) compute empty boxes (interior, exterior, all) of a given point set
c) optimise a point set w.r.t. minimising its dispersion
d) visualise a sequence of a point set using OpenGL while its creation

## Usability

a) easy to build, low maintenance => minimal dependencies
b) executable in UNIX, Mac?, Windows?
c) executable in server environments (console based) (compute intensive tasks)
d) helper scripts to generate publication-oriented figures
e) maybe: bash scripts to run many tasks in parallel
f) proper use of IO pipes and argument passing (standard in UNIX)
g) interoperable with 3rd party toolkits: UTK (using scripts, file IO)

## Dependencies

a) C++11
b) optional: OpenGL, GLFW3

## Building

The following commands assume that the current working directory equals this project's root directory.

Configure the local build environment by running

``
bash configure.sh
``

This pulls internal dependencies, configures them, and pre-builds them if necessary.

Create an out-of-source build directory

``
mkdir build
cd build
``

and prepare the build using

``
cmake ../ -DCMAKE_BUILD_TYPE=Debug
``

or using

``
cmake ../ -DCMAKE_BUILD_TYPE=Release
``

Its output shows which dependencies are used, along with their current version (along with a range of compatible version numbers). If either of the dependencies, GLFW or OpenGL, is missing or is expected to be incompatible, modules depending on these will be omitted from the build. Installing these dependencies following your operating system's standard practices should resolve this issue as long as compatibility is to be expected (technically).

In Linux operating systems, usually run

``
make -j{P}
``

with {P} equals the positive integer of parallel builds, for instance

``
make -j8
``

runs 8 builds in parallel.

In Windows, this build directory contains a Visual Studio solution file to opened.

# License

See the file *license.md* within the directory in which this file resides.

# Contributors

Benjamin Sommer: 2020 - current

