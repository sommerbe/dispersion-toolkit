# Dispersion toolkit

## Purpose

* compute dispersion of a given point set
* compute empty boxes (interior, exterior, all) of a given point set
* optimise a point set w.r.t. minimising its dispersion
* visualise a sequence of a point set

## Documentation

A manpage exists for each executable. After build configuration, these are to be found in ./build/man/, for convenience, or in a chosen directory after installing the build, for instance {CMAKE_BUILD_PREFIX}/share/man/.

A manpage is located in the directory of each C++11 main.cpp source file.

## Usability

* easy to build, low maintenance => minimal dependencies
* executable in UNIX, Mac?, Windows?
* executable in server environments (console based) (compute intensive tasks)
* proper use of IO pipes and argument passing (standard in UNIX)
* interoperable with 3rd party toolkits: UTK (using scripts, file IO)

## Dependencies

### Users

* cmake 3.4 or more recent (tested until 3.19)
* C++11
* visualisation: python3, numpy, matplotlib
* optional (dispoptgs): OpenMP

### Developers

* manual generation: bash, pandoc

## Building

The following commands assume that the parent working directory equals this project's root directory.

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

Its output shows which dependencies are used, along with their current version. If either of the dependencies, GLFW or OpenGL, is missing or is expected to be incompatible, modules depending on these will be omitted from the build. Installing these dependencies following your operating system's standard practices should resolve this issue as long as compatibility is to be expected (technically).

### UNIX

In UNIX operating systems, a build parallelised among {P} threads is started with

``
make -j{P}
``

For instance,

``
make -j8
``

runs 8 builds in parallel.

### Windows

In Windows, this build directory contains a Visual Studio solution file to be opened.

## Examples

The following examples assume the parent working directy to be the above mentioned build directory, and that the build completed without errors.

Compute dispersion of a Fibonacci lattice:

``
./bin/fibonaccilattice --fibonacci-index=10 | ./bin/dispgs --disp
``

Compute point set's cardinality, n, multiplied by dispersion of a Fibonacci lattice:

``
./bin/fibonaccilattice --fibonacci-index=10 | ./bin/dispgs --ndisp
``

Optimise a Fibonacci lattice w.r.t. minimising dispersion using gradient ascent and obtain its dispersion:

``
./bin/fibonaccilattice --fibonacci-index=10 | ./bin/dispoptgs --tau=2e-15 --stepsize=0.01 --iteration-limit=10000 | ./bin/dispgs --ndisp
``

Randomise the Fibonacci lattice using coordinate swapping to emit a point set sequence, compute dispersion of each point set and estimate the inter quartile range statistics with upper and lower whiskers used to generate statistical box plots along with the arithmetic mean:

``
./bin/fibonaccilattice --fibonacci-index=10 | ./bin/cswap --count=1 --repeat=512 | ./bin/dispgs --ndisp | ./bin/confidence --iqr-box --mean
``

In addition, try to minimise dispersion:

``
./bin/fibonaccilattice --fibonacci-index=10 | ./bin/cswap --count=1 --repeat=512 | ./bin/dispoptgs --tau=2e-15 --stepsize=0.01 --iteration-limit=10000 | ./bin/dispgs --ndisp | ./bin/confidence --iqr-box --mean
``

Notice that dispoptgs retrieves a point set sequence greater 1. If the cmake build configuration found OpenMP, dispoptgs optimises each point set with maximum parallelism supported by the actual hardware. Although being optional, using multi threading is highly recommended.

Visualise how points are moved by the dispersion optimisation:

``
./bin/fibonaccilattice --fibonacci-index=10 | ./bin/dispoptgs --tau=2e-15 --stepsize=0.01 --iteration-limit=10000 --pointset-sequence | python ./bin/pss.py
``

While dispoptgs handles point set sequences, using --pointset-sequence to emit point sets during the gradient ascent would result in a sequence of point set sequences, being not supported by pss.py. To be precise, this stream would be equivalent to a longer point set sequence, in which previous point set sequences are stacked after each other in order. Therefore at some frame, pss.py would show points of an unoptimised set.

During this visualisation, each frame may be exported to permanent storage:

``
./bin/fibonaccilattice --fibonacci-index=10 | ./bin/dispoptgs --tau=2e-15 --stepsize=0.01 --iteration-limit=10000 --pointset-sequence | python ./bin/pss.py --image-path='seq-{i}.png' --image-ppi=300
``

The result is a sequence of images, (seq-{0}.png, seq-{1}.png, ..., seq-{n}.png).

### Recommendation

Although IO piping is a powerful feature, offering a wide range of flexibility, a long command with many pipes might be difficult to follow, and at least to debug. Instead, keeping the number of pipes small while splitting the command into multiple commands, and using a script, for instance a bash script, increases both readability and detail of the protocol, academic or technical. For instance, dispoptgs generates a log of how the gradient ascent progresses. This information may become important. However, these comment like logs are lost as soon as its entire output is piped into another program of this toolkit, such as dispgs.

## License

See the file *license.md* within the directory in which this file resides.

## Contributors

Benjamin Sommer: 2020 - current

