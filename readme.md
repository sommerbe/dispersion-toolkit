# Dispersion toolkit

## Introduction

### Purpose

* compute dispersion of a given point set
* compute empty boxes (interior, exterior, all) of a given point set
* optimise a point set w.r.t. minimising its dispersion
* visualise a sequence of a point set

### Usability

* easy to build, low maintenance => minimal dependencies
* executable in UNIX, Mac?, Windows?
* executable in server environments (console based) (compute intensive tasks)
* proper use of IO pipes and argument passing (standard in UNIX)
* interoperable with 3rd party toolkits: UTK (using scripts, file IO)

### Documentation

Each executable accepts the program argument

````
--help or -h
````

to print a very simple list of valid arguments. A detailed explanation is available in the secondary repository

````
dispersion-toolkit-manpages.git
````

### Dependencies

The least amount of requirements in order to build this toolkit to use it are:

* cmake 3.4 or more recent (tested until 3.19)
* C++11
* visualisation (optional): python3, numpy, matplotlib
* mindispgs (optional): OpenMP


## Getting started

### Building on UNIX systems

The following commands assume that the parent working directory equals this project's root directory.

Create an out-of-source build directory
````
mkdir build
cd build
````
and prepare the build using
````
cmake ../ -DCMAKE_BUILD_TYPE=Debug
````
or using

````
cmake ../ -DCMAKE_BUILD_TYPE=Release
````

Its output shows which dependencies are used, along with their current version. If either of the dependencies, GLFW or OpenGL, is missing or is expected to be incompatible, modules depending on these will be omitted from the build. Installing these dependencies following your operating system's standard practices should resolve this issue as long as compatibility is to be expected (technically).

In UNIX operating systems, a build parallelised among {P} threads is started with

````
make -j{P}
````

For instance,

````
make -j8
````

runs 8 builds in parallel.

### Windows systems

In Windows, create an out-of-source build directory and configure it using CMAKE. Build the source, usually with Visual Studio.


## Conventions for developers

Please refer to *conventions.md* within this directory.


## Source code organisation

The source code within the directory *./src* is organised as follows:

Directory | Computational capabilities
--------- | --------------------------
./adapter | transducing between various point set formats
./measure | algorithms computing dispersion measures
./opt     | optimising point sets w.r.t. dispersion measures
./set     | constructing d-dimensional point sets
./stat    | computing statistical descriptors
./vis     | visualising point sets, point set sequences or other results


## Program interoperability

Most programs of this toolkit interoperate using a *d-dimensional point set sequence*. This format is specified in the file *specifications.md* of this directory.

This approach is designed to work (well) in *headless* computing environments (cf. https://en.wikipedia.org/wiki/Headless_software), most notably computing clusters. These programs provide the building blocks in easily constructing more complex tasks through writing scripts, for instance Bash or Z shell.


## Exemplary usages

The following examples assume the parent working directy to be the above mentioned build directory, and that the build completed without errors.

**Dispersion of fibonacci lattice**

Compute dispersion of a Fibonacci lattice, a) by using IO pipes
````
./bin/fibonaccilattice --fibonacci-index=10 | ./bin/dispgs --disp
````
or b) by storing the lattice to the file ``pointset.dat`` and subsequently loading it:
````
./bin/fibonaccilattice --m=10 --o pointset.dat
./bin/dispgs --disp --i pointset.dat
````
or
````
./bin/fibonaccilattice --m=10 > pointset.dat
./bin/dispgs --disp --i pointset.dat
````
or append to an existing point set sequence
````
./bin/fibonaccilattice --m=10 >> pointset-sequence.dat
./bin/dispgs --disp --i pointset.dat
````
Usually, the order of arguments does not matter. In case the meaning of a program is unclear, or in case of wishing to find out available program parameters,
````
./bin/fibonaccilattice -h
````
shows a small description along with available options, and
````
./bin/fibonaccilattice --help
````
additionally shows the meaning the all options and further program requirements.


Compute point set's cardinality, n, multiplied by dispersion of a Fibonacci lattice:

````
./bin/fibonaccilattice --m=10 | ./bin/dispgs --ndisp
````

**Minimise dispersion of a fibonacci lattice**

Optimise a Fibonacci lattice w.r.t. minimising dispersion using gradient ascent and obtain its dispersion:
````
./bin/fibonaccilattice --fibonacci-index=10 | ./bin/mindispgs --tau=2e-15 --stepsize=0.01 --iteration-limit=10000 | ./bin/dispgs --ndisp
````
or with using files:
````
./bin/fibonaccilattice --m=10 --o pts-fibonacci.dat
./bin/mindispgs --tau=2e-15 --stepsize=0.01 --iteration-limit=10000 --i pts-fibonacci-m10.dat --o pts-min-fibonacci-m10.dat
./bin/dispgs --ndisp --i pts-min-fibonacci-m10.dat
````
Storing intermediate results in files is recommended since it keeps parameters and comments of each program or process. In this case, a storing these commands in scripts is recommended, for instance a bash script. IO pipes are supported for convenience, and for this toolkit to be in line with the UNIX system.

**Randomise a lattice and estimate a statistic of dispersion**

Randomise the Fibonacci lattice using coordinate swapping to emit a point set sequence, compute dispersion of each point set and estimate the inter quartile range statistics with upper and lower whiskers used to generate statistical box plots along with the arithmetic mean:
````
./bin/fibonaccilattice --m=10 | ./bin/cswap --count=1 --repeat=512 | ./bin/dispgs --ndisp | ./bin/confidence --iqr-box --mean
````
or with using files:
````
./bin/fibonaccilattice --m=10 --o pts-fibonacci-m10.dat
./bin/cswap --count=1 --repeat=512 --i pts-fibonacci-m10.dat --o pts-cswap.dat
./bin/dispgs --ndisp --i pts-cswap.dat --o disp-cswap.dat
./bin/confidence --iqr-box --mean --i disp-cswap.dat
````

**Estimate minimal dispersion of a randomized lattice**

In addition, try to minimise dispersion:
````
./bin/fibonaccilattice --m=10 | ./bin/cswap --count=1 --repeat=512 | ./bin/mindispgs --tau=2e-15 --stepsize=0.01 --iteration-limit=10000 | ./bin/dispgs --ndisp | ./bin/confidence --iqr-box --mean
````
or with using files:
````
./bin/fibonaccilattice --m=10 --o pts-fibonacci-m10.dat
./bin/cswap --count=1 --repeat=512 --i pts-fibonacci-m10.dat --o pts-cswap.dat
./bin/mindispgs --tau=2e-15 --stepsize=0.01 --iteration-limit=10000 --i pts-cswap.dat --o pts-cswap-mindisp.dat
./bin/dispgs --ndisp --i pts-cswap-mindisp.dat --o disp-cswap-ndisp.dat
./bin/confidence --iqr-box --mean --i disp-cswap-ndisp.dat
````

Notice that mindispgs retrieves a point set sequence greater or equal to 1. If the cmake build configuration found OpenMP, mindispgs optimises each point set with maximum parallelism supported by the actual hardware. Although being optional, using multi threading is highly recommended.

**Visualise minimisation of dispersion**

Visualise how points are moved by the dispersion optimisation:

````
./bin/fibonaccilattice --m=10 | ./bin/mindispgs --tau=2e-15 --stepsize=0.01 --iteration-limit=10000 --pointset-sequence | python ./bin/pss.py
````
or using files:
````
./bin/fibonaccilattice --m=10 --o pts-fibonacci-m10.dat
./bin/mindispgs --tau=2e-15 --stepsize=0.01 --iteration-limit=10000 --i pts-fibonacci-m10.dat --o pts-fibonacci-mindisp.dat
cat pts-fibonacci-mindisp.dat | python ./bin/pss.py
````
where ``cat`` is a UNIX or at least Linux system program to read files to stdin.

While mindispgs handles point set sequences, using ``--pointset-sequence`` to emit point sets during the gradient ascent would result in a sequence of point set sequences, being not supported by pss.py. To be precise, this stream would be equivalent to a longer point set sequence, in which previous point set sequences are stacked after each other in order. Therefore at some frame, pss.py would show points of an unoptimised set.

During this visualisation, each frame may be exported to permanent storage:

````
cat pts-fibonacci-mindisp.dat | python ./bin/pss.py --image-path='seq-{i}.png' --image-ppi=300
````

The result is a sequence of images, 
````
(seq-0.png, seq-1.png, ..., seq-n.png)
````

### Recommendation

Although IO piping is a powerful feature, offering a wide range of flexibility, a long command with many pipes might be difficult to follow, and at least to debug. Instead, keeping the number of pipes small while splitting the command into multiple commands, and using a script, for instance a bash script, increases both readability and detail of the protocol, academic or technical. For instance, mindispgs generates a log of how the gradient ascent progresses. This information may become important. However, these comment like logs are lost as soon as its entire output is piped into another program of this toolkit, such as dispgs.

## License

See the file *license* within the directory in which this file resides.

## Contributors

The file *credits* contains a list of those people who contributed throughout this project's lifetime. For those who (actively) maintain this project, or parts of it, please refer to the file *maintainers*.