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

### Versioning

All releases are versioned according to

````
MAJOR.MINOR.FIX
````

where

* MAJOR reflects a code change **breaking backwards compatibility**,
* MINOR reflects a code change **with backwards compatibility while introducing new features**,
* FIX reflects a code change **with error fixes only**, and therefore remaining backwards compatible without introducing new features,

and MAJOR, MINOR, and FIX are unsigned integers in increasing order.


## File format of a point set sequence

The bounded domain is spanned by a tuple of smallest possible coordinates,
````
(low_0 low_1 ... low_n),
````
and a tuple of greatest possible coordinates
````
(up_0 up_1 ... up_n).
````

An n-dimensional point set P_i with cardinality k, having points p_j with coordinates c_{i,j,s} is serialised with ASCII encoding according to the format
````
#d low_0 low_1 ... low_n up_0 up_1 ... up_n
c_{i,0,0} c_{i,0,1} ... c_{i,0,n}
c_{i,1,0} c_{i,1,1} ... c_{i,1,n}
                     .
                     .
                     .
c_{i,k,0} c_{i,k,1} ... c_{i,k,n}
#eos
````
where columns are separated by an ASCII character named *delimiter*, here
````
' '.
````

Each line of the serialised result may additionally contain a commenting line,
````
# this is a comment.
````

A point set sequence of length m is the concatenation of serialised point sets, i.e.
````
#d low_0 low_1 ... low_n0 up_0 up_1 ... up_n0
c_{0,0,0}  c_{0,0,1}  ... c_{0,0,n0}
                       .
c_{0,k0,0} c_{0,k0,1} ... c_{0,k0,n0}
#eos
#d low_0 low_1 ... low_n1 up_0 up_1 ... up_n1
c_{1,0,0}  c_{1,0,1}  ... c_{1,0,n1}
                       .
c_{1,k1,0} c_{1,k1,1} ... c_{1,k1,n1}
#eos
                       .
                       .
#d low_0 low_1 ... low_nm up_0 up_1 ... up_nm
c_{m,0,0}  c_{m,0,1}  ... c_{m,0,nm}
                       .
c_{m,km,0} c_{m,km,1} ... c_{m,km,nm}
#eos
````


## Examples

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

Benjamin Sommer: 2020 - current

