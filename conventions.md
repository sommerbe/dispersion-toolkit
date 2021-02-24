# Conventions for developers

## Versioning

All releases are versioned according to the *Backus–Naur Form of Semantic Versioning 2.0.0* to be found at https://semver.org/. Basically, this entails the version core to match

````
MAJOR.MINOR.PATCH
````

where

* MAJOR reflects a code change breaking backwards compatibility,
* MINOR reflects a code change with backwards compatibility while introducing new features,
* PATCH reflects a code change with error fixes only, and therefore remaining backwards compatible without introducing new features,

and MAJOR, MINOR, and PATCH are unsigned integers in increasing order.


## Repository size

The probability of erroneously asserting a backwards incompatibility, infered due to the MAJOR version number, with respect to downstream software needs to be reduced as much as possible. This goal is followed strongly with keeping the computational capabilities exposed by this software repository as small as possible. At the same time, the number of repositories of this dispersion toolkit should be relatively small to ease its compilation.

As this toolkit grows in computational capabilities, splitting it into separate repositories becomes mandatory, for instance with respect to point set generation, visualisation, and optimisation.