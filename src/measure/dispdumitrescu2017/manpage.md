% DISP-DUMITRESCU2017(1) 1.2.0 | Dispersion Toolkit Manuals
% Benjamin Sommer
% December 9, 2020

# NAME

disp-dumitrescu2017 - compute dispersion with algorithm of Dumitrescu and Jiang, 2017

# SYNOPSIS

**disp-dumitrescu2017** [**\--i** *FILE*] [**\--o** *FILE*] [**\--disp**] [**\--ndisp**] [**\--count-boxes**] [**\--silent**]

# DESCRIPTION

Computes dispersion, *n* * dispersion, and/or the number of RRs considered of a given point set with cardinality *n*, almost according to Dumitrescu and Jiang, 2017 (https://doi.org/10.1007%2Fs00454-017-9871-1). Instead of using the 2-dimensional algorithm of Aggarwal und Suri, 1987 (https://doi.org/10.1145%2F41958.41988) as described in the publication of Dumitrescu and Jiang, the implementation uses the simpler while computationally more expensive algorithm published by Naamad et. al. 1984 ().

Computational complexity: O(epsilon^(-2) n^(-2/3) n^2),

Memory complexity: O(n).

The measures are written to *standard output*, or to the file given by **\--o** *FILE*.

# OPTIONS

**\--i** *FILE*, **\--i**=*FILE*
:   Retrieves a point set sequence (P_0, P_1, ..., P_i, ..., P_m) from *FILE*. Its absence is substituted by *stdin*. The end of a point set P_i, which equals the line #eos, starts the algorithm to emit the requested measure(s) applied on *P_i*.

**\--o** *FILE*, **\--o**=*FILE*
:   Redirects the measures to *FILE*, opened in overwrite mode (not appending mode). Without *FILE*, results are forwarded to *stdout*. Errors encountered during the program's execution are streamed into *stderr*, and not into either *stdout* or *FILE*.

**\--disp**
:   Computes dispersion of *P_i*.

**\--ndisp**
:   Computes dispersion of *P_i*, multiplied by the cardinality *n* of *P_i*.

**\--count-boxes**
:   Counts all empty boxes of *P_i*, considered as RRs in the publication.

**\--silent**
:   Suppress comments in the output stream, yielding only the computed value. The latter could be the point set or its cardinality.

# RETURN FORMAT

## Dispersion, number of boxes

A point set of cardinality *m* with each axis representing a computed measure.

point set | disp | n*disp | #RR boxes
--- | --- | --- | ---
P_0 | . | . | .
P_1 | . | . | .
... | . | . | .
P_m | . | . | .

Notice that the first column is not returned.

# LIMITATION

The algorithm requires a 3-dimensional point sets sequence.