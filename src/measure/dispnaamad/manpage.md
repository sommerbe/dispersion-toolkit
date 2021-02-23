% DISP-NAAMAD(1) 1.2.0 | Dispersion Toolkit Manuals
% Benjamin Sommer
% December 9, 2020

# NAME

disp-naamad - compute dispersion with algorithm of Naamad et. al. 1984

# SYNOPSIS

**disp-naamad** [**\--i** *FILE*] [**\--o** *FILE*] [**\--disp**] [**\--ndisp**] [**\--count-boxes**] [**\--silent**]

# DESCRIPTION

Computes dispersion, *n* * dispersion, and/or the number of restricted rectangles (RRs) considered of a given point set with cardinality *n*, according to Naamad et. al. 1984 (https://doi.org/10.1016/0166-218X(84)90124-0)

Computational complexity: O(n^2),

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
:   Counts all empty RRs of *P_i*, considered as RRs in the publication.

**\--silent**
:   Suppress comments in the output stream, yielding only the computed value. The latter could be the point set or its cardinality.

# RETURN FORMAT

## Dispersion, number of RRs

A point set of cardinality *m* with each axis representing a computed measure.

point set | disp | n*disp | #RRs
--- | --- | --- | ---
P_0 | . | . | .
P_1 | . | . | .
... | . | . | .
P_m | . | . | .

Notice that the first column is not returned.

# LIMITATION

The algorithm requires a two-dimensional point sets sequence.
