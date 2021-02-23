% DISP-GS(1) 1.2.0 | Dispersion Toolkit Manuals
% Benjamin Sommer
% December 3, 2020

# NAME

disp-gs - compute dispersion using a grow&shrink algorithm

# SYNOPSIS

**disp-gs** [**\--i** *FILE*] [**\--o** *FILE*] [**\--disp**] [**\--ndisp**] [**\--count-boxes**] [**\--silent**]

# DESCRIPTION

Computes dispersion, *n* * dispersion and/or the number of empty boxes of a given point set with cardinality *n* and dimension *d*.

Computational complexity: n * n * d + n * log(n) * d.

Memory complexity: n * BINARY64 + n * d * U64 = n * (d + 1) * 64bit.

The measures are written to *standard output*, or to the file given by **\--o** *FILE*.

# OPTIONS

**\--i** *FILE*, **\--i**=*FILE*
:   Retrieves a point set sequence (P_0, P_1, ..., *P_i*, ..., *P_m*) from *FILE*. Its absence is substituted by *stdin*. The end of a point set *P_i*, which equals the line #eos, starts the algorithm to emit the requested measure(s) applied on *P_i*.

**\--o** *FILE*, **\--o**=*FILE*
:   Redirects the measures to *FILE*, opened in overwrite mode (not appending mode). Without *FILE*, results are forwarded to *stdout*. Errors encountered during the program's execution are streamed into *stderr*, and not into either *stdout* or *FILE*.

**\--disp**
:   Computes dispersion of *P_i*.

**\--ndisp**
:   Computes dispersion of *P_i*, multiplied by the cardinality *n* of *P_i*.

**\--count-boxes**
:   Counts all empty boxes of *P_i*, including interiour and exteriour boxes.

**\--silent**
:   Suppress comments in the output stream, yielding only the computed value. The latter could be the point set or its cardinality.

# RETURN FORMAT

A point set of cardinality *m* with each axis representing a computed measure.

point set | disp | n*disp | #boxes
--- | --- | --- | ---
P_0 | . | . | .
P_1 | . | . | .
... | . | . | .
P_m | . | . | .

Notice that the first column is not returned.