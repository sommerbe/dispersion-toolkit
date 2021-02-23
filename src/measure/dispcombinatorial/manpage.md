% DISP-COMBINATORIAL(1) 1.2.0 | Dispersion Toolkit Manuals
% Benjamin Sommer
% December 2, 2020

# NAME

disp-combinatorial - compute dispersion with a combinatorial algorithm (exhaustive search)

# SYNOPSIS

**disp-combinatorial** [**\--i** *FILE*] [**\--o** *FILE*] [**\--disp**] [**\--ndisp**] [**\--count-boxes**] [**\--boxes**] [**\--interior-boxes**] [**\--greatest-box**] [**\--silent**]

# DESCRIPTION

Computes dispersion, *n* * dispersion, the number of empty boxes, and/or a list of empty boxes of a given point set with cardinality *n*.

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
:   Counts all empty boxes of *P_i*, including interiour and exteriour boxes.

**\--boxes**
:   Feeds all empty boxes of *P_i*, including interiour and exteriour boxes, to output. With the smallest point of a box *B* being (low_0 low_1), i.e. no other point of *B* has smaller coordinates, and with the upper point of *B* being (up_0 up_1), i.e. no other point of *B* has greater coordinates, *B* is formatted as a tuple (low_0 low_1 up_0 up_1) on a distinct line of the output.

**\--interior-boxes**
:   Similar to the option **\--boxes** with the exception of *B* being restricted to integriour boxes. A Box is integriour if it does not intersect with the domain boundary of *P_i*.

**\--greatest-box**
:   Similar to the option **\--boxes** with the exception of *B* being the first greatest empty box found during the exhaustive search.

**\--silent**
:   Suppress comments in the output stream, yielding only the computed value. The latter could be the point set or its cardinality.

# RETURN FORMAT

## Greatest box

A point set of cardinality *m* with each axis representing a coordinate of a box,

point set | low_0 | low_1 | up_0 | up_1
--- | --- | --- | --- | ---
P_0 | . | . | . | .
P_1 | . | . | . | .
... | . | . | . | .
P_m | . | . | . | . ,

where low_i is the smallest coordinate along axis i, and up_i the greatest one. Notice that the first column is not returned.

## Boxes

A point set sequence, where each point set is a list of the tuple (low_0 low_1 up_0 up_1).

## Dispersion, number of boxes

A point set of cardinality *m* with each axis representing a computed measure.

point set | disp | n*disp | #boxes
--- | --- | --- | ---
P_0 | . | . | .
P_1 | . | . | .
... | . | . | .
P_m | . | . | .

Notice that the first column is not returned.

# LIMITATION

The algorithm requires a two-dimensional point sets sequence.