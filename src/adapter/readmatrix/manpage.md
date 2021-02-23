% READ-MATRIX(1) 1.2.0 | Dispersion Toolkit Manuals
% Benjamin Sommer
% December 14, 2020

# NAME

read-matrix - reads a matrix of points

# SYNOPSIS

**read-matrix** [**\--i** *FILE*] [**\--o** *FILE*] [**\--delimiter**=*CHARACTER*] [**\--silent**]

# DESCRIPTION

Reads a matrix of points, where each column corresponds to a distinct dimension, and where each row represents a point.

Every line needs to start with floating point values, or it may be entirely empty (nothing else).

# OPTIONS

**\--i** *FILE*, **\--i**=*FILE*
:   Retrieves a point set from *FILE*. Its absence is substituted by *stdin*.

**\--o** *FILE*, **\--o**=*FILE*
:   Redirects the point set sequence to *FILE*, opened in overwrite mode (not appending mode). Without *FILE*, results are forwarded to *stdout*. Errors encountered during the program's execution are streamed into *stderr*, and not into either *stdout* or *FILE*.

**\--delimiter**=*CHARACTER*
:   A point set's coordinates are separated by *CHARACTER*, while each point resides on a distinct line. Recommended condition: *CHARACTER* = \' \', *CHARACTER* = \'\\t\'.

**\--silent**
:   Suppress comments in the output stream, yielding only the computed value. The latter could be the point set or its cardinality.

# LIMITATION

The algorithm requires a point set, not a point set sequence.