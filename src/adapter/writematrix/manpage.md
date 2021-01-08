% WRITEMATRIX(1) 1.1.0 | Dispersion Toolkit Manuals
% Benjamin Sommer
% December 18, 2020

# NAME

writematrix - writes as a matrix of points

# SYNOPSIS

**writematrix** [**\-ts**] [**\--i** *FILE*] [**\--o** *FILE*] [**\--transpose**] [**\--delimiter**=*CHARACTER*] [**\--silent**]

# DESCRIPTION

Writes as a matrix of points, where each column corresponds to a distinct dimension, and where each row represents a point.

# OPTIONS

**\--i** *FILE*, **\--i**=*FILE*
:   Retrieves a point set from *FILE*. Its absence is substituted by *stdin*.

**\--o** *FILE*, **\--o**=*FILE*
:   Redirects the point set sequence to *FILE*, opened in overwrite mode (not appending mode). Without *FILE*, results are forwarded to *stdout*. Errors encountered during the program's execution are streamed into *stderr*, and not into either *stdout* or *FILE*.

**\--transpose**, **\-t**
:   Transpose the matrix before writing to output stream.

**\--delimiter**=*CHARACTER*
:   A point set's coordinates are separated by *CHARACTER*, while each point resides on a distinct line. Recommended condition: *CHARACTER* = \' \', *CHARACTER* = \'\\t\'.

**\--silent**, **\-s**
:   Suppress comments in the output stream, yielding only the computed value. The latter could be the point set or its cardinality.

# LIMITATION

The algorithm requires a point set, not a point set sequence.