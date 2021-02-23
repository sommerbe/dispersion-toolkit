% UTK-ADAPTER(1) 1.2.0 | Dispersion Toolkit Manuals
% Benjamin Sommer
% December 1, 2020

# NAME

read-utk - transform a point set sequence from UTK framework

# SYNOPSIS

**read-utk** [**\--i** *FILE*] [**\--o** *FILE*] [**\--delimiter**=*CHARACTER*] [**\--silent**]

# DESCRIPTION

Converts a point set sequence in the .dat format as produced by the UTK framework, https://utk-team.github.io/utk, into the format of a point set sequence defined by the present dispersion toolkit.

# OPTIONS

**\--i** *FILE*, **\--i**=*FILE*
:   Retrieves a point set sequence (P_0, P_1, ..., *P_i*, ..., P_m) from *FILE*. Its absence is substituted by *stdin*. The end of a point set *P_i* starts the conversion to emit the equivalent stream of *P_i*.

**\--o** *FILE*, **\--o**=*FILE*
:   Redirects the point set sequence to *FILE*, opened in overwrite mode (not appending mode). Without *FILE*, results are forwarded to *stdout*. Errors encountered during the program's execution are streamed into *stderr*, and not into either *stdout* or *FILE*.

**\--delimiter**=*CHARACTER*
:   A point set's coordinates are separated by *CHARACTER*, while each point resides on a distinct line. Recommended condition: *CHARACTER* = \' \', *CHARACTER* = \'\\t\'.

**\--silent**
:   Suppress comments in the output stream, yielding only the computed value. The latter could be the point set or its cardinality.

# LIMITATION

The algorithm requires a point set sequence in .dat as produced by the UTK framework. File formats .edat, .bin and .ebin remain unsupported.