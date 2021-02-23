% SET-CSWAP(1) 1.2.0 | Dispersion Toolkit Manuals
% Benjamin Sommer
% November 30, 2020

# NAME

set-cswap - swap coordinates of a given point set

# SYNOPSIS

**set-cswap** [**\--i** *FILE*] [**\--o** *FILE*]  [**\--count**|**\--c**=*INTEGER*] [**\--axis**=-1|*INTEGER*] [**\--percentage**|**\--p**=*BINARY64*] [**\--repeat**|**\--r**=*INTEGER*] [**\--silent**]

# DESCRIPTION

Swaps pairs of coordinates pseudo-randomly using an alternative Mersenne Twister, implemented as std::mt19937_64 in the C++11 standard. The pairs a chosen based on a uniform integer distribution.

The axis along to swap may be chosen randomly, using *-1*, or deterministically using *INTEGER* as values to the argument **\--axis**.

The stochastic algorithm enables the computation of statistics, for instance with respect to dispersion of stochastic point sets. To this end, a *point set sequence* is emitted with **\--repeat**.

The resulting lattice is written to *standard output*, or to the file given by **\--o** *FILE*.

# OPTIONS

**\--i** *FILE*, **\--i**=*FILE*
:   Retrieves a point set, not a point set sequence, from *FILE*. Its absence is substituted by *stdin*. The end of a point set, which equals the line #eos, starts the algorithm.

**\--o** *FILE*, **\--o**=*FILE*
:   Redirects the computed results to *FILE*, opened in overwrite mode (not appending mode). Without *FILE*, results are forwarded to *stdout*. Errors encountered during the program's execution are streamed into *stderr*, and not into either *stdout* or *FILE*.

**\--count**=*INTEGER*, **\--c**=*INTEGER*, --c *INTEGER*
:   The algorithm will swap *INTEGER* pairs of coordinates. By default, *INTEGER*=1.

**\--percentage**=*BINARY64*, **\--p**=*BINARY64*, --p *BINARY64*
:   The algorithm will swap *n* pairs of coordinates, with *n* = *BINARY64* * number of points.

**\--axis**=-1, **\--axis**=*INTEGER*,  **\--axis** *INTEGER*
:   The algorithm will swap coordinates along the axis *INTEGER*. In case *INTEGER*=-1, the axis is chosen with a uniformly distributed pseudo-random number based on std::mt19937_64.

**\--repeat**=*INTEGER*, **\--r**=*INTEGER*,  **\--repeat** *INTEGER*
:   A point set sequence with exactly *INTEGER* i.i.d. elements will be emitted. By default, *INTEGER*=1.

**\--delimiter**=*CHARACTER*
:   Define a delimiter character to separate coordinates of a point in the resulting output. A usual *CHARACTER* could be \' \' or \'\\t\', for instance.

**\--silent**
:   Suppress comments in the output stream, yielding only the computed value. The latter could be the point set or its cardinality.

# LIMITATION

The algorithm requires two-dimensional point sets.