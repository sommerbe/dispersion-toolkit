% DISPOPTGS(1) 1.0.0 | Dispersion Toolkit Manuals
% Benjamin Sommer
% November 30, 2020

# NAME

dispoptgs - a gradient ascent to reduce dispersion based on grow&shrink strategy

# SYNOPSIS

**dispoptgs** [**\--i** *FILE*] [**\--o** *FILE*] [**\--iteration-limit**=*INTEGER*] [**\--tau**=*BINARY64*] [**\--stepsize**=*BINARY64*] [**\--delimiter**=*CHARACTER*] [**\--no-pointset**] [**\--compute-iterations**] [**\--pointset-sequence**] [**\--silent**]

# DESCRIPTION

Seeks to reduce dispersion of a point set by using gradient ascent on local dispersion. *Local dispersion* at a point p is the area of the greatest empty box out of all empty boxes of the point set which contain p.

The axis aligned gradient of local dispersion at each point p of the point set equals the vector with which p is moved iteratively during the ascent. This vector is scaled with a step size.

The ascent terminates as soon as the greatest gradient magnitude falls below a threshold, or if upon reaching an iteration limit.

# OPTIONS

**\--i** *FILE*, **\--i**=*FILE*
:   Retrieves a point set sequence (P_0, P_1, ..., *P_i*, ..., P_m) from *FILE*. Its absence is substituted by *stdin*. The end of a point set *P_i*, which equals the line #eos, starts the algorithm to emit the requested measure(s) applied on *P_i*.

**\--o** *FILE*, **\--o**=*FILE*
:   Redirects point set(s) to *FILE*, opened in overwrite mode (not appending mode). Without *FILE*, results are forwarded to *stdout*. Errors encountered during the program's execution are streamed into *stderr*, and not into either *stdout* or *FILE*.

**\--iteration-limit**=*INTEGER*, **\--c**=*INTEGER*
:   The ascent terminates as soon as reaching *INTEGER* iterations. Boundary condition: *INTEGER* > 0.

**\--tau**=*BINARY64*, **\--t**=*BINARY64*
:   The ascent terminates as soon as the greatest gradient magnitude falls below *BINARY64*. Boundary condition: 0 < *BINARY64* << 1, and ideally close to machine zero.

**\--stepsize**=*BINARY64*, **\--dt**=*BINARY64*
:   In each iteration, the gradients are scaled with *BINARY64*, with which elements of the point set are moved. Recommended condition: 0 < *BINARY64* < 1.

**\--delimiter**=*CHARACTER*
:   A point set's coordinates are separated by *CHARACTER*, while each point resides on a distinct line. Recommended condition: *CHARACTER* = \' \', *CHARACTER* = \'\\t\'.

**\--no-pointset**
:   No point set is emitted. This option is only useful (and valid) in combination with **\--compute-iterations**.

**\--compute-iterations**
:   Feed the number of actual iterations needed to the chosen output. This number is returned before the optimised point set.

**\--pointset-sequence**
:   Feed each intermediate point set, of each iteration, to the output. This option is required to visualise the modifications made, for instance using **pssy.py**.
   
**\--silent**
:   Suppress comments in the output stream, yielding only the computed value. The latter could be the point set or its cardinality.

# LIMITATION

The algorithm requires a two-dimensional point set sequence.