% MINDISP-GS(1) 1.2.0 | Dispersion Toolkit Manuals
% Benjamin Sommer
% December 14, 2020

# NAME

mindisp-gs - a gradient ascent to reduce dispersion based on grow&shrink strategy

# SYNOPSIS

**mindisp-gs** [**\--i** *FILE*] [**\--o** *FILE*] [**\--iteration-limit**=*INTEGER*] [**\--tau**=*BINARY64*] [**\--stepsize**=*BINARY64*] [**\--delimiter**=*CHARACTER*] [**\--no-pointset**] [**\--compute-iterations**] [**\--pointset-sequence**] [**\--silent**]

# DESCRIPTION

Seeks to reduce dispersion of a point set *P* with cardinality *n* and dimension *d* by using gradient ascent on local dispersion. *Local dispersion* at a point p is the area of the greatest empty box out of all empty boxes of this point set which contain p.

The axis aligned gradient of local dispersion at each point p of *P* equals the vector with which p is moved iteratively during the ascent. This vector is scaled with a step size.

The ascent terminates after *t* iterations, as soon as the greatest gradient magnitude falls below the threshold *tau*, or upon reaching an iteration limit.

Computational complexity: t * (n * n * d + n * log(n) * d).

Memory complexity: 2 * n * BINARY64 + n * d * U64 = n * (d + 2) * 64bit.

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

# RETURNS

* A *d* dimensional point set with cardinality *n* resulting from the last iteration. 
* A sequence of *d* dimensional point sets each with cardinality *n* resulting from each iteration in the presence of option **\--pointset-sequence**.