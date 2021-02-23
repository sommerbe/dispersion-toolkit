% STAT-CONFIDENCE(1) 1.2.0 | Dispersion Toolkit Manuals
% Benjamin Sommer
% November 30, 2020

# NAME

stat-confidence - estimate confidence intervals, median and arithmetic mean

# SYNOPSIS

**stat-confidence** [**\--i** *FILE*] [**\--o** *FILE*] [**\--percentiles**=*BINARY64* *BINARY64* ...] [**\--2sigma**] [**\--iqr**] [**\--iqr-box**] [**\--mean**] [**\--silent**]

# DESCRIPTION

Estimates a confidence statistics of a given point set, usually values of a measure.

Supported statistics:

* percentiles
* inter quartile range (IQR): Q1, Q2, Q3
* extended IQR for box plots: Q1 - 1.5 IQR, Q1, Q2, Q3, Q3 + 1.5 IQR
* two-sigma rule: -2 sigma, Q2, +2 sigma
* arithmetic mean

These statistics are estimated independently for each dimension *d* of the point set, with the values of the statistics being coordinates of points, along the axis *d*. These values are ordered as shown above.

# OPTIONS

**\--i** *FILE*, **\--i**=*FILE*
:   Retrieves a point set *P* with dimension *n*, not a point set sequence, from *FILE*. Its absence is substituted by *stdin*. The end of a point set, which equals the line #eos, starts the algorithm. The coordinates are interpreted as *BINARY64*.

**\--o** *FILE*, **\--o**=*FILE*
:   Redirects the computed results to *FILE*, opened in overwrite mode (not appending mode). Without *FILE*, results are forwarded to *stdout*. Errors encountered during the program's execution are streamed into *stderr*, and not into either *stdout* or *FILE*.

**\--mean**
:   Additionally compute the arithmetic mean along each axis.

**\--silent**
:   Suppress comments in the output stream, yielding only the computed value. The latter could be the point set or its cardinality.

## Statistics

Only one of these options may be chosen for each invocation. The values of the statistic are fed to the output in order, for each axis.

**\--percentiles**=*BINARY64* *BINARY64* ..., **\--p**=*BINARY64* *BINARY64* ...
:   A list of percentiles *BINARY64* separated by whitespace. Boundary condition: 0 <= *BINARY64* <= 1. 

**\--2sigma**
:   The two-sigma rule with the median, i.e. -2 sigma, Q2, +2 sigma. Emitted percentiles: 0.02275, 0.5, 1 - 0.02275.

**\--iqr**
:   Values to compute the IQR = Q3 - Q1, along with the median Q2. Emitted percentiles: Q1 = 0.25, Q2 = 5, Q3 = 0.75.

**\--iqr-box**
:   In addition to **\--iqr**, both upper and lower whiskers used for stistical box plots, are computed. Common practice employs whiskers to determine statistical outliers. Emitted values: Q1 - 1.5 IQR, Q1, Q2, Q3, Q3 + 1.5 IQR.

# RETURN FORMAT

A point set of dimension *n* with each axis representing the tuple of *m* statistical descriptors S_i, including percentiles and derived quantities,

descriptors | axis_0 | axis_1 | ... | axis_n
--- | --- | --- | --- | ---
S_0 | . | . | ... | .
S_1 | . | . | ... | .
... | . | . | ... | .
S_m | . | . | ... | . .

Therein, axis_i relates to the statistics describing axis_i of the given point set *P*. Notice that the first column is not returned.

# LIMITATION

The algorithm requires any-dimensional point sets, while point set sequences remain unsupported.