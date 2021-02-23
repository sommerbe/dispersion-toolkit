% PSS.PY(1) 1.2.0 | Dispersion Toolkit Manuals
% Benjamin Sommer
% November 30, 2020

# NAME

pss.py - visualise a point set sequence

# SYNOPSIS

**pss.py** [**\--i** *FILE*] [**\--delay**=*BINARY64*] [**\--image-path**=*PATTERN*] [**\--image-ppi**=*BINARY64*] [**\--silent**]

# DESCRIPTION

Visualises points of a point set sequence, how they change w.r.t. the first point set of this sequence.

# OPTIONS

**\--i** *FILE*, **\--i**=*FILE*
:   Retrieves a point set sequence from *FILE*. Its absence is substituted by *stdin*. The end of a point set, which equals the line #eos, updates the figure.

**\--delay**=*BINARY64*
:   The animation speed is slowed down by introducing a delay of *BINARY64* seconds between frame updates. A frame update corresponds to the end of point set within the sequence.

**\--image-path**=*PATTERN*
:   Each frame may be stored to a file as a realisation of formatting *PATTERN*. This *PATTERN* is required to contain *{i}* to be substituted by the current point set index within the sequence. Examples:

    * seq-{i}.png
    * seq-{i}.pgf

**\--image-ppi**=*BINARY64*
:   The resolution of the resulting images, pixel per inch (PPI) unit. Default value: 300.

**\--silent**
:   Suppress comments in the output stream, yielding only the computed value. The latter could be the point set or its cardinality.

# DYNAMIC DEPENDENCIES

* python3
* numpy
* matplotlib

# LIMITATION

The algorithm requires a two-dimensional point set sequence.