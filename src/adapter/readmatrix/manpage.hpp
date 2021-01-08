#pragma once
#include <string>

/* generated by ./script/configure-file.py using template
 * src/adapter/readmatrix/manpage.hpp.in with data src/adapter/readmatrix/manpage.1.pre.
 * Manual edits will be OVERWRITTEN UNCONDITIONALLY. */

namespace dptk {
std::string manpage
  = R"V0G0N(READMATRIX(1)             Dispersion Toolkit Manuals             READMATRIX(1)

NAME
       readmatrix - reads a matrix of points

SYNOPSIS
       readmatrix [--i FILE] [--o FILE] [--delimiter=CHARACTER] [--silent]

DESCRIPTION
       Reads  a  matrix of points, where each column corresponds to a distinct
       dimension, and where each row represents a point.

       Every line needs to start with floating point values, or it may be  en‐
       tirely empty (nothing else).

OPTIONS
       --i FILE, --i=FILE
              Retrieves  a point set from FILE.  Its absence is substituted by
              stdin.

       --o FILE, --o=FILE
              Redirects the point set sequence to FILE,  opened  in  overwrite
              mode  (not appending mode).  Without FILE, results are forwarded
              to stdout.  Errors encountered during  the  program’s  execution
              are streamed into stderr, and not into either stdout or FILE.

       --delimiter=CHARACTER
              A point set’s coordinates are separated by CHARACTER, while each
              point resides on a distinct line.  Recommended condition:  CHAR‐
              ACTER = ' ', CHARACTER = '\t'.

       --silent
              Suppress  comments  in the output stream, yielding only the com‐
              puted value.  The latter could be the point set or its cardinal‐
              ity.

LIMITATION
       The algorithm requires a point set, not a point set sequence.

AUTHORS
       Benjamin Sommer.

1.1.0                          December 14, 2020                 READMATRIX(1)
)V0G0N";
}
