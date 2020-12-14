#pragma once
#include <string>

/* generated by ./script/configure-file.py using template src/set/fibonacci/manpage.hpp.in with data src/set/fibonacci/manpage.1.pre. Manual edits will be OVERWRITTEN UNCONDITIONALLY. */

namespace dptk
{
 std::string manpage = R"V0G0N(FIBONACCILATTICE(1)       Dispersion Toolkit Manuals       FIBONACCILATTICE(1)

NAME
       fibonaccilattice - compute a Fibonacci lattice

SYNOPSIS
       fibonaccilattice  --fibonacci-index|--i=INTEGER  [--o FILE] [--compute-
       fibonacci-number|--cardinality] [--no-pointset] [--delimiter=CHARACTER]
       [--silent]

DESCRIPTION
       Computes  the  Fibonacci lattice given the Fibonacci index --fibonacci-
       index INTEGER, which needs to be >= 3.

       Unless the option --no-pointset is  given,  the  resulting  lattice  is
       written to standard output, or to the file given by --o FILE.

MANDATORY ARGUMENTS
       --fibonacci-index=INTEGER, --i=INTEGER, --i INTEGER
              The INTEGER equals the Fibonacci index m of the Fibonacci number
              F_m.  The cardinality of the resulting  point  set  equals  F_m.
              Boundary condition: m > 2.

OPTIONS
       --o FILE, --o=FILE
              Redirects the computed results to FILE, opened in overwrite mode
              (not appending mode).  Without FILE, results  are  forwarded  to
              stdout.   Errors  encountered during the program’s execution are
              streamed into stderr, and not into either stdout or FILE.

       --compute-fibonacci-number, --cardinality
              Computes the point set’s cardinality and feeds it to the  output
              stream.   In the case of the Fibonacci lattice, this cardinality
              equals the Fibonacci number F_m, where m is the Fibonacci  index
              (see above).

       --delimiter=CHARACTER
              Define  a delimiter character to separate coordinates of a point
              in the resulting output.  A usual CHARACTER  could  be  '  '  or
              '\t', for instance.

       --silent
              Suppress  comments  in the output stream, yielding only the com‐
              puted value.  The latter could be the point set or its cardinal‐
              ity.

LIMITATION
       The Fibonacci lattice will be two-dimensional, only.

AUTHORS
       Benjamin Sommer.

1.0.0                          November 29, 2020           FIBONACCILATTICE(1)
)V0G0N";
}
