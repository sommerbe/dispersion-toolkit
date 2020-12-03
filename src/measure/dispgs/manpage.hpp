#pragma once
#include <string>

/* generated by ./script/configure-file.py using template
 * src/measure/dispgs/manpage.hpp.in with data src/measure/dispgs/manpage.1.pre. Manual
 * edits will be OVERWRITTEN UNCONDITIONALLY. */

namespace dptk {
std::string manpage
  = R"V0G0N(DISPGS(1)                 Dispersion Toolkit Manuals                 DISPGS(1)

NAME
       dispgs - compute dispersion using a grow&shrink algorithm

SYNOPSIS
       dispgs   [--i  FILE]  [--o  FILE]  [--disp]  [--ndisp]  [--count-boxes]
       [--silent]

DESCRIPTION
       Computes dispersion, n * dispersion and/or the number of empty boxes of
       a given point set with cardinality n and dimension d.

       Computational complexity: n * n * d + n * log(n) * d.

       Memory complexity: n * BINARY64 + n * d * U64 = n * (d + 1) * 64bit.

       The  measures  are  written to standard output, or to the file given by
       --o FILE.

OPTIONS
       --i FILE, --i=FILE
              Retrieves a point set sequence (P_0, P_1, ...,  P_i,  ...,  P_m)
              from  FILE.   Its absence is substituted by stdin.  The end of a
              point set P_i, which equals the line #eos, starts the  algorithm
              to emit the requested measure(s) applied on P_i.

       --o FILE, --o=FILE
              Redirects  the  measures  to FILE, opened in overwrite mode (not
              appending mode).  Without FILE, results are forwarded to stdout.
              Errors  encountered  during the program’s execution are streamed
              into stderr, and not into either stdout or FILE.

       --disp Computes dispersion of P_i.

       --ndisp
              Computes dispersion of P_i, multiplied by the cardinality  n  of
              P_i.

       --count-boxes
              Counts all empty boxes of P_i, including interiour and exteriour
              boxes.

       --silent
              Suppress comments in the output stream, yielding only  the  com‐
              puted value.  The latter could be the point set or its cardinal‐
              ity.

RETURN FORMAT
       A point set of cardinality m with each  axis  representing  a  computed
       measure.

       point set   disp   n*disp   #boxes
       ───────────────────────────────────
       P_0         .      .        .
       P_1         .      .        .
       ...         .      .        .
       P_m         .      .        .

       Notice that the first column is not returned.

AUTHORS
       Benjamin Sommer.

1.0.0                          December 2, 2020                      DISPGS(1)
)V0G0N";
}
