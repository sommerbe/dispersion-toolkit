#include "../../io/argparse.hpp"
#include "../../io/ipointset.hpp"
#include "../../io/opointset.hpp"
#include "../../io/ostream.hpp"
#include "../../math/arithmetic.hpp"
#include "../../math/pointset.hpp"
#include "manpage.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdlib.h>

namespace dptk {

typedef b64                    prec;
typedef regular_pointset<prec> pointset;

struct program_param
{
  std::string   input;
  std::string   output;
  std::istream* is;
  std::ostream* os;
  u1            silent;
  pointset      percentiles;
  pointset      pts;
  i8            delimiter;
  u1            compute_arithmetic_mean;
  u1            compute_iqr_box;
};

struct problem_param
{
  pointset       percentiles;
  pointset       arithmetic_mean;
  program_param* rt;
};

void percentiles(problem_param* problem)
{
  assert(problem != nullptr);
  assert(problem->rt->pts.size() > 0);

  pointset       pts_axis;
  program_param* rt;

  rt = problem->rt;

  // allocate statistics: percentiles
  // [ ] axis_0 axis_1 .. axis_d
  // p_0 ...
  // p_1        ...
  // ..
  // p_i                  ...
  problem->percentiles.allocate(rt->percentiles.coords.size(), rt->pts.dimensions);

  // allocate statistics: arithmetic mean
  if (rt->compute_arithmetic_mean) {
    problem->arithmetic_mean.allocate(1, rt->pts.dimensions);
  }

  // compute percentiles for each dimension
  for (u64 d = 0; d < rt->pts.dimensions; ++d) {
    // extract points along dimension d (contiguous array needed for sorting)
    rt->pts.extract(d, pts_axis);

    // sort values in ascending order
    std::sort(pts_axis.coords.begin(), pts_axis.coords.end());

    // retrieve percentiles (using rounding to nearest integer)
    for (u64 i = 0; i < rt->percentiles.coords.size(); ++i) {
      prec p = rt->percentiles.coords[i];
      u64  j = math::round(p * (pts_axis.points - 1));

      assert(j < pts_axis.points);

      *problem->percentiles.at(i, d) = pts_axis.coords[j];
    }

    // option: compute arithmetic mean
    // - using incremental sum to increase numerical accuracy
    if (rt->compute_arithmetic_mean) {
      prec* m = problem->arithmetic_mean.at(0, d);

      *m = 0;
      for (u64 i = 0; i < pts_axis.points;) {
        u64 j = i + 1;
        *m *= prec(i) / prec(j);
        *m += pts_axis.coords[i] / prec(j);
        i = j;
      }
    }
  }
}

void compute_iqr_box(problem_param* problem)
{
  assert(problem->rt->compute_iqr_box);
  assert(problem->percentiles.points == 3);

  pointset ext;

  // allocate extended range (Q1 − 1.5 IQR, Q1, Q2=MEDIAN, Q3, Q3 + 1.5 IQR)
  ext.allocate(5, problem->percentiles.dimensions);

  for (u64 d = 0; d < problem->percentiles.dimensions; ++d) {
    // copy Q1, MEDIAN, Q3
    for (u64 i = 0; i < 3; ++i)
      *ext.at(i + 1, d) = *problem->percentiles.at(i, d);

    // compute IQR = Q3 - Q1
    prec iqr = *problem->percentiles.at(2, d) - *problem->percentiles.at(0, d);
    assert(iqr >= 0);

    // compute lower and upper whisker
    *ext.at(0, d) = *problem->percentiles.at(0, d) - 1.5 * iqr;
    *ext.at(4, d) = *problem->percentiles.at(2, d) + 1.5 * iqr;
  }

  // store results
  problem->percentiles = ext;
}

i32 return_results(const program_param& rt, const problem_param& problem)
{
  if (!rt.silent) {
    *rt.os << "# src = " << rt.input << std::endl;
  }

  if (rt.percentiles.coords.size() > 0) {
    putln(rt.os, "# percentiles:", !rt.silent);
    if (!rt.silent) {
      putln(rt.os, "# Q1 − 1.5 IQR", rt.compute_iqr_box);
      for (u64 i = 0; i < rt.percentiles.coords.size(); ++i) {
        *rt.os << "# " << rt.percentiles.coords[i] << std::endl;
      }
      putln(rt.os, "# Q3 + 1.5 IQR", rt.compute_iqr_box);
    }
    putln(rt.os, "# coordinates of points:", !rt.silent);
    putln(rt.os, "# (percentile_0 percentile_0 ... )", !rt.silent);
    putln(rt.os, "# (percentile_1 percentile_1 ... )", !rt.silent);
    putln(rt.os, "# ...", !rt.silent);
    putln(rt.os, "# (percentile_n percentile_n ... )", !rt.silent);
    write_pointset(rt.os, problem.percentiles, rt.delimiter);
  }

  if (rt.compute_arithmetic_mean) {
    putln(rt.os, "# arithmetic mean:", !rt.silent);
    putln(rt.os, "# coordinates of points:", !rt.silent);
    putln(rt.os, "# (mean_axis_0 mean_axis_1 ... mean_axis_n )", !rt.silent);
    write_pointset(rt.os, problem.arithmetic_mean, rt.delimiter);
  }

  return EXIT_SUCCESS;
}

u1 parse_progargs(i32 argc, const i8** argv, program_param& rt)
{
  assert(argc >= 0);

  std::vector<std::string> arg;

  // generate a canonical representation of argument parameters
  // i.e. one without abbreviations (compressed syntax)
  argparse::canonical(argc, argv, arg);

  for (u64 i = 0; i < arg.size(); ++i) {
    const std::string& s = arg[i];

    if (s == "--percentiles" || s == "--p") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing percentile value. Consider using -h or --help.");

      // a tuple of values is a special case of a pointset
      std::stringstream is(arg[++i]);
      rt.percentiles.clear();
      read_pointset(is, rt.percentiles);

      if (rt.percentiles.coords.empty())
        return argparse::err("invalid argument: percentile tuple may not be empty.");

      rt.compute_iqr_box = false;

    } else if (s == "--2sigma") {
      rt.percentiles.allocate(1, 3);
      rt.percentiles.coords[0] = 0.02275;
      rt.percentiles.coords[1] = 0.5;
      rt.percentiles.coords[2] = 1 - 0.02275;
      rt.compute_iqr_box       = false;

    } else if (s == "--iqr") {
      rt.percentiles.allocate(1, 3);
      rt.percentiles.coords[0] = 0.25;
      rt.percentiles.coords[1] = 0.5;
      rt.percentiles.coords[2] = 0.75;
      rt.compute_iqr_box       = false;

    } else if (s == "--iqr-box") {
      rt.percentiles.allocate(1, 3);
      rt.percentiles.coords[0] = 0.25;
      rt.percentiles.coords[1] = 0.5;
      rt.percentiles.coords[2] = 0.75;
      rt.compute_iqr_box       = true;

    } else if (s == "--mean") {
      rt.compute_arithmetic_mean = true;
    } else if (s == "--silent") {
      rt.silent = true;
    } else if (s == "--i") {
      if (++i == arg.size())
        return argparse::err("invalid argument: -i misses a mandatory parameter");
      rt.input = arg[i];
    } else if (s == "--o") {
      if (++i == arg.size())
        return argparse::err("invalid argument: -o misses a mandatory parameter");
      rt.output = arg[i];
    } else if (s == "-h" || s == "--help") {
      std::cout << manpage;
      return false;
    }
  }

  if (!(rt.compute_arithmetic_mean || rt.percentiles.size() > 0)) {
    std::cerr
      << "fatal error: unsupported output option. Consider program option -h or --help."
      << std::endl;
    return false;
  }

  return true;
}

} // namespace dptk

dptk::i32 main(dptk::i32 argc, const dptk::i8** argv)
{
  dptk::problem_param problem;
  dptk::program_param rt;
  dptk::i32           r;

  // default configuration
  // rt.axis        = 0;
  rt.compute_arithmetic_mean = false;
  rt.compute_iqr_box         = false;
  rt.delimiter               = ' ';
  rt.silent                  = false;
  rt.input                   = "-";
  rt.output                  = "-";
  problem.rt                 = &rt;
  r                          = EXIT_SUCCESS;

  rt.percentiles.clear();
  rt.pts.clear();

  // parse arguments
  if (!dptk::parse_progargs(argc, argv, rt)) {
    return EXIT_FAILURE;
  }

  // initialize io streams
  dptk::istream_init(rt.input, rt.is);
  dptk::ostream_init(rt.output, rt.os);

  // retrieve point set
  dptk::read_pointset(*rt.is, rt.pts);

  assert(rt.is != nullptr);
  assert(rt.os != nullptr);

  // estimate percentiles; and optionally the arithmetic mean
  dptk::percentiles(&problem);

  // option: compute iqr box (to feed to statistical box plots)
  if (rt.compute_iqr_box) {
    dptk::compute_iqr_box(&problem);
  }

  // show result
  r = dptk::return_results(rt, problem);

  // clean up (heap allocations)
  dptk::istream_close(rt.is);
  dptk::ostream_close(rt.os);

  return r;
}