#include "../../io/argparse.hpp"
#include "../../io/ipointset.hpp"
#include "../../io/opointset.hpp"
#include "../../io/ostream.hpp"
#include "../../math/arithmetic.hpp"
#include "../../math/pointset.hpp"
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
  std::string           input;
  std::string           output;
  std::istream*         is;
  std::ostream*         os;
  u1                    silent;
  pointset              percentiles;
  pointset              pts;
  std::vector<pointset> graphs;
  i8                    delimiter;
  u1                    del_use_ipts;
  u1                    compute_arithmetic_mean;
  u1                    compute_iqr_box;
  u1                    layout_stacked_graphs;
  prec                  stacked_graphs_deviation_limit;
};

struct problem_param
{
  pointset       percentiles;
  pointset       arithmetic_mean;
  pointset       stats;
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
  problem->percentiles.reset_inf_bound();

  // allocate statistics: arithmetic mean
  if (rt->compute_arithmetic_mean) {
    problem->arithmetic_mean.allocate(1, rt->pts.dimensions);
  } else {
    problem->arithmetic_mean.allocate(0, rt->pts.dimensions);
  }
  problem->arithmetic_mean.reset_inf_bound();

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

void merge_statistics(problem_param& p)
{
  u64 cp;
  u64 o;

  cp = p.percentiles.points + p.arithmetic_mean.points;
  o  = 0;

  p.stats.allocate(cp, p.percentiles.dimensions);
  p.stats.reset_inf_bound();
  p.stats.retrieve_points(o, p.percentiles);
  o += p.percentiles.points;
  p.stats.retrieve_points(o, p.arithmetic_mean);
}

void retrieve_stacked_graphs(program_param& rt, ipointset_read_info& ipts_inf)
{
  u32 common_args;
  u32 common_vals;

  // iterate through point set sequence
  while (!rt.is->eof()) {
    rt.graphs.emplace_back();
    read_pointset(*rt.is, rt.graphs.back(), &ipts_inf);

    if (rt.graphs.back().empty()) {
      rt.graphs.pop_back();
      continue;
    }

    if (rt.graphs.size() == 1) {
      common_args = rt.graphs.back().size();
      common_vals = rt.graphs.back().dimensions;
    }
    if (common_args != rt.graphs.back().size()) {
      std::cerr << "warning: stacked graphs vary in number of arguments (graph #"
                << rt.graphs.size() << std::endl;
    }
    if (common_vals != rt.graphs.back().dimensions) {
      std::cerr << "warning: stacked graphs vary in number of dimensions (graph #"
                << rt.graphs.size() << std::endl;
    }
  }
}

/**
 * check whether arguments of stacked graphs are similar (difference < dev)
 */
u1 check_argsim_stacked_graphs(program_param&               rt,
                               const std::vector<pointset>& graphs,
                               prec                         dev)
{
  u1   eq;
  prec v;
  prec c;

  eq = true;

  for (u64 i = 0; i < graphs[0].points; ++i) {
    v = graphs[0].at(i, 0)[0];
    for (u64 j = 1; j < graphs.size(); ++j) {
      c = graphs[j].at(i, 0)[0];
      if (math::abs(v - c) > dev) {
        eq = false;
        std::cerr << "# graph arguments not similar: " << v << " (arg=" << i
                  << ", graph=0) and " << c << " (arg=" << i << ", graph=" << j << ")"
                  << std::endl;
      }
    }
  }

  putln(rt.os,
        "# checking arguments' similarity of stacked graphs (1 := passed): ",
        eq,
        !rt.silent);

  return eq;
}

void unstack_graphs(program_param& rt)
{
  /**
   * assume: stacked graphs are equal in their shape (# points, # dimensions)
   *
   * algorithm: confidence is computed along axis d (so at(:,d)).
   * stacked graphs:
   *  - confidence is computed along graph[:].at(p, d>0), p and d fixed
   *  - the first coordinate (d=0) of each graph represents the arguments
   */
  if (rt.graphs.empty()) {
    return;
  }
  assert(rt.graphs[0].dimensions > 1);

  u64  num_graphs;
  u64  num_points;
  u64  num_dimensions;
  u64  np;
  b64* pxl;

  num_graphs     = rt.graphs[0].dimensions - 1;
  num_points     = rt.graphs.size();
  np             = rt.graphs[0].points;
  num_dimensions = np * num_graphs;

  rt.pts.clear();
  rt.pts.allocate(num_points, num_dimensions);
  rt.pts.reset_inf_bound();

  for (u64 i = 0; i < num_points; ++i) {
    pxl = rt.pts.at(i, 0);
    for (u64 j = 0; j < num_dimensions; ++j) {
      pxl[j] = *rt.graphs[i].at(j % np, 1 + j / np);
    }
  }
}

void restack_layout_stats(const std::vector<pointset>& graphs, pointset& stats)
{
  if (graphs.empty()) {
    return;
  }

  u64      num_graphs;
  u64      num_points;
  u64      num_dimensions;
  u64      num_stats;
  pointset pts;
  b64*     pxl;

  num_stats      = stats.points;
  num_graphs     = graphs[0].dimensions - 1;
  num_points     = graphs[0].points;
  num_dimensions = 1 + num_stats * num_graphs;

  pts.allocate(num_points, num_dimensions);
  pts.reset_inf_bound();

  for (u64 i = 0; i < num_points; ++i) {
    pxl = pts.at(i, 0);
    // copy arguments
    pxl[0] = graphs[0].at(i, 0)[0];

    // copy statistics
    for (u64 j = 1; j < num_dimensions; ++j) {
      u64 k  = j - 1;
      pxl[j] = *stats.at(k % num_stats, i + num_points * (k / num_stats));
    }
  }

  stats = pts;
}

i32 return_results(const program_param& rt, const problem_param& problem)
{
  // if (rt.percentiles.coords.size() > 0) {
  //   putln(rt.os, "# percentiles:", !rt.silent);
  //   if (!rt.silent) {
  //     putln(rt.os, "# Q1 − 1.5 IQR", rt.compute_iqr_box);
  //     for (u64 i = 0; i < rt.percentiles.coords.size(); ++i) {
  //       *rt.os << "# " << rt.percentiles.coords[i] << std::endl;
  //     }
  //     putln(rt.os, "# Q3 + 1.5 IQR", rt.compute_iqr_box);
  //   }
  //   putln(rt.os, "# coordinates of points:", !rt.silent);
  //   putln(rt.os, "# (percentile_0 percentile_0 ... )", !rt.silent);
  //   putln(rt.os, "# (percentile_1 percentile_1 ... )", !rt.silent);
  //   putln(rt.os, "# ...", !rt.silent);
  //   putln(rt.os, "# (percentile_n percentile_n ... )", !rt.silent);
  //   write_pointset(rt.os, problem.percentiles, rt.delimiter);
  // }

  if (problem.stats.empty()) {
    return EXIT_SUCCESS;
  }

  if (rt.layout_stacked_graphs) {
    putln(rt.os,
          "# statistical coordinates: (argument x (statistics of graph i)_i), where row "
          "x column",
          !rt.silent);
  } else {
    putln(rt.os,
          "# statistical coordinates: (descriptor_i x dimension_d), where row x column",
          !rt.silent);
  }

  write_pointset(rt.os, problem.stats, rt.delimiter);

  // if (rt.compute_arithmetic_mean) {
  //   putln(rt.os, "# arithmetic mean:", !rt.silent);
  //   putln(rt.os, "# coordinates of points:", !rt.silent);
  //   putln(rt.os, "# (mean_axis_0 mean_axis_1 ... mean_axis_n )", !rt.silent);
  //   write_pointset(rt.os, problem.arithmetic_mean, rt.delimiter);
  // }

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

    } else if (s == "--stacked-graphs") {
      rt.layout_stacked_graphs = true;

    } else if (s == "--stacked-graphs-deviation-limit") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing deviation limit of arguments of given stacked "
                             "graphs. Consider using -h or --help.");
      rt.stacked_graphs_deviation_limit = std::stod(arg[++i]);

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
      std::cout << "NAME: estimate confidence intervals, median and arithmetic mean"
                << std::endl;
      std::cout
        << "SYNOPSIS: [--i FILE] [--o FILE] [--percentiles=BINARY64  "
           "BINARY64...] [--2sigma] [--iqr] [--iqr-box] [--mean] "
           "[--stacked-graphs] [--stacked-graphs-deviation-limit=BINARY64] [--silent]"
        << std::endl;
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
  dptk::problem_param       problem;
  dptk::program_param       rt;
  dptk::i32                 r;
  dptk::ipointset_read_info ipts_inf;

  // default configuration
  // rt.axis        = 0;
  rt.compute_arithmetic_mean        = false;
  rt.compute_iqr_box                = false;
  rt.delimiter                      = ' ';
  rt.del_use_ipts                   = true;
  rt.silent                         = false;
  rt.input                          = "-";
  rt.output                         = "-";
  rt.layout_stacked_graphs          = false;
  rt.stacked_graphs_deviation_limit = -1;
  problem.rt                        = &rt;
  r                                 = EXIT_SUCCESS;

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
  if (rt.layout_stacked_graphs) {
    dptk::retrieve_stacked_graphs(rt, ipts_inf);
    dptk::unstack_graphs(rt);
    if (rt.stacked_graphs_deviation_limit >= 0) {
      dptk::check_argsim_stacked_graphs(rt, rt.graphs, rt.stacked_graphs_deviation_limit);
    }
  } else {
    dptk::read_pointset(*rt.is, rt.pts, &ipts_inf);
  }

  dptk::forward_delimiter(rt.del_use_ipts, ipts_inf, rt.delimiter);

  assert(rt.is != nullptr);
  assert(rt.os != nullptr);

  // show parameters
  dptk::putparam(rt.os, "compute IQR box", rt.compute_iqr_box, !rt.silent);
  dptk::putparam(
    rt.os, "compute arithmetic mean", rt.compute_arithmetic_mean, !rt.silent);
  dptk::putparam(rt.os, "delimiter", rt.delimiter, !rt.silent);
  dptk::putparam(rt.os, "source", rt.input, !rt.silent);

  // estimate percentiles; and optionally the arithmetic mean
  dptk::percentiles(&problem);

  // option: compute iqr box (to feed to statistical box plots)
  if (rt.compute_iqr_box) {
    dptk::compute_iqr_box(&problem);
  }

  // merge separate statistical descriptors into a single point set
  dptk::merge_statistics(problem);

  // restack back again in case of having stacked layout
  if (rt.layout_stacked_graphs) {
    dptk::restack_layout_stats(rt.graphs, problem.stats);
  }

  // show result
  r = dptk::return_results(rt, problem);

  // clean up (heap allocations)
  dptk::istream_close(rt.is);
  dptk::ostream_close(rt.os);

  return r;
}