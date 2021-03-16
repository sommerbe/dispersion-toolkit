#include "../../io/argparse.hpp"
#include "../../io/ipointset.hpp"
#include "../../io/opointset.hpp"
#include "../../io/ostream.hpp"
#include "../../math/arithmetic.hpp"
#include "../../math/pointset.hpp"
#include <iomanip>
#include <iostream>
#include <random>
#include <stdlib.h>

namespace dptk {

typedef b64                    prec;
typedef regular_pointset<prec> pointset;

const prec EQUALITY_THRESHOLD = 1e-10;

enum finite_differences : u32
{
  FORWARD_DIFFERENCES = 0,
  BACKWARD_DIFFERENCES
};

struct program_param
{
  std::string   input;
  std::string   output;
  std::istream* is;
  std::ostream* os;
  u1            graph_layout;
  u1            absolute;
  u1            silent;
  u32           diff_method;
  i8            delimiter;
  u1            del_use_ipts;
  pointset      pts;
};

struct problem_param
{
  pointset       pts;
  program_param* rt;
};

void convergance(problem_param* problem)
{
  assert(problem != nullptr);
  assert(problem->rt->pts.size() > 0);
  assert(problem->rt->pts.dimensions > 0);

  // allocate new sequence
  problem->pts = problem->rt->pts;

  problem->pts.reset_inf_bound();

  if (problem->rt->graph_layout) {
    u64 z                        = problem->pts.dimensions;
    problem->pts.domain_bound[0] = problem->rt->pts.domain_bound[0];
    problem->pts.domain_bound[z] = problem->rt->pts.domain_bound[z];
  }

  u64   d;
  prec  delta;
  prec  h;
  prec* prev;
  prec* next;

  // skip first dimension for graph layouts being the argument axis
  d = problem->rt->graph_layout;

  for (; d < problem->pts.dimensions; ++d) {
    for (u64 i = 1; i < problem->pts.points; ++i) {
      prev = problem->pts.at(i - 1, 0);
      next = problem->pts.at(i, 0);

      // argument delta
      h = 1.0;
      if (problem->rt->graph_layout) {
        h = next[0] - prev[0];
      }

      // value delta
      delta = next[d] - prev[d];

      // postprocessing: make absolute
      if (problem->rt->absolute) {
        delta = math::abs(delta);
      }

      // store
      prev[d] = delta / h;
    }
  }

  // remove last point
  problem->pts.remove_last();
}

i32 return_results(const program_param& rt, const problem_param& problem)
{
  // putln(rt.os, "# coordinates of points:", !rt.silent);
  // putln(rt.os, "# (coord_0 coord_1):", !rt.silent);
  write_pointset(rt.os, problem.pts, rt.delimiter);

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

    if (s == "--ffd") {
      rt.diff_method = FORWARD_DIFFERENCES;

    } else if (s == "--graph-layout") {
      rt.graph_layout = true;
    } else if (s == "--absolute") {
      rt.absolute = true;

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
      std::cout << "NAME: swap coordinates of a given point set" << std::endl;
      std::cout << "SYNOPSIS: [--i FILE] [--o FILE] [--ffd] [--graph-layout] "
                   "[--absolute] [--silent]"
                << std::endl;
      return false;
    }
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
  rt.diff_method  = dptk::FORWARD_DIFFERENCES;
  rt.graph_layout = false;
  rt.absolute     = false;
  rt.delimiter    = ' ';
  rt.del_use_ipts = true;
  rt.silent       = false;
  rt.input        = "-";
  rt.output       = "-";
  problem.rt      = &rt;
  r               = EXIT_SUCCESS;

  rt.pts.clear();

  // parse arguments
  if (!dptk::parse_progargs(argc, argv, rt)) {
    return EXIT_FAILURE;
  }

  // initialize io streams
  dptk::istream_init(rt.input, rt.is);
  dptk::ostream_init(rt.output, rt.os);

  // show parameters
  dptk::putparam(rt.os, "finite differences", rt.diff_method, !rt.silent);
  dptk::putparam(rt.os, "absolute", rt.absolute, !rt.silent);
  dptk::putparam(rt.os, "delimiter", rt.delimiter, !rt.silent);
  dptk::putparam(rt.os, "source", rt.input, !rt.silent);

  // retrieve point set
  dptk::read_pointset(*rt.is, rt.pts, &ipts_inf);
  dptk::forward_delimiter(rt.del_use_ipts, ipts_inf, rt.delimiter);

  assert(rt.is != nullptr);
  assert(rt.os != nullptr);

  // compute convergence sequence
  dptk::convergance(&problem);

  // show result
  r = dptk::return_results(rt, problem);

  // clean up (heap allocations)
  dptk::istream_close(rt.is);
  dptk::ostream_close(rt.os);

  return r;
}