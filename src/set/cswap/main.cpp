#include "../../io/argparse.hpp"
#include "../../io/ipointset.hpp"
#include "../../io/opointset.hpp"
#include "../../io/ostream.hpp"
#include "../../math/pointset.hpp"
#include "manpage.hpp"
#include <iomanip>
#include <iostream>
#include <random>
#include <stdlib.h>

namespace dptk {

typedef b64                    prec;
typedef regular_pointset<prec> pointset;

const prec EQUALITY_THRESHOLD = 1e-10;

struct program_param
{
  std::string   input;
  std::string   output;
  std::istream* is;
  std::ostream* os;
  u1            silent;
  u64           axis;
  u1            axis_random;
  u64           count;
  prec          percentage;
  i8            delimiter;
  u1            del_use_ipts;
  u64           repeat;
  pointset      pts;
};

struct problem_param
{
  pointset       pts;
  u64            pts_sequence_index;
  program_param* rt;
};

void swap_coordinates(problem_param* problem, std::mt19937_64& rengine)
{
  assert(problem != nullptr);
  assert(problem->rt->pts.size() > 0);
  assert(problem->rt->axis < problem->rt->pts.dimensions);

  // obtain local copy due to needing independent repetitions
  problem->pts = problem->rt->pts;

  // prepare uniform sampling distribution
  std::uniform_int_distribution<u64> distr(0, problem->pts.size() - 1);
  std::uniform_int_distribution<u64> distr_axis(0, 1);

  u64 i;
  u64 j;
  u64 axis;

  if (!problem->rt->silent) {
    *problem->rt->os << "# sequence index = " << problem->pts_sequence_index << std::endl;
    *problem->rt->os << "# number of swaps = " << problem->rt->count << std::endl;
  }

  // swap coordinates x times, x = problem->rt->count
  for (std::size_t z = 0; z < problem->rt->count; ++z) {
    do {
      i = distr(rengine);
      j = distr(rengine);
    } while (i == j);

    if (problem->rt->axis_random) {
      axis = distr_axis(rengine);
    } else {
      axis = problem->rt->axis;
    }

    std::swap(*problem->pts.at(i, axis), *problem->pts.at(j, axis));

    if (!problem->rt->silent) {
      *problem->rt->os << "# swapped coordinates: i=" << i << "\t, j=" << j
                       << "\t, axis=" << axis << std::endl;
    }
  }
}

i32 return_results(const program_param& rt, const problem_param& problem)
{
  putln(rt.os, "# coordinates of points:", !rt.silent);
  putln(rt.os, "# (coord_0 coord_1):", !rt.silent);
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

    if (s == "--axis") {
      // argparse::require_argval(
      //   arg, i, "missing axis value. Consider using -h or --help.");
      // rt.axis = std::strtol(arg[++i].c_str(), nullptr, 10);
      // argparse::ensure(rt.axis == -1 || rt.axis >= 0,
      //                  "invalid argument: axis = -1 or axis >= 0.");
      if (!argparse::argval(arg, i))
        return argparse::err("missing axis value. Consider using -h or --help.");
      rt.axis_random = arg[++i][0] == '-';
      if (!rt.axis_random)
        rt.axis = std::strtol(arg[i].c_str(), nullptr, 10);

    } else if (s == "--count" || s == "--c") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing count value. Consider using -h or --help.");
      rt.count = std::strtol(arg[++i].c_str(), nullptr, 10);
      if (arg[i][0] == '-' || rt.count == 0)
        return argparse::err("invalid argument: count > 0.");

    } else if (s == "--percentage" || s == "--p") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing percentage value. Consider using -h or --help.");
      rt.percentage = std::strtod(arg[++i].c_str(), nullptr);
      if (rt.percentage <= 0 || rt.percentage > 1)
        return argparse::err("invalid argument: 0 < percentage < 1");

    } else if (s == "--repeat" || s == "--r") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing repeat value. Consider using -h or --help.");
      rt.repeat = std::strtol(arg[++i].c_str(), nullptr, 10);
      if (arg[i][0] == '-' || rt.repeat == 0)
        return argparse::err("invalid argument: repeat > 0.");

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
  rt.axis         = 0;
  rt.axis_random  = false;
  rt.count        = 1;
  rt.repeat       = 1;
  rt.percentage   = 0;
  rt.delimiter    = ' ';
  rt.del_use_ipts = true;
  rt.silent       = false;
  rt.input        = "-";
  rt.output       = "-";
  problem.rt      = &rt;
  r               = EXIT_SUCCESS;

  rt.pts.clear();

  // initialize a pseudo-random number generator
  std::random_device rd;
  std::mt19937_64    rengine(rd());

  // parse arguments
  if (!dptk::parse_progargs(argc, argv, rt)) {
    return EXIT_FAILURE;
  }

  // initialize io streams
  dptk::istream_init(rt.input, rt.is);
  dptk::ostream_init(rt.output, rt.os);

  // show parameters
    dptk::putparam(rt.os, "axis", rt.axis, !rt.silent);
    dptk::putparam(rt.os, "random axis", rt.axis_random, !rt.silent);
    dptk::putparam(rt.os, "count", rt.count, !rt.silent);
    dptk::putparam(rt.os, "repeat", rt.repeat, !rt.silent);
    dptk::putparam(rt.os, "percentage", rt.percentage, !rt.silent);
    dptk::putparam(rt.os, "delimiter", rt.delimiter, !rt.silent);
    dptk::putparam(rt.os, "source", rt.input, !rt.silent);

  // retrieve point set
  dptk::read_pointset(*rt.is, rt.pts, &ipts_inf);
  dptk::forward_delimiter(rt.del_use_ipts, ipts_inf, rt.delimiter);

  assert(rt.pts.dimensions == 2);
  assert(rt.is != nullptr);
  assert(rt.os != nullptr);

  if (rt.axis >= rt.pts.dimensions) {
    dptk::argparse::err(
      "invalid argument: axis to swap exceeds dimensions of given point set");
    return EXIT_FAILURE;
  }

  // convert percentage to count number
  if (rt.percentage > 0) {
    rt.count = rt.percentage * rt.pts.size();
    rt.count = std::max(dptk::u64(1), rt.count);
  }

  // emit multiple realisations (repetitions)
  for (dptk::u64 i = 0; i < rt.repeat && r == EXIT_SUCCESS; ++i) {
    problem.pts_sequence_index = i;

    // scramble coordinates randomly
    dptk::swap_coordinates(&problem, rengine);

    // show result
    r = dptk::return_results(rt, problem);
  }

  // clean up (heap allocations)
  dptk::istream_close(rt.is);
  dptk::ostream_close(rt.os);

  return r;
}