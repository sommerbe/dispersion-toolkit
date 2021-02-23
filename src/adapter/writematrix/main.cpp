#include "../../io/argparse.hpp"
#include "../../io/ipointset.hpp"
#include "../../io/opointset.hpp"
#include "../../io/ostream.hpp"
#include "../../math/pointset.hpp"
#include <iomanip>
#include <iostream>
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
  i8            delimiter;
  u1            transpose;
  u1            del_use_ipts;
};

struct problem_param
{
  pointset       pts;
  program_param* rt;
};

void write_matrix(problem_param& p)
{
  prec* c;

  ensure_precision(p.rt->os, *c);

  for (u64 i = 0; i < p.pts.size(); ++i) {
    c = p.pts.at(i, 0);
    for (u64 j = 0; j < p.pts.dimensions; ++j) {
      if (j > 0) {
        *p.rt->os << p.rt->delimiter;
      }
      *p.rt->os << c[j];
    }
    *p.rt->os << std::endl;
  }
};

void transpose(problem_param& p)
{
  pointset pts;
  prec*    c;

  pts.allocate(p.pts.dimensions, p.pts.points);
  pts.reset_inf_bound();

  for (u64 i = 0; i < p.pts.size(); ++i) {
    c = p.pts.at(i, 0);
    for (u64 j = 0; j < p.pts.dimensions; ++j) {
      *pts.at(j, i) = c[j];
    }
  }

  p.pts = pts;
};

u1 parse_progargs(i32 argc, const i8** argv, program_param& rt)
{
  assert(argc >= 0);

  std::vector<std::string> arg;

  // generate a canonical representation of argument parameters
  // i.e. one without abbreviations (compressed syntax)
  argparse::canonical(argc, argv, arg);

  for (u64 i = 0; i < arg.size(); ++i) {
    const std::string& s = arg[i];

    if (s == "--delimiter") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing delimiter value. Consider using -h or --help.");
      rt.delimiter    = arg[++i][0];
      rt.del_use_ipts = false;

    } else if (s == "--transpose" || s == "-t") {
      rt.transpose = true;

    } else if (s == "--silent" || s == "-s") {
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
      std::cout << "NAME: writes as a matrix of points" << std::endl;
      std::cout << "SYNOPSIS: [-ts]  [--i  FILE]  [--o  FILE]  [--transpose]  "
                   "[--delim‐iter=CHARACTER] [--silent]"
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
  rt.silent       = false;
  rt.delimiter    = ' ';
  rt.del_use_ipts = true;
  rt.transpose    = false;
  rt.input        = "-";
  rt.output       = "-";
  problem.rt      = &rt;
  r               = EXIT_SUCCESS;

  // parse arguments
  if (!dptk::parse_progargs(argc, argv, rt)) {
    return EXIT_FAILURE;
  }

  // initialize io streams
  dptk::istream_init(rt.input, rt.is);
  dptk::ostream_init(rt.output, rt.os);

  assert(rt.is != nullptr);
  assert(rt.os != nullptr);

  // iterate through pointset sequence
  while (!rt.is->eof() && r == EXIT_SUCCESS) {
    // clear pointset
    problem.pts.clear();

    // retrieve point set
    dptk::read_pointset(*rt.is, problem.pts, &ipts_inf);

    // skip empty points
    if (problem.pts.coords.empty()) {
      continue;
    }

    dptk::forward_delimiter(rt.del_use_ipts, ipts_inf, rt.delimiter);

    // skip empty point sets
    if (problem.pts.coords.empty()) {
      continue;
    }

    // option: transpose
    if (rt.transpose) {
      dptk::transpose(problem);
    }

    // show result
    dptk::write_matrix(problem);
  }

  // clean up (heap allocations)
  dptk::istream_close(rt.is);
  dptk::ostream_close(rt.os);

  return r;
}