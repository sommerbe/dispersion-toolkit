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
  u1            transpose;
  u1            del_use_ipts;
  i8            point_delimiter;
  i8            coord_delimiter;
  i8            point_prefix;
  i8            point_suffix;
  i8            set_prefix;
  i8            set_suffix;
};

struct problem_param
{
  pointset       pts;
  program_param* rt;
};

void feed_token(std::ostream* os, i8 token)
{
  if (token == '\0')
    return;
  *os << token;
}

void write_matrix(problem_param& p)
{
  prec*          c;
  program_param* rt;

  rt = p.rt;

  ensure_precision(rt->os, *c);

  feed_token(rt->os, rt->set_prefix);
  for (u64 i = 0; i < p.pts.size(); ++i) {
    feed_token(rt->os, rt->point_prefix);
    c = p.pts.at(i, 0);
    for (u64 j = 0; j < p.pts.dimensions; ++j) {
      if (j > 0) {
        feed_token(rt->os, rt->coord_delimiter);
      }
      *rt->os << c[j];
    }
    feed_token(rt->os, rt->point_suffix);
    feed_token(rt->os, rt->point_delimiter);
  }
  feed_token(rt->os, rt->set_suffix);
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

    if (s == "--mathematica") {
      rt.point_delimiter = ',';
      rt.coord_delimiter = ',';
      rt.point_prefix    = '{';
      rt.point_suffix    = '}';
      rt.set_prefix      = '{';
      rt.set_suffix      = '}';
      rt.del_use_ipts    = false;

    } else if (s == "--csv") {
      rt.point_delimiter = '\n';
      rt.coord_delimiter = ' ';
      rt.point_prefix    = '\0';
      rt.point_suffix    = '\0';
      rt.set_prefix      = '\0';
      rt.set_suffix      = '\0';
      rt.del_use_ipts    = false;

    } else if (s == "--point-delimiter") {
      if (!argparse::argval(arg, i))
        return argparse::err(
          "missing point-delimiter value. Consider using -h or --help.");
      rt.point_delimiter = arg[++i][0];
      rt.del_use_ipts    = false;

    } else if (s == "--coord-delimiter") {
      if (!argparse::argval(arg, i))
        return argparse::err(
          "missing coord-delimiter value. Consider using -h or --help.");
      rt.coord_delimiter = arg[++i][0];
      rt.del_use_ipts    = false;

    } else if (s == "--point-prefix") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing point-prefix value. Consider using -h or --help.");
      rt.point_prefix = arg[++i][0];
      rt.del_use_ipts = false;

    } else if (s == "--point-suffix") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing point-suffix value. Consider using -h or --help.");
      rt.point_suffix = arg[++i][0];
      rt.del_use_ipts = false;

    } else if (s == "--set-prefix") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing set-prefix value. Consider using -h or --help.");
      rt.set_prefix   = arg[++i][0];
      rt.del_use_ipts = false;

    } else if (s == "--set-suffix") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing set-suffix value. Consider using -h or --help.");
      rt.set_suffix   = arg[++i][0];
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
      std::cout << "SYNOPSIS: [-ts] [--i FILE] [--o FILE] [--mathematica] [--csv] "
                   "[--point-delimiter] [--coord-delimiter] [--point-prefix] "
                   "[--point-suffix] [--set-prefix] [--set-suffix] [--transpose] "
                   "[--delimiter=CHARACTER] [--silent]"
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
  rt.silent          = false;
  rt.del_use_ipts    = true;
  rt.transpose       = false;
  rt.input           = "-";
  rt.output          = "-";
  rt.point_delimiter = '\n';
  rt.coord_delimiter = ' ';
  rt.point_prefix    = '\0';
  rt.point_suffix    = '\0';
  rt.set_prefix      = '\0';
  rt.set_suffix      = '\0';
  problem.rt         = &rt;
  r                  = EXIT_SUCCESS;

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

    dptk::forward_delimiter(rt.del_use_ipts, ipts_inf, rt.coord_delimiter);

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