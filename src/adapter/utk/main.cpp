#include "../../io/argparse.hpp"
#include "../../io/ipointset.hpp"
#include "../../io/opointset.hpp"
#include "../../io/ostream.hpp"
#include "../../math/pointset.hpp"
#include "manpage.hpp"
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
};

struct problem_param
{
  pointset       pts;
  program_param* rt;
};

void read_utk_dat(std::istream& in, pointset& out)
{
  std::string ln;
  i8*         b;
  i8*         e;
  i8*         y;
  i8*         z;
  u64         d;
  b64         c;

  while (std::getline(in, ln) && !in.fail()) {
    if (ln.empty()) {
      continue;
    }
    if (ln == "#") {
      break;
    }
    if (ln[0] == '#') {
      continue;
    }
    z = (i8*)ln.data();
    y = &ln[ln.size()];
    b = z;
    d = 0;
    while (b != y) {
      if (std::isspace(*b) != 0) {
        ++b;
        continue;
      }
      c = std::strtod(b, &e);
      out.coords.push_back(c);
      ++d;
      b = e;
    }
    // add point
    if (d > 0) {
      ++out.points;
    }

    // check: dimension match
    if (out.dimensions == 0) {
      out.dimensions = d;
    } else if (out.dimensions != d) {
      std::string m = "pointset contains points of variable dimensions violating the "
                      "regular pointset condition. ";
      m += "expected dimension = " + std::to_string(out.dimensions);
      m += "read dimension = " + std::to_string(d);
      m += "line = " + ln;
      m += ".";
      throw std::runtime_error(m);
    }
  }
};

i32 return_results(const program_param& rt, const problem_param& problem)
{
  if (!rt.silent)
    *rt.os << "# src=" << rt.input << std::endl;

  putln(rt.os, "# coordinates of points:", !rt.silent);
  putln(rt.os, "# (coord_0 coord_1 ... coord_n):", !rt.silent);
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

    if (s == "--delimiter") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing delimiter value. Consider using -h or --help.");
      rt.delimiter = arg[++i][0];

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
    } else if (s == "-h") {
      std::cout << extract_range(manpage, "NAME", "OPTIONS");
      std::cout << "Option --help expands this manual." << std::endl;
      return false;
    } else if (s == "--help") {
      std::cout << manpage;
      return false;
    }
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
  rt.silent    = false;
  rt.delimiter = ' ';
  rt.input     = "-";
  rt.output    = "-";
  problem.rt   = &rt;
  r            = EXIT_SUCCESS;

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
    dptk::read_utk_dat(*rt.is, problem.pts);

    // skip empty point sets
    if (problem.pts.coords.empty()) {
      continue;
    }

    // show result
    r = dptk::return_results(rt, problem);
  }

  // clean up (heap allocations)
  dptk::istream_close(rt.is);
  dptk::ostream_close(rt.os);

  return r;
}