#include "../../io/argparse.hpp"
#include "../../io/ipointset.hpp"
#include "../../io/opointset.hpp"
#include "../../io/ostream.hpp"
#include "../../math/pointset.hpp"
#include <cctype>
#include <iomanip>
#include <iostream>
#include <sstream>
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
  i8            point_delimiter;
  i8            coord_delimiter;
  i8            point_prefix;
  i8            point_suffix;
  i8            set_prefix;
  i8            set_suffix;
  i8            comment_prefix;
  u1            domain_boundary_unity;
  pointset      domain_boundary;
};

struct problem_param
{
  pointset       pts;
  program_param* rt;

  u64 input_line;
};

enum parser_state : u32
{
  outside_set = 0,
  outside_point,
  inside_point
};

u1 known_token(i8 head, const std::string& candidates)
{
  for (u64 i = 0; i < candidates.size(); ++i) {
    if (candidates[i] == head)
      return true;
  }
  return false;
}

u1 known_token_opt(i8 head, const std::string& candidates)
{
  if (candidates.empty())
    return true;
  return known_token(head, candidates);
}

void print_syntax_error(i8 head, problem_param& p)
{
  std::string ln;
  std::getline(*p.rt->is, ln);

  std::cerr << "syntax error during reading matrix at character '" << head << "' on line "
            << p.input_line << ". The remain input line is: " << ln << std::endl;

  p.rt->is->setstate(std::ios_base::failbit);
}

void parse_number(std::string& number, problem_param& p)
{
  b64 c;

  c = std::strtod(number.c_str(), nullptr);
  p.pts.coords.push_back(c);
  number.clear();
}

void ensure_dimension(problem_param& p)
{
  if (p.pts.dimensions > 0)
    return;
  p.pts.dimensions = p.pts.coords.size();
}

void read_matrix(std::istream& in, problem_param& p)
{
  assert(p.rt->point_delimiter != '\0');
  assert(p.rt->coord_delimiter != '\0');
  assert(p.rt->point_suffix != '\0' || p.rt->coord_delimiter != p.rt->point_delimiter);
  assert(p.rt->coord_delimiter != p.rt->point_suffix);

  std::string    ln;
  i8             head;
  u64            d;
  parser_state   state;
  u1             get_next = true;
  std::string    number;
  program_param* rt = p.rt;

  p.pts.dimensions = 0;
  state            = outside_set;
  p.input_line     = 0;

  while (in.good()) {
    if (get_next) {
      head = in.get();

      // increment line counter
      if (head == '\n') {
        ++p.input_line;
      }
    }
    get_next = true;

    // skip comment lines
    if (head == rt->comment_prefix) {
      std::getline(in, ln);
      ++p.input_line;
      continue;
    }

    // parse point
    if (inside_point == state) {
      // check: end of point
      if (head == rt->point_suffix
          || (head == rt->point_delimiter && rt->point_suffix == '\0')) {
        state    = outside_point;
        get_next = head == rt->point_suffix;

        parse_number(number, p);
        ensure_dimension(p);

        ++p.pts.points;
        ++d;

        // assert(d == p.pts.dimensions);
        if (d != p.pts.dimensions) {
          return print_syntax_error(head, p);
        }

        continue;
      }

      // check: end of coordinate
      if (head == rt->coord_delimiter) {
        parse_number(number, p);
        ++d;
        continue;
      }

      number.push_back(head);
    }

    if (outside_point == state) {
      if (head == rt->point_prefix
          || (head == rt->point_delimiter && rt->point_prefix == '\0')) {
        state = inside_point;
        d     = 0;
        continue;
      } else if (head == rt->set_suffix) {
        state = outside_set;
      } else if (rt->point_prefix == '\0' && p.pts.coords.empty()) {
        state    = inside_point;
        get_next = false;
        d        = 0;
        continue;
      }
    }

    // start point set
    if (outside_set == state) {
      if (head == rt->set_prefix || rt->set_prefix == '\0') {
        state    = outside_point;
        get_next = rt->set_prefix != '\0';
        continue;
      } else if (head == rt->set_suffix) {
        break;
      }
    }
  }
};

i32 return_results(const program_param& rt, const problem_param& problem)
{
  // check: given domain boundary and actual read dimension
  if (problem.pts.dimensions * 2 != problem.pts.domain_bound.size()) {
    std::cerr << "error: specified domain boundary mismatches the dimension of the point "
                 "set read from input. "
              << "Point set dimension = " << problem.pts.dimensions
              << ". Domain boundaries = " << problem.pts.domain_bound.size()
              << " (should be " << (problem.pts.dimensions * 2) << ")." << std::endl;
    return EXIT_FAILURE;
  }

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
      if (!argparse::retrieve(arg, i, rt.delimiter))
        return argparse::err("missing delimiter value. Consider using -h or --help.");

    } else if (s == "--domain-boundary") {
      if (!argparse::retrieve(arg, i, rt.domain_boundary.coords))
        return argparse::err("missing subdomain value. Consider using -h or --help.");
      if (rt.domain_boundary.coords.size() % 2 != 0)
        return argparse::err(
          "invalid domain-boundary value. The size of this list is not a multiple of 2.");

    } else if (s == "--domain-boundary-unity") {
      rt.domain_boundary_unity = true;

    } else if (s == "--mathematica") {
      rt.point_delimiter = ',';
      rt.coord_delimiter = ',';
      rt.point_prefix    = '{';
      rt.point_suffix    = '}';
      rt.set_prefix      = '{';
      rt.set_suffix      = '}';

    } else if (s == "--csv") {
      rt.point_delimiter = '\n';
      rt.coord_delimiter = ' ';
      rt.point_prefix    = '\0';
      rt.point_suffix    = '\0';
      rt.set_prefix      = '\0';
      rt.set_suffix      = '\0';

    } else if (s == "--point-delimiter") {
      if (!argparse::retrieve(arg, i, rt.point_delimiter))
        return argparse::err(
          "missing point-delimiter value. Consider using -h or --help.");

    } else if (s == "--coord-delimiter") {
      if (!argparse::retrieve(arg, i, rt.coord_delimiter))
        return argparse::err(
          "missing coord-delimiter value. Consider using -h or --help.");

    } else if (s == "--point-prefix") {
      if (!argparse::retrieve(arg, i, rt.point_prefix))
        return argparse::err("missing point-prefix value. Consider using -h or --help.");

    } else if (s == "--point-suffix") {
      if (!argparse::retrieve(arg, i, rt.point_suffix))
        return argparse::err("missing point-suffix value. Consider using -h or --help.");

    } else if (s == "--set-prefix") {
      if (!argparse::retrieve(arg, i, rt.set_prefix))
        return argparse::err("missing set-prefix value. Consider using -h or --help.");

    } else if (s == "--set-suffix") {
      if (!argparse::retrieve(arg, i, rt.set_suffix))
        return argparse::err("missing set-suffix value. Consider using -h or --help.");

    } else if (s == "--silent") {
      rt.silent = true;

    } else if (s == "--i") {
      if (!argparse::retrieve(arg, i, rt.input))
        return argparse::err("invalid argument: -i misses a mandatory parameter");

    } else if (s == "--o") {
      if (!argparse::retrieve(arg, i, rt.output))
        return argparse::err("invalid argument: -o misses a mandatory parameter");

    } else if (s == "-h" || s == "--help") {
      std::cout << "NAME: reads a matrix of points" << std::endl;
      std::cout << "SYNOPSIS: [--i FILE] [--o FILE] --domain-boundary LIST_OF_NUMBERS | "
                   "--domain-boundary-unity "
                   "[--mathematica] [--csv] [--point-delimiter=CHARACTER] "
                   "[--coord-delimiter=CHARACTER] "
                   "[--point-prefix=CHARACTER] [--point-suffix=CHARACTER] "
                   "[--set-prefix=CHARACTER] [--set-suffix=CHARACTER] "
                   "[--delimiter=CHARACTER] [--silent]"
                << std::endl;
      return false;
    }
  }

  if (!(rt.domain_boundary_unity || !rt.domain_boundary.coords.empty())) {
    return argparse::err(
      "Missing command option. Either specify --domain-boundary LIST_OF_NUMBERS or "
      "--domain-boundary-unity. Consider using -h or --help.");
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
  rt.silent                = false;
  rt.delimiter             = ' ';
  rt.input                 = "-";
  rt.output                = "-";
  rt.point_delimiter       = '\n';
  rt.coord_delimiter       = ' ';
  rt.point_prefix          = '\0';
  rt.point_suffix          = '\0';
  rt.set_prefix            = '\0';
  rt.set_suffix            = '\0';
  rt.comment_prefix        = '#';
  rt.domain_boundary_unity = false;
  problem.rt               = &rt;
  r                        = EXIT_SUCCESS;

  // parse arguments
  if (!dptk::parse_progargs(argc, argv, rt)) {
    return EXIT_FAILURE;
  }

  // initialize io streams
  dptk::istream_init(rt.input, rt.is);
  dptk::ostream_init(rt.output, rt.os);

  assert(rt.is != nullptr);
  assert(rt.os != nullptr);

  dptk::putparam(rt.os, "source", rt.input, !rt.silent);

  // iterate through pointset sequence
  while (rt.is->good() && r == EXIT_SUCCESS) {
    // clear pointset
    problem.pts.clear();

    // retrieve point set
    dptk::read_matrix(*rt.is, problem);

    // skip empty point sets
    if (problem.pts.coords.empty()) {
      continue;
    }

    // apply specified domain boundary
    if (rt.domain_boundary_unity) {
      problem.pts.reset_unity_bound();
    } else {
      problem.pts.domain_bound = rt.domain_boundary.coords;
    }

    // deduct argument
    problem.pts.arguments.push_back(problem.pts.size());

    // show result
    r = dptk::return_results(rt, problem);
  }

  // clean up (heap allocations)
  dptk::istream_close(rt.is);
  dptk::ostream_close(rt.os);

  return r;
}