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
  std::string   output;
  std::ostream* os;
  u1            silent;
  u1            compute_cardinality;
  u1            compute_pointset;
  u32           fibonacci_index;
  i8            delimiter;
};

struct problem_param
{
  pointset       pts;
  u64            fibonacci_number;
  program_param* rt;
};

void fibonacci(u64 kmax, u64& fkm2, u64& fkm1, u64& fkm0)
{
  fkm2 = 1;
  fkm1 = 1;
  fkm0 = fkm2 + fkm1;

  for (u64 k = 4; k <= kmax; ++k) {
    fkm2 = fkm1;
    fkm1 = fkm0;
    fkm0 = fkm2 + fkm1;
  }
};

void init_domain(pointset& pts)
{
  pts.domain_bound = { 0., 0., 1., 1. };
}

void fibonacci_lattice(problem_param* problem)
{
  assert(problem != nullptr);
  assert(problem->rt->fibonacci_index > 2);

  u64   Fm;
  u64   Fm1;
  u64   Fm2;
  prec  Fm_inv;
  prec* pt;
  prec  integral_part;

  // store problem domain
  init_domain(problem->pts);

  // compute last 3 Fibonacci numbers
  fibonacci(problem->rt->fibonacci_index, Fm2, Fm1, Fm);

  problem->fibonacci_number = Fm;

  // no need to compute lattice as requested
  if (!(problem->rt->compute_pointset || problem->rt->compute_cardinality)) {
    return;
  }

  Fm_inv = prec(1) / prec(Fm);

  // allocate number of points in two-dimensions
  problem->pts.allocate(Fm, 2);

  for (u64 k = 0; k < Fm; ++k) {
    pt    = problem->pts.at(k, 0);
    pt[0] = k * Fm_inv;
    pt[1] = std::modf((k * Fm2) * Fm_inv, &integral_part);
  }
};

i32 return_results(const program_param& rt, const problem_param& problem)
{
  if (rt.compute_cardinality) {
    putln(rt.os, "# computed set's cardinality (== Fibonacci number):", !rt.silent);
    putln(rt.os, problem.fibonacci_number);
  }
  if (rt.compute_pointset) {
    putln(rt.os, "# coordinates of points:", !rt.silent);
    putln(rt.os, "# (coord_0 coord_1):", !rt.silent);
    write_pointset(rt.os, problem.pts, rt.delimiter);
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

    if (s == "--fibonacci-index" || s == "--i") {
      if (!argparse::argval(arg, i))
        return argparse::err(
          "missing fibnacci-index value. Consider using -h or --help.");
      rt.fibonacci_index = std::strtol(arg[++i].c_str(), nullptr, 10);
      if (rt.fibonacci_index < 3)
        return argparse::err(
          "invalid argument: value to --fibonacci-index needs to be >= 3.");

    } else if (s == "--delimiter") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing delimiter value. Consider using -h or --help.");
      rt.delimiter = arg[++i][0];

    } else if (s == "--cardinality" || s == "--compute-fibonacci-number") {
      rt.compute_cardinality = true;
    } else if (s == "--no-pointset") {
      rt.compute_pointset = false;
    } else if (s == "--silent") {
      rt.silent = true;
    } else if (s == "--o") {
      if (++i == arg.size()) {
        std::cerr << "invalid argument: -o misses a mandatory parameter" << std::endl;
        return false;
      }
      rt.output = arg[i];
    } else if (s == "-h") {      
      std::cout << extract_range(manpage, "NAME", "MANDATORY");
      std::cout << "Option --help expands this manual." << std::endl;
      return false;
    } else if (s == "--help") {
      std::cout << manpage;
      return false;
    }
  }

  if (!(rt.compute_cardinality || rt.compute_pointset)) {
    std::cerr
      << "fatal error: unsupported output option. Consider program option -h or --help."
      << std::endl;
    return false;
  }

  if (rt.fibonacci_index < 3) {
    return argparse::err(
          "missing argument: value to --fibonacci-index needs to be >= 3.");
  }

  return true;
}

} // namespace dptk

dptk::i32 main(dptk::i32 argc, const dptk::i8** argv)
{
  dptk::problem_param problem;
  dptk::program_param rt;

  // default configuration
  rt.compute_cardinality   = false;
  rt.compute_pointset      = true;
  rt.silent                = false;
  rt.output                = "-";
  rt.delimiter             = ' ';
  rt.fibonacci_index = 0;
  problem.rt               = &rt;
  problem.fibonacci_number = 0;

  problem.pts.clear();

  // parse arguments
  if (!dptk::parse_progargs(argc, argv, rt)) {
    return EXIT_FAILURE;
  }

  // initialize io streams
  dptk::ostream_init(rt.output, rt.os);
  assert(rt.os != nullptr);

  // compute dispersion
  dptk::fibonacci_lattice(&problem);

  // show result
  dptk::i32 r = dptk::return_results(rt, problem);

  // clean up (heap allocations)
  dptk::ostream_close(rt.os);

  return r;
}
