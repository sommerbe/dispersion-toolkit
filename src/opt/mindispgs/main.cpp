#include "../../index/pointset_dsorted.hpp"
#include "../../io/argparse.hpp"
#include "../../io/ipointset.hpp"
#include "../../io/opointset.hpp"
#include "../../io/ostream.hpp"
#include "../../math/pointset.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>

namespace dptk {

typedef b64                    prec;
typedef regular_pointset<prec> pointset;

const prec EQUALITY_THRESHOLD = 1e-10;

struct program_param
{
  std::string   input;
  std::string   output;
  std::ostream* os;
  std::istream* is;
  u1            silent;
  u1            compute_pointset;
  u1            compute_sequence_size;
  u1            compute_pointset_sequence;
  prec          tau;
  prec          dt;
  u64           iteration_limit;
  i8            delimiter;
  u1            del_use_ipts;
};

struct problem_param
{
  pointset               pts;
  u64                    pts_end_idx;
  pointset_dsorted_index idx;
  std::stringstream*     os;
  const program_param*   rt;
  dptk::i32              r;
};

struct grow_shrink_param
{
  u64               base_idx;
  i64               head_idx;
  i64               end_idx;
  std::vector<prec> bound;
  u32               axis;
  i32               dir;
  problem_param*    problem;
};

struct ascent_param
{
  prec              tau;
  prec              dt;
  u64               i;
  u64               i_end;
  prec              disp;
  std::vector<prec> grad;
  prec              grad_len;
};

void msg(const ascent_param& a, const problem_param* p)
{
  assert(p->pts.size() > 0);

  if (a.i % 1000 != 0 && a.i > 10)
    return;

  if (!p->rt->silent) {
    *p->os << "# " << a.i << " grad=" << std::scientific << std::setprecision(16)
           << a.grad_len << " ndisp=" << std::scientific << std::setprecision(16)
           << (p->pts.size() * a.disp) << std::endl;
  }
}

void msg_param(const ascent_param& a, const problem_param* p)
{
  if (!p->rt->silent) {
    *p->os << "# i_end=" << std::scientific << a.i_end << std::endl
           << "# tau=" << a.tau << std::endl
           << "# dt=" << a.dt << std::endl
           << "# points=" << p->pts.size() << std::endl;
  }
}

prec grow_shrink(grow_shrink_param& p)
{
  assert(p.base_idx >= 0);
  assert(p.base_idx < p.problem->pts.size());
  assert(p.base_idx != p.end_idx);

  u64   b;
  u64   k;
  prec* sb;
  prec* sk;
  prec  ldisp = 0;
  prec  z;
  u64   d;

  b  = p.problem->idx.at(p.base_idx, p.axis);
  sb = &p.problem->pts[b];
  d  = p.problem->pts.dimensions;

  while (p.head_idx != p.end_idx) {
    // lookup index of head
    k  = p.problem->idx.at(p.head_idx, p.axis);
    sk = &p.problem->pts[k];

    // prepare for next iteration
    p.head_idx += p.dir;

    // skip out-of-bounds (along axis xs != p.axis)
    u1 outside = false;
    for (u64 xs = 0; xs < d && !outside; ++xs) {
      if (xs == p.axis)
        continue;
      outside |= sk[xs] < p.bound[xs] || sk[xs] > p.bound[d + xs];
    }
    if (outside) {
      continue;
    }

    // grow the bound (along axis=p.axis)
    z = std::abs(sk[p.axis] - sb[p.axis]);
    for (u64 xs = 0; xs < d; ++xs) {
      if (xs == p.axis)
        continue;
      z *= p.bound[d + xs] - p.bound[xs];
    }
    ldisp = std::max(z, ldisp);

    // shrink the bound (along axis xs != p.axis)
    for (u64 xs = 0; xs < d; ++xs) {
      if (xs == p.axis)
        continue;
      if (sk[xs] > sb[xs]) {
        p.bound[d + xs] = sk[xs];
      } else if (sk[xs] < sb[xs]) {
        p.bound[xs] = sk[xs];
      } else {
        p.bound[xs] = sk[xs];
      }
    }
  }

  // include the box aligned to the edge
  k = ((p.dir + 1) / 2) * d + p.axis;
  z = std::abs(p.problem->pts.domain_bound[k] - sb[p.axis]);
  for (u64 xs = 0; xs < d; ++xs) {
    if (xs == p.axis)
      continue;
    z *= p.bound[d + xs] - p.bound[xs];
  }
  ldisp = std::max(z, ldisp);

  return ldisp;
};

prec sided_local_disp(u64 pt, grow_shrink_param& p)
{
  assert(p.dir == -1 || p.dir == 1);

  p.bound = p.problem->pts.domain_bound;

  p.base_idx = p.problem->idx.search(pt, p.axis);
  p.head_idx = p.base_idx + p.dir;

  if (p.dir > 0)
    p.end_idx = p.problem->pts.size();
  else
    p.end_idx = -1;

  return grow_shrink(p);
};

// void local_disp(u64 pt, problem_param* problem, prec& disp)
// {
//   assert(problem != nullptr);

//   prec              ldisp;
//   grow_shrink_param gsp;

//   gsp.problem = problem;

//   for (i32 axis = 0; axis < problem->pts.dimensions; ++axis) {
//     for (i32 dir = 0; dir < 3; dir += 2) {
//       gsp.axis = axis;
//       gsp.dir  = dir - 1;
//       ldisp    = sided_local_disp(pt, gsp);
//       disp     = std::max(ldisp, disp);
//       assert(ldisp >= 0);
//     }
//   }
// };

void gradient_local_disp(std::size_t    pt,
                         problem_param* problem,
                         prec&          disp,
                         prec*          grad,
                         prec&          lgrad_max)
{
  assert(problem != nullptr);

  prec              ldisp;
  prec              lgrad;
  prec              lgrad_abs;
  grow_shrink_param gsp;

  disp        = 0;
  lgrad_max   = 0;
  gsp.problem = problem;

  // set gradient to zero
  for (u32 i = 0; i < problem->pts.dimensions; ++i)
    grad[i] = 0;

  // grow and shrink along each axis forward and backward
  for (u64 axis = 0; axis < problem->pts.dimensions; ++axis) {
    lgrad = 0;
    for (i32 dir = 0; dir < 3; dir += 2) {
      gsp.axis = axis;
      gsp.dir  = dir - 1;
      ldisp    = sided_local_disp(pt, gsp);
      disp     = std::max(ldisp, disp);
      lgrad += gsp.dir * ldisp;
      assert(ldisp >= 0);
    }
    lgrad_abs = std::abs(lgrad);
    if (lgrad_abs > lgrad_max) {
      for (u32 i = 0; i < problem->pts.dimensions; ++i)
        grad[i] = 0;
      grad[axis] = lgrad;
      lgrad_max  = lgrad_abs;
    }
  }
};

void gradient_ascent(ascent_param& ga, problem_param* problem)
{
  assert(problem != nullptr);
  assert(problem->pts.size() > 0);
  assert(problem->pts.size() <= INT64_MAX);
  assert(problem->pts.domain_bound.size() == 2 * problem->pts.dimensions);
  assert(ga.grad.size() == problem->pts.dimensions);

  pointset pts;
  prec*    bd;
  prec*    pt;
  u64      d;
  prec     lgrad_max;

  ga.i                 = 0;
  pts                  = problem->pts;
  bd                   = &problem->pts.domain_bound[0];
  ga.grad_len          = ga.tau;
  problem->pts_end_idx = problem->pts.size() - 1;
  d                    = problem->pts.dimensions;

  while (ga.grad_len >= ga.tau && ga.i < ga.i_end) {
    // stream intermediate pointsets to ostream
    if (problem->rt->compute_pointset_sequence) {
      write_pointset(problem->rt->os, problem->pts, problem->rt->delimiter);
    }

    // rebuild index
    sort(problem->pts, problem->idx);

    // iterate: update each points using gradient ascent
    ga.grad_len = 0;
    ga.disp     = 0;
    for (u64 j = 0; j < problem->pts.size(); ++j) {
      pt = pts.at(j, 0);

      // compute gradient
      gradient_local_disp(j, problem, ga.disp, &ga.grad[0], lgrad_max);

      for (u64 x = 0; x < d; ++x) {
        // move point by iteration scheme
        pt[x] += ga.dt * ga.grad[x];

        // check: out-of-boundary
        // pt[x] = std::fmod(pt[0] + bd[1], bd[1]);
        pt[x] = std::max(bd[x], std::min(pt[x], bd[d + x]));
      }

      // accumulate maximum gradient magnitude as termination criteria
      // - axis-aligned gradients, so max( along each axis) is justified
      ga.grad_len = std::max(ga.grad_len, lgrad_max);
    }
    // move temporary points
    problem->pts = pts;

    // finish up
    ++ga.i;
    msg(ga, problem);
  }

  if (!problem->rt->silent)
    *problem->os << "# runs=" << ga.i << std::endl;
};

i32 return_results(const program_param& rt,
                   const problem_param& problem,
                   const ascent_param&  ga)
{
  if (rt.compute_sequence_size) {
    putln(problem.os, "# sequence size of gradient ascent steps:", !rt.silent);
    putln(problem.os, ga.i);
  }

  if (rt.compute_pointset || rt.compute_pointset_sequence) {
    putln(problem.os, "# coordinates of points:", !rt.silent);
    putln(problem.os, "# (coord_0 coord_1):", !rt.silent);
    write_pointset(problem.os, problem.pts, rt.delimiter);
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

    if (s == "--iteration-limit" || s == "--c") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing iterations value. Consider using -h or --help.");
      rt.iteration_limit = std::strtol(arg[++i].c_str(), nullptr, 10);

    } else if (s == "--tau" || s == "--t") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing tau value. Consider using -h or --help.");
      rt.tau = std::strtod(arg[++i].c_str(), nullptr);

    } else if (s == "--stepsize" || s == "--dt") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing stepsize value. Consider using -h or --help.");
      rt.dt = std::strtod(arg[++i].c_str(), nullptr);

    } else if (s == "--delimiter") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing delimiter value. Consider using -h or --help.");
      rt.delimiter    = arg[++i][0];
      rt.del_use_ipts = false;

    } else if (s == "--pointset-sequence") {
      rt.compute_pointset_sequence = true;
    } else if (s == "--no-pointset") {
      rt.compute_pointset = false;
    } else if (s == "--compute-iterations") {
      rt.compute_sequence_size = true;
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
      std::cout
        << "NAME: a gradient ascent to reduce dispersion based on grow&shrink strategy"
        << std::endl;
      std::cout
        << "SYNOPSIS: [--i FILE] [--o FILE] [--iteration-limit=INTEGER] [--tau=BINARY64] "
           "[--stepsize=BINARY64] [--delimiter=CHARACTER] [--no-pointset] "
           "[--compute-iterations] [--pointset-sequence] [--silent]"
        << std::endl;

      std::cout << "# CURRENT PARAMETERS #" << std::endl;
      std::cout << "--iteration-limit=" << rt.iteration_limit << std::endl
                << "--tau=" << rt.tau << std::endl
                << "--stepsize=" << rt.dt << std::endl
                << "--delimiter='" << rt.delimiter << "'" << std::endl;
      return false;
    }
  }

  if (!(rt.compute_pointset || rt.compute_sequence_size
        || rt.compute_pointset_sequence)) {
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
  dptk::ascent_param               ga;
  dptk::program_param              rt;
  std::vector<dptk::problem_param> problems;
  dptk::i32                        r;
  dptk::ipointset_read_info        ipts_inf;

  // default configuration
  rt.tau                       = 2e-15;
  rt.iteration_limit           = 100000;
  rt.dt                        = 0.1;
  rt.silent                    = false;
  rt.compute_pointset_sequence = false;
  rt.compute_pointset          = true;
  rt.compute_sequence_size     = false;
  rt.input                     = "-";
  rt.output                    = "-";
  rt.delimiter                 = ' ';
  rt.del_use_ipts              = true;

  // parse arguments
  if (!dptk::parse_progargs(argc, argv, rt)) {
    return EXIT_FAILURE;
  }

  // prepare ascent param
  ga.tau   = rt.tau;
  ga.i_end = rt.iteration_limit;
  ga.dt    = rt.dt;

  // initialize io streams
  dptk::istream_init(rt.input, rt.is);
  dptk::ostream_init(rt.output, rt.os);

  assert(rt.is != nullptr);
  assert(rt.os != nullptr);

  // show parameters
  dptk::putparam(rt.os, "tau", rt.tau, !rt.silent);
  dptk::putparam(rt.os, "iteration limit", rt.iteration_limit, !rt.silent);
  dptk::putparam(rt.os, "stepsize, dt", rt.dt, !rt.silent);
  dptk::putparam(
    rt.os, "compute pointet sequence", rt.compute_pointset_sequence, !rt.silent);
  dptk::putparam(rt.os, "compute pointset", rt.compute_pointset, !rt.silent);
  dptk::putparam(rt.os, "compute sequence size", rt.compute_sequence_size, !rt.silent);
  dptk::putparam(rt.os, "delimiter", rt.delimiter, !rt.silent);
  dptk::putparam(rt.os, "source", rt.input, !rt.silent);

  // iterate through point set sequence
  while (!rt.is->eof()) {
    // allocate problem
    dptk::problem_param problem;

    problem.rt = &rt;

    problem.pts.clear();

    // retrieve point set
    dptk::read_pointset(*rt.is, problem.pts, &ipts_inf);
    dptk::forward_delimiter(rt.del_use_ipts, ipts_inf, rt.delimiter);

    // skip empty point sets
    if (problem.pts.coords.empty()) {
      continue;
    }

    assert(problem.pts.dimensions == 2);

    // allocate thread local output stream
    // - a shared ostream yields race conditions or
    // - performance penalty due to using mutex's and locks
    problem.os = new std::stringstream;

    // allocate index
    problem.idx.allocate(problem.pts.size(), problem.pts.dimensions);

    // store problem for parallel work
    problems.push_back(problem);
  }

// parallel iterate through point set sequence
// - due to MSVC 19, the variable i is declared as dptk::i64 (signed). Notice, however, that the OpemMP 5.1 specification allows
// "A variable of a signed or unsigned integer type." (https://www.openmp.org/spec-html/5.1/openmpsu45.html#x70-700002.11.1).
#pragma omp parallel for
  for (dptk::i64 i = 0; i < problems.size(); ++i) {
    // need local copies to ensure performance
    dptk::problem_param* p   = &problems[i];
    dptk::ascent_param   gap = ga;

    // allocate vectors
    gap.grad.resize(p->pts.dimensions);

    // compute optimisation
    dptk::gradient_ascent(gap, p);

    // show result
    p->r = dptk::return_results(rt, *p, gap);
  }

  // finalise parallel iterate
  r = EXIT_SUCCESS;
  for (dptk::u64 i = 0; i < problems.size(); ++i) {
    dptk::problem_param* p = &problems[i];

    // redirect local output buffers to program output
    *rt.os << p->os->str();

    // collect exit code (worst case)
    if (p->r == EXIT_FAILURE) {
      r = p->r;
    }

    // clean up heap allocations
    delete p->os;
  }

  // clean up (heap allocations)
  dptk::istream_close(rt.is);
  dptk::ostream_close(rt.os);

  return r;
}
