#include "../../index/pointset_dsorted.hpp"
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
  std::ostream* os;
  std::istream* is;
  u1            silent;
  u1            compute_disp;
  u1            compute_ndisp;
  u1            compute_boxcount;
};

struct problem_param
{
  pointset               pts;
  u64                    pts_end_idx;
  pointset_dsorted_index idx;
  prec                   domain_bound[2];
  u64                    box_count;
  program_param*         rt;
};

struct grow_shrink_param
{
  u64            base_idx;
  i64            head_idx;
  i64            end_idx;
  prec           bound[2];
  u32            axis;
  i32            dir;
  problem_param* problem;
};

struct dispersion_param
{
  prec disp;
};

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
  u32   o;

  b  = p.problem->idx.at(p.base_idx, p.axis);
  sb = &p.problem->pts[b];
  o  = 1 - p.axis;

  while (p.head_idx != p.end_idx) {
    // lookup index of head
    k  = p.problem->idx.at(p.head_idx, p.axis);
    sk = &p.problem->pts[k];

    // prepare for next iteration
    p.head_idx += p.dir;

    // skip out-of-bounds (along axis=o)
    if (sk[o] > p.bound[1] || sk[o] < p.bound[0])
      continue;

    // grow the bound (along axis=p.axis)
    z     = (p.bound[1] - p.bound[0]) * std::abs(sk[p.axis] - sb[p.axis]);
    ldisp = std::max(z, ldisp);

    // update counting
    if (p.problem->rt->compute_boxcount) {
      ++p.problem->box_count;
    }

    // shrink the bound (potentially; along axis=o)
    if (sk[o] > sb[o]) {
      p.bound[1] = sk[o];
    } else if (sk[o] < sb[o]) {
      p.bound[0] = sk[o];
    } else {
      // proceed above
      k          = p.head_idx;
      z          = p.bound[0];
      p.bound[0] = sk[o];
      ldisp      = std::max(grow_shrink(p), ldisp);

      // proceed below
      p.head_idx = k;
      p.bound[0] = z;
      p.bound[1] = sk[o];
      ldisp      = std::max(grow_shrink(p), ldisp);

      break;
    }
  }

  // include the box aligned to the edge
  z = (p.bound[1] - p.bound[0]);
  z *= std::abs(p.problem->domain_bound[(p.dir + 1) / 2] - sb[p.axis]);
  ldisp = std::max(z, ldisp);

  // update counting
  if (p.problem->rt->compute_boxcount) {
    ++p.problem->box_count;
  }

  return ldisp;
};

prec sided_local_disp(u64 pt, grow_shrink_param& p)
{
  assert(p.dir == -1 || p.dir == 1);

  p.bound[0] = p.problem->domain_bound[0];
  p.bound[1] = p.problem->domain_bound[1];

  p.base_idx = p.problem->idx.search(pt, p.axis);
  p.head_idx = p.base_idx + p.dir;

  if (p.dir > 0)
    p.end_idx = p.problem->pts.size();
  else
    p.end_idx = -1;

  return grow_shrink(p);
};

void local_disp(u64 pt, problem_param* problem, prec& disp)
{
  assert(problem != nullptr);
  assert(problem->pts.dimensions == 2);

  prec              ldisp;
  grow_shrink_param gsp;

  gsp.problem = problem;

  for (i32 axis = 0; axis < 2; ++axis) {
    for (i32 dir = 0; dir < 3; dir += 2) {
      gsp.axis = axis;
      gsp.dir  = dir - 1;
      ldisp    = sided_local_disp(pt, gsp);
      disp     = std::max(ldisp, disp);
      assert(ldisp >= 0);
    }
  }
};

void dispersion(dispersion_param& ga, problem_param* problem)
{
  assert(problem != nullptr);
  assert(problem->pts.size() > 0);
  assert(problem->domain_bound[0] == 0);
  assert(problem->pts.size() <= INT64_MAX);

  prec  percentile[3];
  prec* bd;

  bd                   = problem->domain_bound;
  problem->pts_end_idx = problem->pts.size() - 1;
  ga.disp              = 0;

  // rebuild index
  sort(problem->pts, problem->idx);

  // compute local dispersion and update global dispersion
  for (u64 j = 0; j < problem->pts.size(); ++j) {
    local_disp(j, problem, ga.disp);
  }
};

i32 return_results(const program_param&    rt,
                   const problem_param&    problem,
                   const dispersion_param& gap)
{
  if (!rt.silent)
    *rt.os << "# src=" << rt.input << std::endl;

  if (rt.compute_disp) {
    putln(rt.os, "# computed dispersion:", !rt.silent);
    putlnsci(rt.os, gap.disp, 16);
  }
  if (rt.compute_ndisp) {
    dptk::prec ndisp = problem.pts.size() * gap.disp;
    putln(rt.os, "# computed n*dispersion:", !rt.silent);
    putlnsci(rt.os, ndisp, 16);
  }
  if (rt.compute_boxcount) {
    putln(rt.os, "# counted number of all empty boxes:", !rt.silent);
    putln(rt.os, problem.box_count);
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

    if (s == "--disp") {
      rt.compute_disp = true;
    } else if (s == "--ndisp") {
      rt.compute_ndisp = true;
    } else if (s == "--count-boxes") {
      rt.compute_boxcount = true;
    } else if (s == "--silent") {
      rt.silent = true;
    } else if (s == "-i") {
      if (++i == arg.size())
        return argparse::err("invalid argument: -i misses a mandatory parameter");
      rt.input = arg[i];
    } else if (s == "-o") {
      if (++i == arg.size())
        return argparse::err("invalid argument: -o misses a mandatory parameter");
      rt.output = arg[i];
    } else if (s == "-h" || s == "--help") {
      std::cout << "# NAME #" << std::endl
                << "" << argv[0] << " - compute dispersion using grow&shrink algorithm"
                << std::endl
                << std::endl;
      std::cout << "# SYNOPSIS #" << std::endl;
      std::cout << "" << argv[0]
                << " [-i FILE] [-o FILE] [--disp] [--ndisp] [--count-boxes] [--silent]"
                << std::endl
                << std::endl;
      std::cout << "# DESCRIPTION #" << std::endl;
      std::cout
        << "Computes dispersion, n*dispersion and/or number of empty boxes found (close "
           "approximation from below), in this order, of a given point set sequence "
           "using -i FILE option. If -i "
           "FILE option is missing, standard input is assumed. The result will be "
           "written to standard output, or to the file given by -o FILE. The option "
           "--silent suppresses comments, yielding only the computed value."
        << std::endl
        << std::endl;
      std::cout << "# LIMITATION #" << std::endl;
      std::cout << "Given point set must be two-dimensional." << std::endl;
      return false;
    }
  }

  if (!(rt.compute_disp || rt.compute_ndisp || rt.compute_boxcount)) {
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
  dptk::dispersion_param gap;
  dptk::problem_param    problem;
  dptk::program_param    rt;
  dptk::i32              r;

  // default configuration
  rt.compute_boxcount     = false;
  rt.compute_disp         = false;
  rt.compute_ndisp        = false;
  rt.silent               = false;
  rt.input                = "-";
  rt.output               = "-";
  problem.rt              = &rt;
  problem.domain_bound[0] = 0;
  problem.domain_bound[1] = 1;
  problem.box_count       = 0;
  r                       = EXIT_SUCCESS;

  // parse arguments
  if (!dptk::parse_progargs(argc, argv, rt)) {
    return EXIT_FAILURE;
  }

  // initialize io streams
  dptk::istream_init(rt.input, rt.is);
  dptk::ostream_init(rt.output, rt.os);

  // iterate through pointset sequence
  while (!rt.is->eof() && r == EXIT_SUCCESS) {
    // clear pointset
    problem.pts.clear();

    // retrieve point set
    dptk::read_pointset(*rt.is, problem.pts);

    // skip empty points
    if (problem.pts.coords.empty()) {
      continue;
    }

    assert(problem.pts.dimensions == 2);
    assert(rt.is != nullptr);
    assert(rt.os != nullptr);

    // allocate
    problem.idx.allocate(problem.pts.size(), problem.pts.dimensions);

    // compute dispersion
    dptk::dispersion(gap, &problem);

    // show result
    r = dptk::return_results(rt, problem, gap);
  }

  // clean up (heap allocations)
  dptk::istream_close(rt.is);
  dptk::ostream_close(rt.os);

  return r;
}