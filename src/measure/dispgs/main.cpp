#include "../../index/pointset_dsorted.hpp"
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
  std::ostream* os;
  std::istream* is;
  u1            silent;
  i8            delimiter;
  u1            del_use_ipts;
  u1            compute_disp;
  u1            compute_ndisp;
  u1            compute_boxcount;
};

struct problem_measures
{
  prec disp;
  prec ndisp;
  u64  boxcount;
  u64  gs_diverge_count = 0;
};

struct problem_param
{
  pointset               pts;
  u64                    pts_end_idx;
  pointset_dsorted_index idx;
  problem_measures*      measures;
  program_param*         rt;
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

    // update counting
    if (p.problem->rt->compute_boxcount) {
      ++p.problem->measures->boxcount;
    }

    // shrink the bound (along axis xs != p.axis)
    for (u64 xs = 0; xs < d; ++xs) {
      if (xs == p.axis)
        continue;
      if (sk[xs] > sb[xs]) {
        p.bound[d + xs] = sk[xs];
      } else if (sk[xs] < sb[xs]) {
        p.bound[xs] = sk[xs];
      } else {
        ++p.problem->measures->gs_diverge_count;
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

  // update counting
  if (p.problem->rt->compute_boxcount) {
    ++p.problem->measures->boxcount;
  }

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

void local_disp(u64 pt, problem_param* problem, prec& disp)
{
  assert(problem != nullptr);

  prec              ldisp;
  grow_shrink_param gsp;

  gsp.problem = problem;

  for (i32 axis = 0; axis < problem->pts.dimensions; ++axis) {
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
  // assert(problem->domain_bound[0] == 0);
  assert(problem->pts.size() <= INT64_MAX);

  problem->pts_end_idx = problem->pts.size() - 1;
  ga.disp              = 0;

  // rebuild index
  sort(problem->pts, problem->idx);

  // compute local dispersion and update global dispersion
  for (u64 j = 0; j < problem->pts.size(); ++j) {
    local_disp(j, problem, ga.disp);
  }
};

void return_bound(const program_param& rt)
{
  i8 ndel = ' ';

  *rt.os << "#d";

  if (rt.compute_disp) {
    *rt.os << ndel << 0;
    ndel = rt.delimiter;
  }
  if (rt.compute_ndisp) {
    *rt.os << ndel << 0;
    ndel = rt.delimiter;
  }
  if (rt.compute_boxcount) {
    *rt.os << ndel << 0;
    ndel = rt.delimiter;
  }

  if (rt.compute_disp) {
    *rt.os << ndel << INFINITY;
  }
  if (rt.compute_ndisp) {
    *rt.os << ndel << INFINITY;
  }
  if (rt.compute_boxcount) {
    *rt.os << ndel << INFINITY;
  }
  *rt.os << std::endl;
}

i32 return_results(const program_param&                       rt,
                   const std::vector<dptk::problem_measures>& measures)
{
  if (!rt.silent) {
    *rt.os << "# point set sequence size = " << measures.size() << std::endl;
  }

  if (rt.compute_disp || rt.compute_ndisp || rt.compute_boxcount) {
    if (!rt.silent) {
      for (u64 i = 0; i < measures.size(); ++i) {
        *rt.os << "# gs diverges = " << measures[i].gs_diverge_count << std::endl;
      }
    }
    if (!rt.silent) {
      i8 ndel = '(';
      *rt.os << "# ";
      if (rt.compute_disp) {
        *rt.os << ndel << "dispersion";
        ndel = ',';
      }
      if (rt.compute_ndisp) {
        *rt.os << ndel << "n*dispersion";
        ndel = ',';
      }
      if (rt.compute_boxcount) {
        *rt.os << ndel << "number of boxes";
      }
      *rt.os << ")" << std::endl;
    }

    return_bound(rt);

    for (u64 i = 0; i < measures.size(); ++i) {
      if (rt.compute_disp) {
        putsci(rt.os, measures[i].disp, 16);
        if (rt.compute_ndisp || rt.compute_boxcount)
          *rt.os << rt.delimiter;
      }
      if (rt.compute_ndisp) {
        putsci(rt.os, measures[i].ndisp, 16);
        if (rt.compute_boxcount)
          *rt.os << rt.delimiter;
      }
      if (rt.compute_boxcount) {
        *rt.os << measures[i].boxcount;
      }
      *rt.os << std::endl;
    }

    write_pointset_eos(rt.os);
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
  dptk::dispersion_param              gap;
  dptk::problem_param                 problem;
  dptk::program_param                 rt;
  dptk::i32                           r;
  dptk::ipointset_read_info           ipts_inf;
  std::vector<dptk::problem_measures> measures;

  // default configuration
  rt.compute_boxcount = false;
  rt.compute_disp     = false;
  rt.compute_ndisp    = false;
  rt.delimiter        = ' ';
  rt.del_use_ipts     = true;
  rt.silent           = false;
  rt.input            = "-";
  rt.output           = "-";
  problem.rt          = &rt;
  r                   = EXIT_SUCCESS;

  // parse arguments
  if (!dptk::parse_progargs(argc, argv, rt)) {
    return EXIT_FAILURE;
  }

  // initialize io streams
  dptk::istream_init(rt.input, rt.is);
  dptk::ostream_init(rt.output, rt.os);

  assert(rt.is != nullptr);
  assert(rt.os != nullptr);

  if (!rt.silent)
    *rt.os << "# src=" << rt.input << std::endl;

  // iterate through pointset sequence
  while (!rt.is->eof() && r == EXIT_SUCCESS) {
    // clear pointset
    problem.pts.clear();

    // retrieve point set
    dptk::read_pointset(*rt.is, problem.pts, &ipts_inf);

    // skip empty point sets
    if (problem.pts.coords.empty()) {
      continue;
    }

    assert(problem.pts.dimensions == 2);
    dptk::forward_delimiter(rt.del_use_ipts, ipts_inf, rt.delimiter);

    // allocate index
    problem.idx.allocate(problem.pts.size(), problem.pts.dimensions);

    // allocate measures
    measures.resize(measures.size() + 1);
    problem.measures           = &measures.back();
    problem.measures->boxcount = 0;

    // compute dispersion
    dptk::dispersion(gap, &problem);

    // store measurements
    problem.measures->disp  = gap.disp;
    problem.measures->ndisp = problem.pts.size() * gap.disp;
  }

  // show result
  r = dptk::return_results(rt, measures);

  // clean up (heap allocations)
  dptk::istream_close(rt.is);
  dptk::ostream_close(rt.os);

  return r;
}