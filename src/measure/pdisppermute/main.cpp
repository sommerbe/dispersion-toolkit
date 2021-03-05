#include "../../io/argparse.hpp"
#include "../../io/ipointset.hpp"
#include "../../io/opointset.hpp"
#include "../../io/ostream.hpp"
#include "../../math/arithmetic.hpp"
#include "../../math/pointset.hpp"
#include <float.h>
#include <iomanip>
#include <iostream>
#include <stdlib.h>

namespace dptk {

typedef b64                    prec;
typedef regular_pointset<prec> pointset;

const prec EQUALITY_THRESHOLD = 1e-10;

struct hyperbox
{
  // order: (coordinates of lower point = low, coordinates of upper point = up)
  // i.e.: (low_0, low_1, ..., low_d, up_0, up_1, ..., up_d) for d+1 dimensions
  std::vector<prec> coords;

  // the area of the hyperbox
  prec area;
};

struct program_param
{
  std::string   input;
  std::string   output;
  std::ostream* os;
  std::istream* is;
  u64           p;
  u1            silent;
  i8            delimiter;
  u1            del_use_ipts;
  u1            compute_pdisp;
  u1            compute_npdisp;
  u1            debug_permutations;
  u1            precompute_distances;
};

struct problem_measures
{
  prec disp;
  prec ndisp;
};

struct problem_param
{
  pointset          pts;
  prec              domain_bound[2];
  std::vector<u64>  permute;
  std::vector<prec> distance_matrix;
  problem_measures* measures;
  program_param*    rt;
};

void print_coords(std::ostream* os, const hyperbox& box, u8 del = ' ')
{
  ensure_precision(os, box.coords[0]);
  for (u32 i = 0; i < box.coords.size(); ++i) {
    if (i > 0) {
      *os << del;
    }
    *os << box.coords[i];
  }
  *os << std::endl;
}

void print_permute_indices(std::ostream* os, const std::vector<u64>& p, u8 del = ',')
{
  *os << "# ";
  for (u64 i = 0; i < p.size(); ++i) {
    if (i > 0) {
      *os << del;
    }
    *os << p[i];
  }
  *os << std::endl;
}

prec compute_distance(prec* a, prec* b, u32 dimensions)
{
  prec v = 0.0;
  prec s;

  for (u32 d = 0; d < dimensions; ++d) {
    s = a[d] - b[d];
    v += s * s;
  }
  v = std::sqrt(v);
  return v;
}

u64 distance_matrix_index(u64 n, u64 i, u64 j)
{
  return i * n + j;
}

void precompute_distance_matrix(problem_param* pb)
{
  u64   n = pb->pts.size();
  u64   k;
  prec* a;
  prec* b;
  prec  d;

  pb->distance_matrix.resize(n * n);

  for (u64 i = 0; i < n; ++i) {
    for (u64 j = i + 1; j < n; ++j) {
      a                      = pb->pts.at(i, 0);
      b                      = pb->pts.at(j, 0);
      d                      = compute_distance(a, b, pb->pts.dimensions);
      k                      = distance_matrix_index(n, i, j);
      pb->distance_matrix[k] = d;
      k                      = distance_matrix_index(n, j, i);
      pb->distance_matrix[k] = d;
    }
  }
}

void pdispersion_permute(problem_param* pb)
{
  assert(pb != nullptr);
  assert(pb->pts.size() > 0);

  std::vector<u64> permute;
  u64              n = 4; // pb->pts.size();
  u64              p = 2; // pb->rt->p;
  u64              e = n - p;

  // allocate an array of permutation indices
  permute.resize(p);

  // initialise permutation list with (0,1,2,...,p-1) indices
  for (u64 i = 0; i < permute.size(); ++i) {
    permute[i] = i;
  }

  // move permute indices (increasing order) up to end of point set array
  // so 4 pts: (0,1) > (0,2) > (0,3) > (1,2) > (1,3) > (2,3)
  while (true) {
    //

    // debug
    if (pb->rt->debug_permutations)
      print_permute_indices(pb->rt->os, permute);

    // all permute indices are on end of point set array
    if (permute[0] == e) {
      break;
    }
  }
};

void pdispersion_permute_stack_leaf(problem_param* pb)
{
  prec  dist = DBL_MAX;
  prec* a;
  prec* b;
  prec  d;
  u64   e = pb->permute.size() - 1;
  u64   n = pb->pts.size();
  u64   k;

  if (pb->rt->precompute_distances) {
    for (u64 i = 0; i < e; ++i) {
      for (u64 j = i + 1; j <= e; ++j) {
        k    = distance_matrix_index(n, pb->permute[i], pb->permute[j]);
        d    = pb->distance_matrix[k];
        dist = math::min(d, dist);
      }
    }
  } else {
    for (u64 i = 0; i < e; ++i) {
      for (u64 j = i + 1; j <= e; ++j) {
        a    = pb->pts.at(pb->permute[i], 0);
        b    = pb->pts.at(pb->permute[j], 0);
        d    = compute_distance(a, b, pb->pts.dimensions);
        dist = math::min(d, dist);
      }
    }
  }

  pb->measures->disp = math::max(dist, pb->measures->disp);

  // debug
  if (pb->rt->debug_permutations)
    print_permute_indices(pb->rt->os, pb->permute);
}

void pdispersion_permute_stack_loop(problem_param* pb, u32 pos, u32 i)
{
  pb->permute[pos] = i;

  if (pos > 0) {
    for (u64 j = i + 1; j < pb->pts.size(); ++j) {
      pdispersion_permute_stack_loop(pb, pos - 1, j);
    }
  } else {
    pdispersion_permute_stack_leaf(pb);
  }
}

void pdispersion_permute_stack(problem_param* pb)
{
  assert(pb != nullptr);
  assert(pb->pts.size() > 0);

  u64 p = math::min(pb->rt->p, pb->pts.size());

  pb->measures->disp = 0;
  pb->permute.resize(p);

  if (pb->rt->precompute_distances) {
    precompute_distance_matrix(pb);
  }

  for (u64 i = 0; i < pb->pts.size(); ++i) {
    pdispersion_permute_stack_loop(pb, p - 1, i);
  }

  pb->measures->ndisp = pb->pts.size() * pb->measures->disp;
}

void return_bound(const program_param& rt)
{
  i8 ndel = ' ';

  *rt.os << "#d";

  put_header_column(rt.os, 0, ndel, rt.delimiter, rt.compute_pdisp);
  put_header_column(rt.os, 0, ndel, rt.delimiter, rt.compute_npdisp);

  put_header_column(rt.os, INFINITY, ndel, rt.delimiter, rt.compute_pdisp);
  put_header_column(rt.os, INFINITY, ndel, rt.delimiter, rt.compute_npdisp);

  *rt.os << std::endl;
}

i32 return_partial_results(const program_param& rt, const problem_param& problem)
{

  return EXIT_SUCCESS;
}

i32 return_partial_results(const program_param&                       rt,
                           const std::vector<dptk::problem_measures>& measures)
{
  putparam(rt.os, "point set sequence size", measures.size(), !rt.silent);

  // dispersion, boxcount
  if (rt.compute_pdisp || rt.compute_npdisp) {
    if (!rt.silent) {
      *rt.os << "# ";
      i8 ndel = '(';
      put_header_column(rt.os, "dispersion", ndel, ',', rt.compute_pdisp);
      put_header_column(rt.os, "n*dispersion", ndel, ',', rt.compute_npdisp);
      *rt.os << ")" << std::endl;
    }

    return_bound(rt);

    for (u64 i = 0; i < measures.size(); ++i) {
      if (rt.compute_pdisp) {
        putsci(rt.os, measures[i].disp, 16);
        if (rt.compute_npdisp)
          *rt.os << rt.delimiter;
      }
      if (rt.compute_npdisp) {
        putsci(rt.os, measures[i].ndisp, 16);
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
      rt.compute_pdisp = true;
    } else if (s == "--ndisp") {
      rt.compute_npdisp = true;
    } else if (s == "--debug-permute") {
      rt.debug_permutations = true;
    } else if (s == "--p") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing p value. Consider using -h or --help.");
      rt.p = std::strtol(arg[++i].c_str(), nullptr, 10);
      if (arg[i][0] == '-' || rt.p < 2)
        return argparse::err("invalid argument: p > 1.");

    } else if (s == "--silent") {
      rt.silent = true;
    } else if (s == "--i") {
      if (++i == arg.size()) {
        std::cerr << "invalid argument: -i misses a mandatory parameter" << std::endl;
        return false;
      }
      rt.input = arg[i];
    } else if (s == "--o") {
      if (++i == arg.size()) {
        std::cerr << "invalid argument: -o misses a mandatory parameter" << std::endl;
        return false;
      }
      rt.output = arg[i];
    } else if (s == "-h" || s == "--help") {
      std::cout
        << "NAME: compute p-dispersion with a permutation algorithm (exhaustive search)"
        << std::endl;
      std::cout
        << "SYNOPSIS: [--i FILE] [--o FILE] [--p=2] [--disp] [--ndisp] [--debug-permute] [--silent]"
        << std::endl;
      return false;
    }
  }

  if (!(rt.compute_pdisp || rt.compute_npdisp)) {
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
  dptk::problem_param                 problem;
  dptk::program_param                 rt;
  dptk::i32                           r;
  dptk::ipointset_read_info           ipts_inf;
  std::vector<dptk::problem_measures> measures;

  // default configuration
  rt.compute_pdisp        = false;
  rt.compute_npdisp       = false;
  rt.delimiter            = ' ';
  rt.del_use_ipts         = true;
  rt.silent               = false;
  rt.input                = "-";
  rt.output               = "-";
  rt.p                    = 2;
  rt.debug_permutations   = false;
  rt.precompute_distances = true;
  problem.rt              = &rt;
  problem.domain_bound[0] = 0;
  problem.domain_bound[1] = 1;
  r                       = EXIT_SUCCESS;

  // parse arguments
  if (!dptk::parse_progargs(argc, argv, rt)) {
    return EXIT_FAILURE;
  }

  // initialize io streams
  dptk::istream_init(rt.input, rt.is);
  dptk::ostream_init(rt.output, rt.os);

  assert(rt.is != nullptr);
  assert(rt.os != nullptr);

  // show parameters
  dptk::putparam(rt.os, "compute disp", rt.compute_pdisp, !rt.silent);
  dptk::putparam(rt.os, "compute n*disp", rt.compute_npdisp, !rt.silent);
  dptk::putparam(rt.os, "p", rt.p, !rt.silent);
  dptk::putparam(rt.os, "delimiter", rt.delimiter, !rt.silent);
  dptk::putparam(rt.os, "source", rt.input, !rt.silent);

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

    assert(problem.pts.dimensions == 2);
    dptk::forward_delimiter(rt.del_use_ipts, ipts_inf, rt.delimiter);

    // allocate measures
    measures.resize(measures.size() + 1);
    problem.measures = &measures.back();

    // compute dispersion
    // dptk::pdispersion_permute(&problem);
    dptk::pdispersion_permute_stack(&problem);

    // store measurements
    problem.measures->ndisp = problem.pts.size() * problem.measures->disp;

    // show result
    r = dptk::return_partial_results(rt, problem);
  }

  // show result
  r = dptk::return_partial_results(rt, measures);

  // clean up (heap allocations)
  dptk::istream_close(rt.is);
  dptk::ostream_close(rt.os);

  return r;
}
