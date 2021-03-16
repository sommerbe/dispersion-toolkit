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
  u1            silent;
  i8            delimiter;
  u1            del_use_ipts;
  u1            compute_disp;
  u1            compute_ndisp;
  u1            compute_boxcount;
  // u1            compute_boxes;
  // u1            compute_box_max;
  // u1            compute_box_interior;
};

struct problem_measures
{
  prec     disp;
  prec     ndisp;
  u64      boxcount;
  hyperbox box_max;
};

struct problem_param
{
  pointset               pts;
  pointset_dsorted_index idx;
  std::vector<hyperbox>  boxes;
  problem_measures*      measures;
  program_param*         rt;
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

void dispersion_naamad(problem_param* p)
{
  assert(p != nullptr);
  assert(p->pts.size() > 0);
  assert(p->pts.domain_bound.size() == 4);
  assert(p->measures != nullptr);
  assert(p->measures->boxcount == 0);

  prec*    p0;
  prec*    p1;
  prec*    p2;
  prec*    pj;
  hyperbox box;
  hyperbox box_max;
  prec     dist;
  prec     dist_max;
  prec     extent_up;
  prec     extent_low;

  // sort points along 0-axis
  // sort points along 1-axis
  sort(p->pts, p->idx);

  // find greatest empty distance along 0-axis
  p0       = p->pts.at(p->idx.at(0, 0), 0);
  dist_max = p0[0] - p->pts.domain_low(0);
  for (u64 i = 1; i < p->pts.size(); ++i) {
    p1   = p->pts.at(p->idx.at(i, 0), 0);
    dist = p1[0] - p0[0];
    assert(dist >= 0.);
    dist_max = std::max(dist_max, dist);
    p0       = p1;
  }
  dist_max     = std::max(dist_max, p->pts.domain_up(0) - p0[0]);
  box_max.area = dist_max * p->pts.domain_extent(1);
  p->measures->boxcount += p->pts.size() + 1;

  // find greatest empty rectangle (quasi interior) iterating along 0-axis
  // - sorted along 1-axis descending order (aka reverse ascending order)
  // - p0 is above p1 along 1-axis
  // - extent_low = T_l; extent_up = T_r (in the official paper)
  // - p0[0] = X_i, p1[0] = X_j; p0[1] = Y_i, p1[1] = Y_j (official paper)
  for (u64 i = p->pts.size(); i > 0;) {
    p0         = p->pts.at(p->idx.at(--i, 1), 0);
    extent_low = p->pts.domain_low(0);
    extent_up  = p->pts.domain_up(0);
    assert(extent_up >= extent_low);
    for (u64 j = i; j > 0;) {
      p1 = p->pts.at(p->idx.at(--j, 1), 0);
      if (p1[0] <= extent_low || p1[0] >= extent_up) {
        continue;
      }
      dist = p0[1] - p1[1];
      assert(dist >= 0.);
      box.area     = (extent_up - extent_low) * dist;
      box_max.area = std::max(box_max.area, box.area);

      // shrink extent along 0-axis
      if (p1[0] > p0[0]) {
        extent_up = p1[0];
      } else {
        extent_low = p1[0];
      }
      ++p->measures->boxcount;
    }
    dist         = p0[1] - p->pts.domain_low(1);
    box.area     = (extent_up - extent_low) * dist;
    box_max.area = std::max(box_max.area, box.area);
    ++p->measures->boxcount;
  }

  // find greatest empty rectangle which is
  // a) bounded by point i
  // b) bounded by top domain boundary
  // - sorted along 0-axis (for efficiency)
  // - find next (0-axis) neighbour of point i being above i (along 1-axis)
  for (u64 i = 0; i < p->pts.size(); ++i) {
    p0 = p->pts.at(p->idx.at(i, 0), 0);
    p1 = nullptr;
    p2 = nullptr;

    // find next right neighour (greater 0-axis coord)
    for (u64 j = i + 1; j < p->pts.size(); ++j) {
      pj = p->pts.at(p->idx.at(j, 0), 0);
      if (pj[0] > p0[0] && pj[1] > p0[1]) {
        p2 = pj;
        break;
      }
    }
    if (p2 == nullptr) {
      p2 = &p->pts.domain_bound[p->pts.domain_idx_up(0)];
    }

    // find next left neighbour (smaller 0-axis coord)
    for (u64 j = i; j > 0;) {
      pj = p->pts.at(p->idx.at(--j, 0), 0);
      if (pj[0] < p0[0] && pj[1] > p0[1]) {
        p1 = pj;
        break;
      }
    }
    if (p1 == nullptr) {
      p1 = &p->pts.domain_bound[p->pts.domain_idx_low(0)];
    }

    // accumulate empty box area
    box.area = (p2[0] - p1[0]) * (p->pts.domain_up(1) - p0[1]);
    assert(box.area >= 0.);
    box_max.area = std::max(box_max.area, box.area);
    ++p->measures->boxcount;
  }

  // update computation tasks
  p->measures->disp = box_max.area;

  // if (p->rt->compute_box_max) {
  //   p->measures->box_max = box_max;
  // }
};



i32 return_results(const program_param&                       rt,
                   const std::vector<dptk::problem_measures>& measures)
{
  putparam(rt.os, "point set sequence size", measures.size(), !rt.silent);

  // greatest box
  // if (rt.compute_box_max) {
  //   putln(rt.os, "# first greatest empty box:", !rt.silent);
  //   putln(rt.os,
  //         "# (low_0, low_1, ..., low_d, up_0, up_1, ..., up_d) for d+1 dimensions:",
  //         !rt.silent);

  //   for (u64 i = 0; i < measures.size(); ++i) {
  //     print_coords(rt.os, measures[i].box_max);
  //   }
  //   write_pointset_eos(rt.os);
  // }

  // dispersion, boxcount
  if (rt.compute_disp || rt.compute_ndisp || rt.compute_boxcount) {
    if (!rt.silent) {
      *rt.os << "# ";
      i8 ndel = '(';
      put_header_column(rt.os, "dispersion", ndel, ',', rt.compute_disp);
      put_header_column(rt.os, "n*dispersion", ndel, ',', rt.compute_ndisp);
      put_header_column(rt.os, "number of boxes", ndel, ',', rt.compute_boxcount);
      *rt.os << ")" << std::endl;
    }

    pointset pts;

    pts.clear();
    pts.append_domain_bound(0, INFINITY, rt.compute_disp);
    pts.append_domain_bound(0, INFINITY, rt.compute_ndisp);
    pts.append_domain_bound(0, INFINITY, rt.compute_boxcount);

    write_pointset_header(rt.os, pts, rt.delimiter);

    for (u64 i = 0; i < measures.size(); ++i) {
      if (rt.compute_disp) {
        putsci(rt.os, measures[i].disp, 16);
        put(rt.os, rt.delimiter, rt.compute_ndisp || rt.compute_boxcount);
      }
      if (rt.compute_ndisp) {
        putsci(rt.os, measures[i].ndisp, 16);
        put(rt.os, rt.delimiter, rt.compute_boxcount);
      }
      if (rt.compute_boxcount) {
        *rt.os << measures[i].boxcount;
      }
      *rt.os << std::endl;
    }
    
    write_pointset_footer(rt.os, pts);
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
      // } else if (s == "--interior-boxes") {
      //   rt.compute_box_interior = true;
      // } else if (s == "--greatest-box") {
      //   rt.compute_box_max = true;
      // } else if (s == "--boxes") {
      //   rt.compute_boxes = true;
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
      std::cout << "NAME: compute dispersion with algorithm of Naamad et al. 1984"
                << std::endl;
      std::cout << "SYNOPSIS: [--i  FILE]  [--o FILE] [--disp] [--ndisp] [--count-boxes] "
                   "[--silent]"
                << std::endl;
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
  dptk::program_param                 rt;
  dptk::i32                           r;
  dptk::ipointset_read_info           ipts_inf;
  std::vector<dptk::problem_measures> measures;
  std::vector<dptk::problem_param>    problems;

  // default configuration
  rt.compute_boxcount = false;
  rt.compute_disp     = false;
  rt.compute_ndisp    = false;
  // rt.compute_box_interior = false;
  // rt.compute_boxes        = false;
  // rt.compute_box_max      = false;
  rt.delimiter    = ' ';
  rt.del_use_ipts = true;
  rt.silent       = false;
  rt.input        = "-";
  rt.output       = "-";
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

  // show parameters
  dptk::putparam(rt.os, "compute disp", rt.compute_disp, !rt.silent);
  dptk::putparam(rt.os, "compute n*disp", rt.compute_ndisp, !rt.silent);
  dptk::putparam(rt.os, "compute number of boxes", rt.compute_boxcount, !rt.silent);
  dptk::putparam(rt.os, "delimiter", rt.delimiter, !rt.silent);
  dptk::putparam(rt.os, "source", rt.input, !rt.silent);

  // iterate through pointset sequence
  while (!rt.is->eof() && r == EXIT_SUCCESS) {
    // allocate problem
    dptk::problem_param problem;

    problem.rt = &rt;

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

    // allocate index
    problem.idx.allocate(problem.pts.size(), problem.pts.dimensions);

    problems.push_back(problem);
  }

  // allocate measures
  measures.resize(problems.size());

  // post allocation (due to pointer invalidation on dynamic vector)
  for (dptk::u64 i = 0; i < problems.size(); ++i) {
    // linking
    problems[i].measures = &measures[i];

    // default values
    problems[i].measures->boxcount = 0;
  }

// parallel compute dispersion
#pragma omp parallel for
  for (dptk::i64 i = 0; i < problems.size(); ++i) {
    // compute dispersion
    dptk::dispersion_naamad(&problems[i]);

    // store measurements
    problems[i].measures->ndisp = problems[i].pts.size() * problems[i].measures->disp;
  }

  // show result
  r = dptk::return_results(rt, measures);

  // clean up (heap allocations)
  dptk::istream_close(rt.is);
  dptk::ostream_close(rt.os);

  return r;
}
