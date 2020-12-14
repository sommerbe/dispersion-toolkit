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
  u1            compute_boxes;
  u1            compute_box_max;
  u1            compute_box_interior;
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
  pointset              pts;
  prec                  domain_bound[2];
  std::vector<hyperbox> boxes;
  problem_measures*     measures;
  program_param*        rt;
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

void set_rect_boundary(pointset& pts,
                       i32       index,
                       prec*     box_corner,
                       u32       axis,
                       prec      boundary_default)
{
  if (index == -1)
    box_corner[axis] = boundary_default;
  else
    box_corner[axis] = *pts.at(index, axis);
}

bool is_inside(const hyperbox& box, u32 dimensions, const prec* point)
{
  assert(u64(dimensions) * 2 == box.coords.size());

  for (u32 d = 0; d < dimensions; ++d) {
    if (point[d] <= box.coords[d] || point[d] >= box.coords[d + dimensions])
      return false;
  }

  return true;
}

prec compute_area(const hyperbox& box, u32 dimensions)
{
  assert(u64(dimensions) * 2 == box.coords.size());

  prec area = 1.0;

  for (u32 d = 0; d < dimensions; ++d) {
    assert(box.coords[d + dimensions] >= box.coords[d]);
    area *= box.coords[d + dimensions] - box.coords[d];
  }
  // exclude possibility of bit-error math
  // - logically, the previous assert should prevent this assert
  assert(area >= 0);
  return area;
}

void dispersion_combinatorial(problem_param* p)
{
  assert(p != nullptr);
  assert(p->pts.size() > 0);
  assert(p->domain_bound[0] == 0);
  assert(p->pts.size() <= INT64_MAX);

  const prec p_zero = 0.0;
  const prec p_one  = 1.0;

  prec     percentile[3];
  prec*    bd;
  hyperbox box;
  hyperbox box_max;
  i64      n;
  bool     inside;
  prec*    r_low;
  prec*    r_up;

  box.coords.resize(2 * p->pts.dimensions);
  p->boxes.clear();

  p->measures->disp     = 0;
  p->measures->boxcount = 0;
  box_max.area          = 0;
  n                     = p->pts.size();
  r_low                 = &box.coords[0];
  r_up                  = &box.coords[p->pts.dimensions];

  for (i32 d0low = -1; d0low < n; ++d0low) {
    set_rect_boundary(p->pts, d0low, r_low, 0, p_zero);
    for (i32 d0up = -1; d0up < n; ++d0up) {
      set_rect_boundary(p->pts, d0up, r_up, 0, p_one);
      // check: boundary condition: order of rect coordinates
      if (r_up[0] <= r_low[0]) // || d0up == d0low)
        continue;

      // dimension: 1
      for (i32 d1low = -1; d1low < n; ++d1low) {
        // check boundary condition: constraint by dim 0 to left/right
        if (d1low >= 0
            && ((d0low >= 0
                 && (*p->pts.at(d1low, 0) < *p->pts.at(d0low, 0)
                     || *p->pts.at(d1low, 1) > *p->pts.at(d0low, 1)))
                || (d0up >= 0
                    && (*p->pts.at(d1low, 0) > *p->pts.at(d0up, 0)
                        || *p->pts.at(d1low, 1) > *p->pts.at(d0up, 1)))))
          continue;

        set_rect_boundary(p->pts, d1low, r_low, 1, p_zero);

        for (i32 d1up = -1; d1up < n; ++d1up) {
          // check boundary condition: constraint by dim 0 to left/right
          if (d1up >= 0
              && ((d0low >= 0
                   && (*p->pts.at(d1up, 0) < *p->pts.at(d0low, 0)
                       || *p->pts.at(d1up, 1) < *p->pts.at(d0low, 1)))
                  || (d0up >= 0
                      && (*p->pts.at(d1up, 0) > *p->pts.at(d0up, 0)
                          || *p->pts.at(d1up, 1) < *p->pts.at(d0up, 1)))
                  || (d1low >= 0 && *p->pts.at(d1up, 1) <= *p->pts.at(d1low, 1))))
            continue;

          set_rect_boundary(p->pts, d1up, r_up, 1, p_one);

          // check: boundary condition: order of rect coordinates
          if (r_up[1] <= r_low[1]) // || d1up == d1low)
            continue;

          // check: emptiness condition of rectangular
          inside = false;
          for (i32 s = 0; s < n && !inside; ++s) {
            inside = is_inside(box, p->pts.dimensions, p->pts.at(s, 0));
          }
          if (inside)
            continue;

          // update: area
          box.area = compute_area(box, p->pts.dimensions);

          // update: greatest empty rectangular
          if (box.area > box_max.area) {
            box_max = box;
          }

          // update computation tasks
          if (p->rt->compute_boxes
              && (!p->rt->compute_box_interior
                  || (d0low >= 0 && d0up >= 0 && d1low >= 0 && d1up >= 0))) {
            p->boxes.push_back(box);
          }

          if (p->rt->compute_boxcount) {
            ++p->measures->boxcount;
          }
        }
      }
    }
  }

  // update computation tasks
  p->measures->disp = box_max.area;

  if (p->rt->compute_box_max) {
    p->measures->box_max = box_max;
  }
};

void return_bound(const program_param& rt)
{
  i8 ndel = ' ';

  *rt.os << "#d";

  put_header_column(rt.os, 0, ndel, rt.delimiter, rt.compute_disp);
  put_header_column(rt.os, 0, ndel, rt.delimiter, rt.compute_ndisp);
  put_header_column(rt.os, 0, ndel, rt.delimiter, rt.compute_boxcount);

  put_header_column(rt.os, INFINITY, ndel, rt.delimiter, rt.compute_disp);
  put_header_column(rt.os, INFINITY, ndel, rt.delimiter, rt.compute_ndisp);
  put_header_column(rt.os, INFINITY, ndel, rt.delimiter, rt.compute_boxcount);

  *rt.os << std::endl;
}

i32 return_partial_results(const program_param& rt, const problem_param& problem)
{
  if (rt.compute_boxes) {
    if (rt.compute_box_interior) {
      putln(rt.os, "# all interior empty box:", !rt.silent);
    } else {
      putln(rt.os, "# all empty box:", !rt.silent);
    }
    putln(rt.os,
          "# (low_0, low_1, ..., low_d, up_0, up_1, ..., up_d) for d+1 dimensions:",
          !rt.silent);
    for (u64 i = 0; i < problem.boxes.size(); ++i) {
      print_coords(rt.os, problem.boxes[i]);
    }
    write_pointset_eos(rt.os);
  }

  return EXIT_SUCCESS;
}

i32 return_partial_results(const program_param&                       rt,
                           const std::vector<dptk::problem_measures>& measures)
{
  putparam(rt.os, "point set sequence size", measures.size(), !rt.silent);

  // greatest box
  if (rt.compute_box_max) {
    putln(rt.os, "# first greatest empty box:", !rt.silent);
    putln(rt.os,
          "# (low_0, low_1, ..., low_d, up_0, up_1, ..., up_d) for d+1 dimensions:",
          !rt.silent);

    for (u64 i = 0; i < measures.size(); ++i) {
      print_coords(rt.os, measures[i].box_max);
    }
    write_pointset_eos(rt.os);
  }

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
    } else if (s == "--interior-boxes") {
      rt.compute_box_interior = true;
    } else if (s == "--greatest-box") {
      rt.compute_box_max = true;
    } else if (s == "--boxes") {
      rt.compute_boxes = true;
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
      std::cout << manpage;
      return false;
    }
  }

  if (!(rt.compute_disp || rt.compute_ndisp || rt.compute_boxcount || rt.compute_box_max
        || rt.compute_boxes)) {
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
  rt.compute_boxcount     = false;
  rt.compute_disp         = false;
  rt.compute_ndisp        = false;
  rt.compute_box_interior = false;
  rt.compute_boxes        = false;
  rt.compute_box_max      = false;
  rt.delimiter            = ' ';
  rt.del_use_ipts         = true;
  rt.silent               = false;
  rt.input                = "-";
  rt.output               = "-";
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
    dptk::putparam(rt.os, "compute disp", rt.compute_disp, !rt.silent);
    dptk::putparam(rt.os, "compute n*disp", rt.compute_ndisp, !rt.silent);
    dptk::putparam(rt.os, "compute number of boxes", rt.compute_boxcount, !rt.silent);
    dptk::putparam(rt.os, "compute interior boxes", rt.compute_box_interior, !rt.silent);
    dptk::putparam(rt.os, "compute boxes", rt.compute_boxes, !rt.silent);
    dptk::putparam(rt.os, "compute greatest box", rt.compute_box_max, !rt.silent);
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
    problem.measures           = &measures.back();
    problem.measures->boxcount = 0;

    // compute dispersion
    dptk::dispersion_combinatorial(&problem);

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
