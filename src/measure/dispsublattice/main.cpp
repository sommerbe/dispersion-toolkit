#include "../../index/pointset_dsorted.hpp"
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

template<typename any>
struct hyperbox
{
  // order: (coordinates of lower point = low, coordinates of upper point = up)
  // i.e.: (low_0, low_1, ..., low_d, up_0, up_1, ..., up_d) for d+1 dimensions
  std::vector<any> coords;

  // the area of the hyperbox
  prec area;
};

struct program_param
{
  std::string    input;
  std::string    output;
  std::ostream*  os;
  std::istream*  is;
  u1             silent;
  i8             delimiter;
  u1             del_use_ipts;
  u1             compute_disp;
  u1             compute_ndisp;
  u1             compute_boxcount;
  u1             compute_boxes;
  u1             compute_box_max;
  u1             compute_box_areas;
  u1             compute_box_coords;
  u1             layout_graph;
  u1 debug;
  prec           box_area_min;
  prec           box_area_max;
  hyperbox<prec> subdomain;
  u64            refptidx;
};

struct problem_measures
{
  prec           disp;
  prec           ndisp;
  u64            boxcount;
  hyperbox<prec> box_max;
};

struct problem_param
{
  pointset                    pts;
  prec                        domain_bound[2];
  std::vector<hyperbox<prec>> boxes;
  problem_measures*           measures;
  const program_param*        rt;

  // internal
  pointset_dsorted_index psort_idx;
  hyperbox<prec>         box;
  u64                    boxspan_ref;
  u64                    boxspan_target;
  u32                    box_exteriour;
  u64 exclude_d[2];
};

void print_coords(std::ostream*         os,
                  const hyperbox<prec>& box,
                  u1                    pcoords,
                  u1                    parea,
                  u8                    del = ' ')
{
  if (!pcoords && !parea)
    return;

  ensure_precision(os, box.coords[0]);

  if (pcoords) {
    for (u32 i = 0; i < box.coords.size(); ++i) {
      if (i > 0) {
        *os << del;
      }
      *os << box.coords[i];
    }
  }

  if (parea) {
    if (pcoords) {
      *os << del;
    }
    *os << box.area;
  }

  *os << std::endl;
}

prec compute_area(const hyperbox<prec>& box, u32 dimensions)
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

u1 predicate_area(prec area, prec range_min, prec range_max)
{
  return range_min <= area && area < range_max;
}

/**
 * checks whether point is inside the hyperbox, excluding its boundary
 */
bool is_inside(const hyperbox<prec>& box, u32 dimensions, const prec* point)
{
  assert(u64(dimensions) * 2 == box.coords.size());

  for (u32 d = 0; d < dimensions; ++d) {
    if (point[d] <= box.coords[d] || point[d] >= box.coords[d + dimensions])
      return false;
  }
  return true;
}

/**
 * checks whether point is inside the hyperbox, excluding its boundary
 */
bool is_inside(const hyperbox<prec>& box, u32 dimensions, u32 exclude_d, const prec* point)
{
  assert(u64(dimensions) * 2 == box.coords.size());

  for (u32 d = 0; d < dimensions; ++d) {
    if (d == exclude_d)
      continue;
    if (point[d] <= box.coords[d] || point[d] >= box.coords[d + dimensions])
      return false;
  }
  return true;
}

/**
 * checks whether point is outside the hyperbox and its boundary.
 *
 * note: is_outside(.) != !is_inside(.)
 */
bool is_outside(const hyperbox<prec>& box, u32 dimensions, const prec* point)
{
  assert(u64(dimensions) * 2 == box.coords.size());

  for (u32 d = 0; d < dimensions; ++d) {
    if (point[d] < box.coords[d] || point[d] > box.coords[d + dimensions])
      return true;
  }
  return false;
}

void extend_hyperbox_axis(problem_param* p, u64 d)
{
  // extend hyperbox p.box maximally along axis index d
  if (d < p->pts.dimensions) {
    // extend to below
    if (d != p->exclude_d[0]) {
      u32 z = 1;
      b64 y = p->pts.domain_low(d);
      for (u64 i=0; i<p->pts.size(); ++i) {
        prec* pt = p->pts.at(i,0);
        if (pt[d] >= p->box.coords[d] || !is_inside(p->box, p->pts.dimensions, d, pt)) {
          continue;
        }
        y = math::max(pt[d], y);
        z = 0;
      }
      if (z == 0) {
        p->box.coords[d] = y;
      }
      p->box_exteriour += z;
    }

    // extend to above
    if (d != p->exclude_d[1]) {
      u32 z = 1;
      u64 di = d + p->pts.dimensions;
      b64 y = p->pts.domain_up(d);
      for (u64 i=0; i<p->pts.size(); ++i) {
        prec* pt = p->pts.at(i,0);
        if (pt[d] <= p->box.coords[di] || !is_inside(p->box, p->pts.dimensions, d, pt)) {
          continue;
        }
        y = math::min(pt[d], y);
        z = 0;
      }
      if (z == 0) {
        p->box.coords[di] = y;
      }
      p->box_exteriour += z;
    }

    // find indices of both box corner points (ref, target)
    // u64 refd = p->psort_idx.search(p->boxspan_ref, d);
    // u64 trgd = p->psort_idx.search(p->boxspan_target, d);
    // u64 il   = std::min(refd, trgd);
    // u64 iu   = std::max(refd, trgd);
    // // check: need to be extendable to below and to above
    // if (il == 0 || iu + 1 == p->psort_idx.stride) {
    //   ++p->box_exteriour;
    //   return;
    // }
    // // extend to below and to above
    // p->box.coords[d]                     = *p->pts.at(p->psort_idx.at(il - 1, d), d);
    // p->box.coords[d + p->pts.dimensions] = *p->pts.at(p->psort_idx.at(iu + 1, d), d);

    // extend along next axis
    return extend_hyperbox_axis(p, d + 1);
  }

  // check: being interiour hyperbox
  // assert(p->box_exteriour == 0);

  // check: emptiness condition of hyperbox
  for (u64 i = 0; i < p->pts.size(); ++i) {
    assert(!is_inside(p->box, p->pts.dimensions, p->pts.at(i, 0)));
  }

  // compute area of hyperbox (axes-aligned)
  p->box.area = compute_area(p->box, p->pts.dimensions);

  // determine greatest hyperbox
  if (p->box.area > p->measures->box_max.area) {
    p->measures->box_max = p->box;
  }

  // protocol hyperboxes
  if (p->rt->compute_boxes
      && predicate_area(p->box.area, p->rt->box_area_min, p->rt->box_area_max)) {
    p->boxes.push_back(p->box);
  }

  // count number of empty hyperboxes
  if (p->rt->compute_boxcount) {
    ++p->measures->boxcount;
  }
}

void dispersion_subdomain_recursive(problem_param* p)
{
  assert(p != nullptr);
  assert(p->pts.size() > 0);
  assert(p->rt->subdomain.coords.size() == 2 * p->pts.dimensions);
  assert(p->rt->refptidx < p->pts.size());

  p->box.coords.resize(2 * p->pts.dimensions);

  p->box_exteriour      = 0;
  p->box.area           = 0;
  p->measures->box_max  = p->box;
  p->measures->disp     = 0;
  p->measures->boxcount = 0;
  p->boxes.clear();

  prec* pt[2];
  hyperbox<prec> box_base;

  // preprocess a sorted list of points, along each axis independently
  // - sorting halfs the quadratic complexity during each recursive step
  p->psort_idx.allocate(p->pts.size(), p->pts.dimensions);
  sort(p->pts, p->psort_idx);

  // reference point index
  p->boxspan_ref = p->rt->refptidx;

  // find target point within subdomain
  for (u64 i = 0; i < p->pts.size(); ++i) {
    // skip null hyperbox
    if (i == p->boxspan_ref) {
      continue;
    }
    // skip points being outside specified subdomain
    if (is_outside(p->rt->subdomain, p->pts.dimensions, p->pts.at(i, 0))) {
      continue;
    }
    // construct a temporary hyperbox
    pt[0] = p->pts.at(i, 0);
    pt[1] = p->pts.at(p->boxspan_ref, 0);
    for (u64 d = 0; d < p->pts.dimensions; ++d) {
      p->box.coords[d]                     = math::min(pt[0][d], pt[1][d]);
      p->box.coords[d + p->pts.dimensions] = math::max(pt[0][d], pt[1][d]);
    }
    // check: emptyness condition
    u64 j;
    for (j = 0; j < p->pts.size(); ++j) {
      if (j == i || j == p->boxspan_ref) {
        continue;
      }
      if (is_inside(p->box, p->pts.dimensions, p->pts.at(j, 0))) {
        break;
      }
    }
    if (j < p->pts.size()) {
      continue;
    }    
    putparam(p->rt->os, "considering box candidate", p->box.coords, p->rt->debug);

    // recursive search of hyperboxes (interiour hyperboxes only)
    // - start with axis 0 -> d-1; at d, record hyperbox
    p->boxspan_target = i;
    box_base = p->box;
    for (u64 d0 = 0; d0 < p->pts.dimensions; ++d0) {
      for (u64 d1 = 0; d1 < p->pts.dimensions; ++d1) {
        putparam(p->rt->os, "exclude axes", std::vector<u64>({d0,d1}), p->rt->debug);
        p->box.coords = box_base.coords;
        p->exclude_d[0] = d0;
        p->exclude_d[1] = d1;
        extend_hyperbox_axis(p, 0);
        putparam(p->rt->os, "found empty box", p->box.coords, p->rt->debug);
      }
    }
  }
  putparam(p->rt->os, "greatest box", p->measures->box_max.coords, p->rt->debug);

  // compute dispersion
  p->measures->disp = p->measures->box_max.area;
}

i32 return_partial_results(const program_param& rt, const problem_param& problem)
{
  if (rt.compute_boxes) {
    putln(rt.os, "# all empty box:", !rt.silent);
    putln(rt.os,
          "# (low_0, low_1, ..., low_d, up_0, up_1, ..., up_d) for d+1 dimensions:",
          !rt.silent);
    for (u64 i = 0; i < problem.boxes.size(); ++i) {
      print_coords(rt.os,
                   problem.boxes[i],
                   rt.compute_box_coords,
                   rt.compute_box_areas,
                   rt.delimiter);
    }
    write_pointset_eos(rt.os);
  }

  putln(rt.os,
        "# encountered exteriour boxes: ",
        problem.box_exteriour,
        problem.box_exteriour > 0);

  return EXIT_SUCCESS;
}

i32 return_partial_results(const program_param&                       rt,
                           const std::vector<dptk::problem_measures>& measures,
                           const std::vector<dptk::problem_param>&    problems)
{
  putparam(rt.os, "point set sequence size", measures.size(), !rt.silent);

  // greatest box
  if (rt.compute_box_max) {
    putln(rt.os, "# first greatest empty box:", !rt.silent);
    putln(rt.os,
          "# (low_0, low_1, ..., low_d, up_0, up_1, ..., up_d) for d+1 dimensions:",
          !rt.silent);

    for (u64 i = 0; i < measures.size(); ++i) {
      print_coords(rt.os,
                   measures[i].box_max,
                   rt.compute_box_coords,
                   rt.compute_box_areas,
                   rt.delimiter);
    }
    write_pointset_eos(rt.os);
  }

  // dispersion, boxcount
  if (rt.compute_disp || rt.compute_ndisp || rt.compute_boxcount) {
    if (!rt.silent) {
      *rt.os << "# ";
      i8 ndel = '(';
      put_header_column(rt.os, "argument", ndel, ',', rt.layout_graph);
      put_header_column(rt.os, "dispersion", ndel, ',', rt.compute_disp);
      put_header_column(rt.os, "n*dispersion", ndel, ',', rt.compute_ndisp);
      put_header_column(rt.os, "number of boxes", ndel, ',', rt.compute_boxcount);
      *rt.os << ")" << std::endl;
    }

    pointset pts;

    pts.clear();
    pts.append_domain_bound(0, INFINITY, rt.layout_graph);
    pts.append_domain_bound(0, INFINITY, rt.compute_disp);
    pts.append_domain_bound(0, INFINITY, rt.compute_ndisp);
    pts.append_domain_bound(0, INFINITY, rt.compute_boxcount);

    write_pointset_header(rt.os, pts, rt.delimiter);

    for (u64 i = 0; i < measures.size(); ++i) {
      if (rt.layout_graph) {
        if (problems[i].pts.arguments.empty()) {
          *rt.os << "0";
        } else {
          putsci(rt.os, problems[i].pts.arguments[0], 16);
        }
        *rt.os << rt.delimiter;
      }
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
    } else if (s == "--greatest-box") {
      rt.compute_box_max = true;
    } else if (s == "--boxes") {
      rt.compute_boxes = true;
    } else if (s == "--box-area") {
      rt.compute_box_areas = true;
    } else if (s == "--no-box-coordinates") {
      rt.compute_box_coords = false;
    } else if (s == "--debug") {
      rt.debug = true;

    } else if (s == "--box-area-min") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing box-area-min value. Consider using -h or --help.");
      rt.box_area_min = std::strtod(arg[++i].c_str(), nullptr);

    } else if (s == "--box-area-max") {
      if (!argparse::argval(arg, i))
        return argparse::err("missing box-area-max value. Consider using -h or --help.");
      rt.box_area_max = std::strtod(arg[++i].c_str(), nullptr);

    } else if (s == "--graph-layout") {
      rt.layout_graph = true;

    } else if (s == "--subdomain") {
      if (!argparse::retrieve(arg, i, rt.subdomain.coords) || rt.subdomain.coords.empty())
        return argparse::err("missing subdomain value. Consider using -h or --help.");
      if (rt.subdomain.coords.size() % 2 != 0)
        return argparse::err(
          "invalid subdomain value. The size of this list is not a multiple of 2.");

    } else if (s == "--refpt-index") {
      if (!argparse::retrieve(arg, i, rt.refptidx))
        return argparse::err("missing refpt-index value. Consider using -h or --help.");

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
        << "NAME: compute dispersion with a combinatorial algorithm  (exhaustive search)"
        << std::endl;
      std::cout << "SYNOPSIS: [--i FILE] [--o FILE] --disp|--ndisp --subdomain BINARY64 "
                   ".. BINARY64 --refpt-index INTEGER  [--count-boxes] "
                   "[--boxes] [--greatest-box] [--box-area] "
                   "[--no-box-coordinates] [--box-area-min=BINARY64] "
                   "[--box-area-max=BINARY64] [--graph-layout] [--silent]"
                << std::endl;
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

  if (rt.subdomain.coords.empty()) {
    std::cerr << "missing option: --subdomain" << std::endl;
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
  rt.compute_boxcount   = false;
  rt.compute_disp       = false;
  rt.compute_ndisp      = false;
  rt.compute_boxes      = false;
  rt.compute_box_max    = false;
  rt.compute_box_areas  = false;
  rt.compute_box_coords = true;
  rt.layout_graph       = false;
  rt.delimiter          = ' ';
  rt.del_use_ipts       = true;
  rt.silent             = false;
  rt.box_area_min       = DBL_MIN;
  rt.box_area_max       = DBL_MAX;
  rt.input              = "-";
  rt.output             = "-";
  rt.refptidx           = 0;
  rt.debug = false;
  r                     = EXIT_SUCCESS;

  rt.subdomain.coords.clear();

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
  dptk::putparam(rt.os, "compute boxes", rt.compute_boxes, !rt.silent);
  dptk::putparam(rt.os, "compute greatest box", rt.compute_box_max, !rt.silent);
  dptk::putparam(rt.os, "subdomain", rt.subdomain.coords, !rt.silent);
  dptk::putparam(rt.os, "reference index", rt.refptidx, !rt.silent);
  dptk::putparam(rt.os, "predicate box area min", rt.box_area_min, !rt.silent);
  dptk::putparam(rt.os, "predicate box area max", rt.box_area_max, !rt.silent);
  dptk::putparam(rt.os, "stream box area", rt.compute_box_areas, !rt.silent);
  dptk::putparam(rt.os, "stream box coord", rt.compute_box_coords, !rt.silent);
  dptk::putparam(rt.os, "delimiter", rt.delimiter, !rt.silent);
  dptk::putparam(rt.os, "source", rt.input, !rt.silent);

  // iterate through pointset sequence
  while (!rt.is->eof() && r == EXIT_SUCCESS) {
    // allocate problem
    dptk::problem_param problem;

    problem.rt              = &rt;
    problem.domain_bound[0] = 0;
    problem.domain_bound[1] = 1;

    // clear pointset
    problem.pts.clear();

    // retrieve point set
    dptk::read_pointset(*rt.is, problem.pts, &ipts_inf);

    // skip empty points
    if (problem.pts.coords.empty()) {
      continue;
    }

    dptk::forward_delimiter(rt.del_use_ipts, ipts_inf, rt.delimiter);

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
    // dptk::dispersion_combinatorial(&problems[i]);
    dptk::dispersion_subdomain_recursive(&problems[i]);

    // store measurements
    problems[i].measures->ndisp = problems[i].pts.size() * problems[i].measures->disp;
  }

  // show sequence of interiour/exteriour/all/.. boxes
  for (dptk::u64 i = 0; i < problems.size(); ++i) {
    dptk::return_partial_results(rt, problems[i]);
  }

  // show sequence of dispersion, n*dispersion, box counts, greatest boxes
  r = dptk::return_partial_results(rt, measures, problems);

  // clean up (heap allocations)
  dptk::istream_close(rt.is);
  dptk::ostream_close(rt.os);

  return r;
}
