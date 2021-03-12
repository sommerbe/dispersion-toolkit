#pragma once
#include "types.hpp"
#include <cassert>
#include <vector>

namespace dptk {

template<typename prec>
struct regular_pointset
{
  u64               points;
  u64               dimensions;
  std::vector<prec> coords;

  // format bound: (low_0 ... low_(d-1) up_0 ... up_(d-1))
  std::vector<prec> domain_bound;

  void allocate(u64 num_points, u64 num_dimensions);
  void clear();

  void reset_inf_bound();

  u64 size() const;
  u1  empty() const;

  inline prec*       at(u64 pidx, u64 axis);
  inline const prec* at(u64 pidx, u64 axis) const;
  inline prec&       operator[](u64 pidx);

  void extract(u64 axis, regular_pointset<prec>& pts) const;

  void retrieve_points(u64 poffset, const regular_pointset<prec>& other);

  prec domain_extent(u64 axis) const;
  prec domain_low(u64 axis) const;
  prec domain_up(u64 axis) const;

  u64 domain_idx_low(u64 axis) const;
  u64 domain_idx_up(u64 axis) const;
};

template<typename prec>
u64 regular_pointset<prec>::size() const
{
  return points;
}

template<typename prec>
u1 regular_pointset<prec>::empty() const
{
  return coords.empty();
}

template<typename prec>
void regular_pointset<prec>::allocate(u64 num_points, u64 num_dimensions)
{
  points     = num_points;
  dimensions = num_dimensions;
  coords.resize(num_points * num_dimensions);
}

template<typename prec>
inline prec* regular_pointset<prec>::at(u64 pidx, u64 axis)
{
  assert(pidx < coords.size());
  assert(axis < dimensions);
  assert(pidx * dimensions + axis < coords.size());
  return &coords[pidx * dimensions + axis];
}

template<typename prec>
inline const prec* regular_pointset<prec>::at(u64 pidx, u64 axis) const
{
  assert(axis < dimensions);
  assert(pidx * dimensions + axis < coords.size());
  return &coords[pidx * dimensions + axis];
}

template<typename prec>
inline prec& regular_pointset<prec>::operator[](u64 pidx)
{
  return *at(pidx, 0);
}

template<typename prec>
void regular_pointset<prec>::clear()
{
  dimensions = 0;
  points     = 0;
  coords.clear();
}

template<typename prec>
void regular_pointset<prec>::reset_inf_bound()
{
  domain_bound.resize(2 * dimensions);
  for (u64 i = 0; i < dimensions; ++i) {
    domain_bound[i]              = -INFINITY;
    domain_bound[i + dimensions] = +INFINITY;
  }
}

template<typename prec>
void regular_pointset<prec>::extract(u64 axis, regular_pointset<prec>& pts) const
{
  pts.allocate(points, 1);
  for (u64 i = 0; i < points; ++i) {
    pts.coords[i] = *at(i, axis);
  }
}

template<typename prec>
void regular_pointset<prec>::retrieve_points(u64                           poffset,
                                             const regular_pointset<prec>& other)
{
  assert(dimensions == other.dimensions);

  for (u64 i = 0; i < other.points; ++i) {
    for (u64 j = 0; j < other.dimensions; ++j) {
      *at(i + poffset, j) = *other.at(i, j);
    }
  }
}

template<typename prec>
prec regular_pointset<prec>::domain_extent(u64 axis) const
{
  return domain_bound[dimensions + axis] - domain_bound[axis];
}

template<typename prec>
prec regular_pointset<prec>::domain_low(u64 axis) const
{
  return domain_bound[axis];
}

template<typename prec>
prec regular_pointset<prec>::domain_up(u64 axis) const
{
  return domain_bound[dimensions + axis];
}

template<typename prec>
u64 regular_pointset<prec>::domain_idx_low(u64 axis) const
{
  return axis;
}

template<typename prec>
u64 regular_pointset<prec>::domain_idx_up(u64 axis) const
{
  return dimensions + axis;
}

} // namespace dptk
