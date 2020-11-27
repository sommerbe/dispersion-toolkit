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

  void allocate(u64 num_points, u64 num_dimensions);
  void clear();

  u64 size() const;

  inline prec*       at(u64 pidx, u64 axis);
  inline const prec* at(u64 pidx, u64 axis) const;
  inline prec&       operator[](u64 pidx);
};

template<typename prec>
u64 regular_pointset<prec>::size() const
{
  return points;
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

} // namespace dptk
