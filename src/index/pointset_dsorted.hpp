#pragma once
#include "../math/pointset.hpp"
#include <algorithm>

namespace dptk {

struct pointset_dsorted_index
{
  u64              stride;
  u64              dimensions;
  std::vector<u64> indices;

  void allocate(u64 n, u64 d);

  u64& at(u64 i, u64 axis);
  u64  search(u64 val, u64 axis);

  std::vector<u64>::iterator begin(u64 axis);
  std::vector<u64>::iterator end(u64 axis);
};

template<typename prec>
void sort(const regular_pointset<prec>& pts, pointset_dsorted_index& idx, u64 axis);

template<typename prec>
void sort(const regular_pointset<prec>& pts, pointset_dsorted_index& idx);

//
// implementation
//

/**
 * Assigns indices of points stored in pts, and sorts them, both along given axis.
 */
template<typename prec>
void sort(const regular_pointset<prec>& pts, pointset_dsorted_index& idx, u64 axis)
{
  assert(axis < pts.dimensions);
  assert(axis < idx.dimensions);
  assert(pts.size() == idx.stride);

  for (u64 i = 0; i < pts.size(); ++i)
    idx.at(i, axis) = i;

  std::sort(idx.begin(axis), idx.end(axis), [&pts, axis](u64 a, u64 b) {
    return *pts.at(a, axis) < *pts.at(b, axis);
  });
}

/**
 * Assigns indices of all points, and sorts them, along all axes.
 */
template<typename prec>
void sort(const regular_pointset<prec>& pts, pointset_dsorted_index& idx)
{
  for (u64 axis = 0; axis < pts.dimensions; ++axis)
    sort(pts, idx, axis);
}

} // namespace dptk
