#include "pointset_dsorted.hpp"

namespace dptk {
void pointset_dsorted_index::allocate(u64 n, u64 d)
{
  stride     = n;
  dimensions = d;
  indices.resize(n * d);
}

u64& pointset_dsorted_index::at(u64 i, u64 axis)
{
  assert(i >= 0 && i < stride);
  return indices[i + axis * stride];
}

u64 pointset_dsorted_index::search(u64 val, u64 axis)
{
  u64 begin = axis * stride;
  u64 end   = (axis + 1) * stride;
  for (u64 i = begin; i < end; ++i) {
    if (indices[i] == val)
      return i - begin;
  }
  return end;
}

std::vector<u64>::iterator pointset_dsorted_index::begin(u64 axis)
{
  return indices.begin() + axis * stride;
}

std::vector<u64>::iterator pointset_dsorted_index::end(u64 axis)
{
  return indices.begin() + (axis + 1) * stride;
}

} // namespace dptk
