#include "arithmetic.hpp"
#include <algorithm>

namespace dptk {
namespace math {

b64 round(b64 v)
{
  return std::round(v);
}

b64 min(b64 a, b64 b)
{
  return std::min(a, b);
}

b64 max(b64 a, b64 b)
{
  return std::max(a, b);
}

b64 abs(b64 a)
{
  return std::abs(a);
}

} // namespace math
} // namespace dptk
