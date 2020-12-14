#include "ostream.hpp"

namespace dptk {
void ensure_precision(std::ostream* os, const b64& data)
{
  *os << std::scientific << std::setprecision(16);
}

std::string extract_range(const std::string& s,
                          const std::string& begin,
                          const std::string& end)
{
  u64 b;
  u64 e;

  b = s.find(begin);
  if (b == std::string::npos)
    b = 0;
  e = s.find(end, b);

  return s.substr(b, e - b);
}

} // namespace dptk
