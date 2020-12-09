#include "ostream.hpp"

namespace dptk {
void ensure_precision(std::ostream* os, const b64& data)
{
  *os << std::scientific << std::setprecision(16);
}

} // namespace dptk
