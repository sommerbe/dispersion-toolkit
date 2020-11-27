#pragma once

#include "../math/types.hpp"
#include <cassert>
#include <iomanip>
#include <ostream>

namespace dptk {

template<typename d>
void putln(std::ostream* os, d& data, u1 predicate = true);

template<typename d>
void putlnsci(std::ostream* os, d& data, u32 precision, u1 predicate = true);

void ensure_precision(std::ostream* os, const b64& data);

//
// implementation
//

void ensure_precision(std::ostream* os, const b64& data)
{
  *os << std::scientific << std::setprecision(16);
}

template<typename d>
void putln(std::ostream* os, d& data, u1 predicate)
{
  assert(os != nullptr);

  if (predicate) {
    *os << data << std::endl;
  }
}

template<typename d>
void putlnsci(std::ostream* os, d& data, u32 precision, u1 predicate)
{
  assert(os != nullptr);

  if (predicate) {
    *os << std::scientific << std::setprecision(precision) << data << std::endl;
  }
}

} // namespace dptk
