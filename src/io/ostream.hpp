#pragma once

#include "../math/types.hpp"
#include <cassert>
#include <iomanip>
#include <ostream>
#include <string>

namespace dptk {

template<typename d>
void put(std::ostream* os, d& data, u1 predicate = true);

template<typename d>
void putln(std::ostream* os, d& data, u1 predicate = true);

template<typename d>
void putlnsci(std::ostream* os, d& data, u32 precision, u1 predicate = true);

template<typename d>
void putsci(std::ostream* os, d& data, u32 precision, u1 predicate = true);

void ensure_precision(std::ostream* os, const b64& data);

template<typename d>
void put_header_column(std::ostream* os,
                       const d&      label,
                       i8&           delimiter_current,
                       i8            delimiter_next,
                       u1            predicate_column = true);

std::string extract_range(const std::string& s, const std::string& begin, const std::string& end);

//
// implementation
//

template<typename d>
void put(std::ostream* os, d& data, u1 predicate)
{
  assert(os != nullptr);

  if (predicate) {
    *os << data;
  }
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

template<typename d>
void putsci(std::ostream* os, d& data, u32 precision, u1 predicate)
{
  assert(os != nullptr);

  if (predicate) {
    *os << std::scientific << std::setprecision(precision) << data;
  }
}

template<typename d>
void put_header_column(std::ostream* os,
                       const d&      label,
                       i8&           delimiter_current,
                       i8            delimiter_next,
                       u1            predicate_column)
{
  if (!predicate_column)
    return;
  *os << delimiter_current << label;
  delimiter_current = delimiter_next;
}

} // namespace dptk
