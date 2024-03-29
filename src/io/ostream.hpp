#pragma once

#include "../math/types.hpp"
#include <cassert>
#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

namespace dptk {

template<typename d>
void put(std::ostream* os, d& data, u1 predicate = true);

template<typename d>
void putln(std::ostream* os, d& data, u1 predicate = true);

template<typename d0, typename d1>
void putln(std::ostream* os, d0& key, d1& value, u1 predicate = true);

template<typename d>
void putlnsci(std::ostream* os, d& data, u32 precision, u1 predicate = true);

template<typename d>
void putsci(std::ostream* os, d& data, u32 precision, u1 predicate = true);

template<typename d>
void putparam(std::ostream*      os,
              const std::string& name,
              const d&           data,
              u1                 predicate = true);

template<typename d>
void putparam(std::ostream*         os,
              const std::string&    name,
              const std::vector<d>& data,
              u1                    predicate = true,
              u8                    delimiter = ' ');

void ensure_precision(std::ostream* os, const b64& data);

template<typename d>
void put_header_column(std::ostream* os,
                       const d&      label,
                       i8&           delimiter_current,
                       i8            delimiter_next,
                       u1            predicate_column = true);

std::string extract_range(const std::string& s,
                          const std::string& begin,
                          const std::string& end);

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

template<typename d0, typename d1>
void putln(std::ostream* os, d0& key, d1& value, u1 predicate)
{
  assert(os != nullptr);

  if (predicate) {
    *os << key << value << std::endl;
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
void putparam(std::ostream* os, const std::string& name, const d& data, u1 predicate)
{
  assert(os != nullptr);

  if (predicate) {
    *os << "# " << name << " = (" << data << ")" << std::endl;
  }
}

template<typename d>
void putparam(std::ostream*         os,
              const std::string&    name,
              const std::vector<d>& data,
              u1                    predicate,
              u8                    delimiter)

{
  assert(os != nullptr);

  if (predicate) {
    *os << "# " << name << " = (";
    for (u64 i = 0; i < data.size(); ++i) {
      if (i > 0) {
        *os << delimiter;
      }
      *os << data[i];
    }
    *os << ")" << std::endl;
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
