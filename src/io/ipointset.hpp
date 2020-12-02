#pragma once

#include "../math/pointset.hpp"
#include <iosfwd>
#include <string>

namespace dptk {

struct ipointset_read_info
{
  i8 delimiter;
};

void istream_init(const std::string& in, std::istream*& is);
void istream_close(std::istream*& is);

void read_pointset(std::istream&          in,
                   regular_pointset<b64>& out,
                   ipointset_read_info*   inf = nullptr);

void forward_delimiter(u1 predicate, const ipointset_read_info& inf, i8& delimiter);

} // namespace dptk
