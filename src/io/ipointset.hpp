#pragma once

#include "../math/pointset.hpp"
#include <iosfwd>
#include <string>

namespace dptk {

void istream_init(const std::string& in, std::istream*& is);
void istream_close(std::istream*& is);

void read_pointset(std::istream& in, regular_pointset<b64>& out);

} // namespace dptk
