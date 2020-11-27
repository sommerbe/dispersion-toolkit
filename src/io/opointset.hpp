#pragma once

#include "../math/pointset.hpp"
#include <iosfwd>
#include <string>

namespace dptk {

void ostream_init(const std::string& in, std::ostream*& out);
void ostream_close(std::ostream*& out);

void write_pointset(std::istream& os, regular_pointset<b64>& pts, i8 del);

} // namespace dptk
