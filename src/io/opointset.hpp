#pragma once

#include "../math/pointset.hpp"
#include <iosfwd>
#include <string>

namespace dptk {

void ostream_init(const std::string& in, std::ostream*& out);
void ostream_close(std::ostream*& out);

void write_label(std::ostream* os, const std::string& label);
void write_vector(std::ostream* os, const std::vector<b64>& list, i8 del);
void write_vector_param(std::ostream*           os,
                        const std::string&      label,
                        const std::vector<b64>& list,
                        i8                      del);

void write_pointset_header(std::ostream* os, const regular_pointset<b64>& pts, i8 del);
void write_pointset(std::ostream* os, const regular_pointset<b64>& pts, i8 del);

void write_pointset_footer(std::ostream* os, const regular_pointset<b64>& pts);
void write_pointset_eos(std::ostream* os);

} // namespace dptk
