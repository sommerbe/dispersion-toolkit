#include "opointset.hpp"
#include "ostream.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

namespace dptk {

void ostream_init(const std::string& in, std::ostream*& os)
{
  if (in == "-") {
    os = &std::cout;
  } else {
    std::ofstream* fs = new std::ofstream(in);
    if (fs->is_open()) {
      os = fs;
    } else {
      os = nullptr;
      delete fs;
    }
  }
}

void ostream_close(std::ostream*& os)
{
  if (os == nullptr)
    return;

  if (os == &std::cout) {
    os = nullptr;
    return;
  }

  std::ofstream* fs = (std::ofstream*)os;
  fs->close();
  delete fs;
  os = nullptr;
}

void write_label(std::ostream* os, const std::string& label)
{
  *os << "#" << label << " ";
}

void write_vector(std::ostream* os, const std::vector<b64>& list, i8 del)
{
  *os << list[0];
  for (u64 i = 1; i < list.size(); ++i) {
    *os << del << list[i];
  }
}

void write_vector_param(std::ostream*           os,
                        const std::string&      label,
                        const std::vector<b64>& list,
                        i8                      del)
{
  write_label(os, label);
  write_vector(os, list, del);
  *os << std::endl;
}

void write_pointset(std::ostream* os, const regular_pointset<b64>& pts, i8 del)
{
  assert(pts.domain_bound.size() > 1);

  ensure_precision(os, pts.domain_bound[0]);

  // write pointset domain
  write_vector_param(os, "domain", pts.domain_bound, del);

  // write pointset argument list
  if (!pts.arguments.empty()) {
    write_vector_param(os, "arg", pts.arguments, del);
  }

  // write coordinates
  if (pts.coords.empty()) {
    return;
  }
  ensure_precision(os, pts.coords[0]);
  for (u64 i = 0, d = 1; i < pts.coords.size(); ++i, ++d) {
    *os << pts.coords[i];
    if (d == pts.dimensions) {
      *os << std::endl;
      d = 0;
    } else {
      *os << del;
    }
  }
  write_pointset_eos(os);
}

void write_pointset_eos(std::ostream* os)
{
  *os << "#eos" << std::endl;
}

} // namespace dptk
