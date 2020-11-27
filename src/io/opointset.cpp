#include "opointset.hpp"
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

void write_pointset(std::ostream& os, regular_pointset<b64>& pts, i8 del)
{
  os << std::scientific << std::setprecision(16);
  for (u64 i = 0, d = 1; i < pts.coords.size(); ++i, ++d) {
    os << pts.coords[i];
    if (d == pts.dimensions) {
      os << std::endl;
      d = 1;
    } else {
      os << del;
    }
  }
  os << "#eos" << std::endl;
}

} // namespace dptk
