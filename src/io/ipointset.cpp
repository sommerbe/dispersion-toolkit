#include "ipointset.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace dptk {

void istream_init(const std::string& in, std::istream*& is)
{
  if (in == "-") {
    is = &std::cin;
  } else {
    std::ifstream* fs = new std::ifstream(in);
    if (fs->is_open()) {
      is = fs;
    } else {
      is = nullptr;
    }
  }
}

void istream_close(std::istream*& is)
{
  if (is == nullptr)
    return;

  if (is == &std::cin) {
    is = nullptr;
    return;
  }

  std::ifstream* fs = (std::ifstream*)is;
  fs->close();
  delete fs;
  is = nullptr;
}

void read_pointset(std::istream& in, regular_pointset<b64>& out)
{
  std::string ln;
  i8*         b;
  i8*         e;
  i8*         y;
  i8*         z;
  u64         d;
  b64         c;

  while (std::getline(in, ln) && !in.fail()) {
    if (ln.empty() || ln[0] == '#') {
      continue;
    }
    z = (i8*)ln.data();
    y = &ln[ln.size()];
    b = z;
    d = 0;
    while (b != y) {
      if (std::isspace(*b) != 0) {
        ++b;
        continue;
      }
      c = std::strtod(b, &e);
      out.coords.push_back(c);
      ++d;
      b = e;
    }
    // add point
    if (d > 0) {
      ++out.points;
    }

    // check: dimension match
    if (out.dimensions == 0) {
      out.dimensions = d;
    } else if (out.dimensions != d) {
      std::string m = "pointset contains points of variable dimensions violating the "
                      "regular pointset condition. ";
      m += "expected dimension = " + std::to_string(out.dimensions);
      m += "read dimension = " + std::to_string(d);
      m += "line = " + ln;
      m += ".";
      throw std::runtime_error(m);
    }
  }
}

} // namespace dptk
