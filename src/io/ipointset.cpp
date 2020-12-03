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

void read_pointset(std::istream& in, regular_pointset<b64>& out, ipointset_read_info* inf)
{
  std::string ln;
  i8*         b;
  i8*         e;
  i8*         y;
  i8*         z;
  u64         d;
  b64         c;

  while (std::getline(in, ln) && !in.fail()) {
    if (ln.empty()) {
      continue;
    }

    // end of point set: #eos
    if (ln == "#eos") {
      break;
    }

    // point set domain: #d low_0 low_1 ... low_(d-1) up_0 ... up_(d-1)
    if (ln.size() > 3 && ln[0] == '#' && ln[1] == 'd' && ln[2] == ' ') {
      b = (i8*)ln.data() + 3;
      y = &ln[ln.size()];
      while (b != y) {
        if (std::isspace(*b) != 0) {
          ++b;
          continue;
        }
        c = std::strtod(b, &e);
        out.domain_bound.push_back(c);
        b = e;
      }
      out.domain_bound.shrink_to_fit();
      continue;
    }

    // comment: # ...
    if (ln[0] == '#') {
      continue;
    }

    // coordinates
    z = (i8*)ln.data();
    y = &ln[ln.size()];
    b = z;
    d = 0;
    while (b != y) {
      if (std::isspace(*b) != 0) {
        // remember delimiter
        if (out.dimensions == 0 && inf != nullptr) {
          inf->delimiter = *b;
        }
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

void forward_delimiter(u1 predicate, const ipointset_read_info& inf, i8& delimiter)
{
  if (!predicate)
    return;
  delimiter = inf.delimiter;
}

} // namespace dptk
