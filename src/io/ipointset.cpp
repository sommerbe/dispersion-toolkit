#include "ipointset.hpp"
#include <cctype>
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
      delete fs;
      throw std::invalid_argument("a non-existing or non-readable file can not be read "
                                  "from, and the following file was defined: "
                                  + in);
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

void read_vector(i8* begin, i8* end, std::vector<b64>& list)
{
  i8* e;
  b64 c;

  while (begin != end) {
    // delimiter between numbers
    if (std::isspace(*begin) != 0) {
      ++begin;
      continue;
    }
    // next number
    c = std::strtod(begin, &e);
    list.push_back(c);
    begin = e;
  }
  list.shrink_to_fit();
}

void read_vector(const std::string& in, u64 offset, std::vector<b64>& list)
{
  i8* b;
  i8* y;

  b = (i8*)in.data() + offset;
  y = (i8*)&in[in.size()];

  read_vector(b, y, list);
}

// static __attribute__((always_inline))
u1 starts_with(const std::string& str, const std::string& prefix)
{
  if (str.size() <= prefix.size() || str[0] != prefix[0]) {
    return false;
  }

  u1 eq = true;
  for (u64 i = 1; eq && i < prefix.size(); ++i) {
    eq &= str[i] == prefix[i];
  }
  return eq;
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
  u1          eos;

  eos = false;

  while (std::getline(in, ln) && !in.fail()) {
    if (ln.empty()) {
      continue;
    }

    // end of point set: #eos
    if (ln == "#eos") {
      eos = true;
      break;
    }

    // point set domain: #d low_0 low_1 ... low_(d-1) up_0 ... up_(d-1)
    // if (ln.size() > 3 && ln[0] == '#' && ln[1] == 'd' && ln[2] == ' ') {
    if (ln.size() > 3 && ln[0] == '#') {
      if (starts_with(ln, "#d ")) {
        read_vector(ln, 3, out.domain_bound);
        continue;
      } else if (starts_with(ln, "#arg ")) {
        read_vector(ln, 5, out.arguments);
      }
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

  // implicitly assume matrix as input
  if (!eos && out.domain_bound.empty()) {
    out.reset_inf_bound();
  }
}

void forward_delimiter(u1 predicate, const ipointset_read_info& inf, i8& delimiter)
{
  if (!predicate)
    return;
  delimiter = inf.delimiter;
}

} // namespace dptk
