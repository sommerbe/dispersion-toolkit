#include "argparse.hpp"
#include <iostream>
#include <stdexcept>
#include <cctype>

namespace dptk {
namespace argparse {

void replace(std::string& str, const std::string& from, const std::string& to)
{
  std::size_t p = 0;
  while (p < str.size() && (p = str.find(from, p)) != str.npos) {
    str.replace(p, from.size(), to);
    p += to.size();
  }
}

void canonical(dptk::i32 argc, const dptk::i8** argv, std::vector<std::string>& args)
{
  args.clear();

  for (i32 i = 1; i < argc; ++i) {
    std::string s = argv[i];

    replace(s, "\\t", "\t");
    replace(s, "\\n", "\n");
    replace(s, "\\f", "\f");
    replace(s, "\\r", "\r");
    // replace(s, "\\'", "\'");
    // replace(s, "\\?", "\?");

    // skip empty elements
    if (s.empty()) {
      continue;
    }

    // expand short options
    // pattern: -abcd... => -a -b -c -d ...
    if (isarg_short(s)) {
      for (u64 j = 1; j < s.size(); ++j) {
        args.push_back("-" + s.substr(j, 1));
      }
      continue;
    }

    // expand long options
    if (isarg_long(s)) {
      u64 p;
      // pattern: --key=value => --key value
      if (s.size() > 4 && (p = s.find_first_of('=', 3)) != std::string::npos) {
        args.push_back("--" + s.substr(2, p - 2));
        args.push_back(s.substr(p + 1));
      } else {
        // pattern: --key value | --key
        // if value does not start with -; it is emitted as normal value below
        args.push_back(s);
      }
      continue;
    }

    // emit normal argument
    args.push_back(s);
  }
}

void require_argval(const std::vector<std::string>& args,
                    u64                             idx,
                    const std::string&              error_message)
{
  if (idx + 1 == args.size()) {
    throw std::invalid_argument(error_message.c_str());
  }
}

void ensure(u1 predicate, const std::string& error_message)
{
  if (!predicate) {
    throw std::invalid_argument(error_message.c_str());
  }
}

u1 isarg_short(const std::string& s)
{
  return s.size() > 1 && s[0] == '-' && s[1] != '-' && !std::isdigit((u8)s[1]) && s[1] != '.';
}

u1 isarg_long(const std::string& s)
{
  return s.size() > 2 && s[0] == '-' && s[1] == '-' && s[2] != '-';
}

u1 isarg_canonical(const std::string& arg)
{
  if (arg.size() == 2 && isarg_short(arg))
    return true;
  if (isarg_long(arg))
    return true;
  return false;
}

u1 argval(const std::vector<std::string>& args, u64 idx)
{
  return idx + 1 < args.size() && !isarg_canonical(args[idx+1]);
}

u1 err(const std::string& msg)
{
  std::cerr << msg << std::endl;
  return false;
}

u1 retrieve(const std::vector<std::string>& args, u64& idx, b64& val)
{
  if (!argparse::argval(args, idx))
    return false;
  
  val = std::strtod(args[++idx].c_str(), nullptr);
  return true;
}

u1 retrieve(const std::vector<std::string>& args, u64& idx, u64& val)
{
  if (!argparse::argval(args, idx))
    return false;
  
  val = std::strtoul(args[++idx].c_str(), nullptr, 10);
  return true;
}

} // namespace argparse
} // namespace dptk