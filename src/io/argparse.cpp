#include "argparse.hpp"

namespace dptk {
namespace argparse {

void canonical(dptk::i32 argc, const dptk::i8** argv, std::vector<std::string>& args)
{
  args.clear();

  for (i32 i = 1; i < argc; ++i) {
    std::string s = argv[i];

    // skip empty elements
    if (s.empty()) {
      continue;
    }

    // expand short options
    // pattern: -abcd... => -a -b -c -d ...
    if (s[0] == '-' && s.size() > 1 && s[1] != '-') {
      for (u64 j = 1; j < s.size(); ++j) {
        args.push_back("-" + s.substr(j, 1));
      }
      continue;
    }

    // expand long options
    if (s.size() > 2 && s[0] == '-' && s[1] == '-' && s[2] != '-') {
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

} // namespace argparse
} // namespace dptk