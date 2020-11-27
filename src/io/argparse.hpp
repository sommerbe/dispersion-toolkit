#pragma once

#include "../math/types.hpp"
#include <string>
#include <vector>

namespace dptk {
namespace argparse {

void canonical(dptk::i32 argc, const dptk::i8** argv, std::vector<std::string>& args);

} // namespace argparse
} // namespace dptk
