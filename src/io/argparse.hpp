#pragma once

#include "../math/types.hpp"
#include <string>
#include <vector>

namespace dptk {
namespace argparse {

void replace(std::string& str, const std::string& from, const std::string& to);

void canonical(dptk::i32 argc, const dptk::i8** argv, std::vector<std::string>& args);

void require_argval(const std::vector<std::string>& args,
                    u64                             idx,
                    const std::string&              error_message);

void ensure(u1 predicate, const std::string& error_message);

u1 argval(const std::vector<std::string>& args, u64 idx);

u1 err(const std::string& msg);

} // namespace argparse
} // namespace dptk
