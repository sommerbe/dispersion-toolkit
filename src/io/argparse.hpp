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

u1 isarg_canonical(const std::string& arg);
u1 isarg_short(const std::string& arg);
u1 isarg_long(const std::string& arg);

u1 argval(const std::vector<std::string>& args, u64 idx);

u1 err(const std::string& msg);

u1 retrieve(const std::vector<std::string>& args, u64& idx, b64& val);
u1 retrieve(const std::vector<std::string>& args, u64& idx, u64& val);

template<typename any>
u1 retrieve(const std::vector<std::string>& args, u64& idx, std::vector<any>& val, u64 count = 0);


template<typename any>
u1 retrieve(const std::vector<std::string>& args, u64& idx, std::vector<any>& val, u64 count)
{
  any v;
  val.clear();
  while ((count == 0 || val.size() < count) && retrieve(args, idx, v)) {
    val.push_back(v);    
  }
  return count == 0 || val.size() == count;
}

} // namespace argparse
} // namespace dptk
