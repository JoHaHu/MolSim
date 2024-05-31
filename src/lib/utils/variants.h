#pragma once

#include <variant>

template<class... Ts>
struct overloaded : Ts... {
  using Ts::operator()...;
};
