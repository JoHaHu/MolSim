#pragma once

#include <variant>

template<class... Ts>
struct overloads : Ts... {
  using Ts::operator()...;
};
