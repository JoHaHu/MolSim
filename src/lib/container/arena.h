#pragma once

#include <algorithm>
#include <cstdlib>
#include <ranges>
#include <vector>

namespace container {

/**
 * A memory arena, to keep references valid when objects gets deleted
 * */
template<typename T>
class arena {

 public:
  struct entry {
    T data;
    bool active;
  };
  explicit arena(size_t size) : allocation(std::vector<entry>()) {
    allocation.reserve(size);
  }

  auto range_entries() -> auto {
    return allocation
        | std::views::filter(&entry::active);
  }

  auto emplace_back(Particle &&particle) -> auto & {
    return allocation.emplace_back(particle, true);
  }

  auto size() -> size_t {
    return std::ranges::fold_left(range_entries(), 0, [](auto acc, auto &c) { return ++acc; });
  }

 private:
  std::vector<entry> allocation;
};

}// namespace container