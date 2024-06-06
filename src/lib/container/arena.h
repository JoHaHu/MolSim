#pragma once

#include <algorithm>
#include <cstdlib>
#include <ranges>
#include <vector>

namespace container {

  /**
   * @brief A memory arena to maintain valid references when objects are deleted.
   *
   * The arena class manages memory for objects, ensuring that references to active objects remain valid even if some objects are removed.
   *
   * @tparam T Type of the objects stored in the arena.
   */
  template<typename T>
  class arena {

  public:
    struct entry {
      T data;
      bool active;
    };

    /**
   * @brief Constructs an arena with a reserved size.
   *
   * @param size Initial reserve size for the arena.
   */
    explicit arena(size_t size) : allocation(std::vector<entry>()) {
      allocation.reserve(size);
    }

    /**
   * @brief Returns a range of active entries.
   *
   * @return A filtered range of active entries.
   */
    auto range_entries() -> auto {
      return allocation
          | std::views::filter(&entry::active);
    }

    /**
   * @brief Adds a new particle to the arena.
   *
   * @param particle The particle to be added.
   * @return Reference to the added entry.
   */
    auto emplace_back(Particle &&particle) -> auto & {
      return allocation.emplace_back(particle, true);
    }

    /**
   * @brief Returns the number of active entries in the arena.
   *
   * @return Number of active entries.
   */
    auto size() -> size_t {
      return std::ranges::fold_left(range_entries(), 0, [](auto acc, auto &c) { return ++acc; });
    }

  private:
    std::vector<entry> allocation;
  };

}// namespace container