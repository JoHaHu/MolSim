#pragma once

#include <memory>
#include <ranges>
#include <vector>

namespace container {

template<std::ranges::forward_range R>
//  requires(std::ranges::view<R>)
class shared_view : public std::ranges::view_interface<shared_view<R>> {
 public:
  shared_view() = default;
  explicit shared_view(std::shared_ptr<R> &&range) : range(std::move(range)) {}


  using iterator = std::ranges::iterator_t<R>;
  using sentinel = std::ranges::sentinel_t<R>;

  auto begin() -> iterator {
    return range.get()->begin();
  }
  auto end() -> sentinel {
    return range.get()->end();
  }

 private:
  std::shared_ptr<R> range;
};

struct _shared : std::ranges::range_adaptor_closure<_shared> {
  template<std::ranges::range R>
    requires(std::ranges::view<R>)
  constexpr auto operator()(R &&range) const {
    return shared_view(std::forward<R>(range));
  }
  template<std::ranges::range R>
  constexpr auto operator()(R &range) const {
    return shared_view(std::ranges::ref_view<R>(range));
  }
};

inline constexpr _shared shared;

static_assert(std::ranges::forward_range<shared_view<std::ranges::ref_view<std::vector<int>>>>);

}// namespace container
