#pragma once

#include <iterator>
#include <ranges>
#include <vector>

namespace container {

template<std::forward_iterator Iter>
class combination_iterator {

 public:
  combination_iterator() = default;
  explicit combination_iterator(Iter start, Iter end) : i1(start), i2(start), end(end) {
    // Sets i2 to the first combination if it's not already at the end.
    // This can happen if the iterator is empty.
    if (i2 != end) {
      i2++;
    }
  }

  using difference_type = std::ptrdiff_t;
  using value_type = std::tuple<typename std::iterator_traits<Iter>::value_type, typename std::iterator_traits<Iter>::value_type>;
  using reference = std::tuple<typename std::iterator_traits<Iter>::reference, typename std::iterator_traits<Iter>::reference>;

  auto operator*() const -> reference {
    return {*i1, *i2};
  }

  auto operator++() -> combination_iterator & {
    if (i2 != end) {
      ++i2;
    }
    if (i2 == end && i1 != end) {
      ++i1;
      if (i1 != end) {
        i2 = std::next(i1);
        if (i2 == end) {
          i1 = end;
        }
      } else {
        i2 = end;
      }
    }
    return *this;
  }

  auto operator++(int) -> combination_iterator {
    auto tmp = *this;
    ++*this;
    return tmp;
  };

  auto operator==(const combination_iterator &iter) const -> bool {
    return i1 == iter.i1 && i2 == iter.i2 && end == iter.end;
  }

 private:
  Iter i1{}, i2{}, end{};
};

static_assert(std::forward_iterator<combination_iterator<std::vector<Particle>::iterator>>);

template<std::ranges::forward_range Range>
  requires(std::ranges::view<Range>)
class combination_view : public std::ranges::view_interface<combination_view<Range>> {

 public:
  combination_view() = default;
  explicit combination_view(Range &&range) : range(std::move(range)) {}

  using iterator = std::ranges::iterator_t<Range>;

  auto begin() -> combination_iterator<iterator> {
    return combination_iterator(range.begin(), range.end());
  }

  auto end() -> combination_iterator<iterator> {
    return combination_iterator(range.end(), range.end());
  }

 private:
  Range range;
};

struct _combination : std::ranges::range_adaptor_closure<_combination> {
  template<std::ranges::range R>
  constexpr auto operator()(R &&range) const {
    return combination_view(std::forward<R>(range));
  }
};

inline constexpr _combination combination;

;
static_assert(std::ranges::forward_range<combination_view<std::ranges::ref_view<std::vector<Particle>>>>);

}// namespace container
