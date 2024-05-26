#pragma once

template<typename T>
concept boundary = requires(T t) {
  { T() };
};
