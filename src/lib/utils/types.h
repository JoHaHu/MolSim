#pragma once
#include <experimental/simd>

namespace stdx {
using namespace std::experimental;
using namespace std::experimental::__proposed;
}// namespace stdx

using double_v = stdx::native_simd<double>;
using double_mask = stdx::native_simd_mask<double>;
using long_v = stdx::native_simd<long>;
using long_mask = stdx::native_simd_mask<long>;
using size_v = stdx::native_simd<size_t>;
