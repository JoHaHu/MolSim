#pragma once
#include <experimental/simd>

namespace stdx {
using namespace std::experimental;
using namespace std::experimental::__proposed;
}// namespace stdx

using double_v = stdx::native_simd<double>;
using double_mask = stdx::native_simd_mask<double>;
using int_mask = stdx::native_simd_mask<int>;
using int_v = stdx::native_simd<int>;
using size_v = stdx::native_simd<size_t>;
