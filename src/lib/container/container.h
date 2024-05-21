#pragma once

#include "Particle.h"
#include <ranges>

template<typename T>
concept ParticleContainer =
    std::ranges::forward_range<T>
    && std::same_as<typename std::ranges::range_value_t<T>, Particle>;
