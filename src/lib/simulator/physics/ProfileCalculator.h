#pragma once

#include "container/Container.h"
#include <array>
#include <fstream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

/**
 * @brief Class for calculating density and velocity profiles along the x-axis.
 */
template<const size_t DIMENSIONS>
class ProfileCalculator {
 public:
  /**
   * @brief Constructor for the ProfileCalculator class.
   * @param num_bins Number of bins to divide the x-axis.
   */
  ProfileCalculator(size_t num_bins)
      : num_bins(num_bins), x_min(std::numeric_limits<double>::max()), x_max(std::numeric_limits<double>::lowest()) {
    density_profile.resize(num_bins, 0);
    velocity_profile.resize(num_bins, std::array<double, DIMENSIONS>{0.0});
    counts.resize(num_bins, 0);
  }

  /**
   * @brief Updates the profiles based on the current particle data.
   * @param particles Container of particles.
   */
  void updateProfiles(container::Container<DIMENSIONS> &particles) {
    // Reset profiles
    std::fill(density_profile.begin(), density_profile.end(), 0);
    std::fill(velocity_profile.begin(), velocity_profile.end(), std::array<double, DIMENSIONS>{0.0});
    std::fill(counts.begin(), counts.end(), 0);

    // Calculate the min and max x values
    x_min = std::numeric_limits<double>::max();
    x_max = std::numeric_limits<double>::lowest();
    particles.linear([this](Particles<DIMENSIONS> &particles, size_t index) {
      double x = particles.positions[0][index];
      if (x < x_min) x_min = x;
      if (x > x_max) x_max = x;
    });

    // Calculate the bin size based on min and max x values
    bin_size = (x_max - x_min) / num_bins;

    // Calculate the profiles
    particles.linear([this](Particles<DIMENSIONS> &particles, size_t index) {
      double x = particles.positions[0][index];
      size_t bin = static_cast<size_t>((x - x_min) / bin_size);
      // Ensure the bin is within valid range
      if (bin < num_bins) {
        density_profile[bin]++;
        for (size_t d = 0; d < DIMENSIONS; ++d) {
          velocity_profile[bin][d] += particles.velocities[d][index];
        }
        counts[bin]++;
      }
    });

    // Calculate averages
    for (size_t bin = 0; bin < num_bins; ++bin) {
      if (counts[bin] > 0) {
        for (size_t d = 0; d < DIMENSIONS; ++d) {
          velocity_profile[bin][d] /= counts[bin];
        }
      }
    }
  }

  /**
   * @brief Writes the profiles to a CSV file.
   * @param filename The name of the CSV file.
   */
  void writeProfilesToCSV(const std::string &filename) const {
    std::ofstream file(filename);
    file << "Bin,Density";
    for (size_t d = 0; d < DIMENSIONS; ++d) {
      file << ",Velocity" << d;
    }
    file << "\n";

    for (size_t bin = 0; bin < num_bins; ++bin) {
      file << bin << "," << density_profile[bin];
      for (size_t d = 0; d < DIMENSIONS; ++d) {
        file << "," << velocity_profile[bin][d];
      }
      file << "\n";
    }
    file.close();
  }

  /**
   * @brief Returns the density profile.
   * @return Reference to the density profile vector.
   */
  const std::vector<int> &getDensityProfile() const {
    return density_profile;
  }

  /**
   * @brief Returns the velocity profile.
   * @return Reference to the velocity profile vector.
   */
  const std::vector<std::array<double, DIMENSIONS>> &getVelocityProfile() const {
    return velocity_profile;
  }

  /**
   * @brief Returns the minimum x-value.
   * @return Minimum x-value.
   */
  double getXMin() const {
    return x_min;
  }

  /**
   * @brief Returns the maximum x-value.
   * @return Maximum x-value.
   */
  double getXMax() const {
    return x_max;
  }

  /**
   * @brief Returns the bin size.
   * @return Bin size.
   */
  double getBinSize() const {
    return bin_size;
  }

 private:
  size_t num_bins;                                             ///< Number of bins to divide the x-axis
  double x_min;                                                ///< Minimum x value
  double x_max;                                                ///< Maximum x value
  double bin_size;                                             ///< Size of each bin
  std::vector<int> density_profile;                            ///< Density profile
  std::vector<std::array<double, DIMENSIONS>> velocity_profile;///< Velocity profile
  std::vector<size_t> counts;                                  ///< Count of particles in each bin
};
