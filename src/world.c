#include "world.h"

/**
 * @brief Uniformly distributes a population between countries.
 *
 * All countries will get <tt>floor(population / num_countries)</tt>
 * individuals, and the first country will also get the remaining <tt>population
 * % num_countries</tt> individuals
 *
 * @param[in] population number of individuals to be distributed
 * @param[in] num_countries number of countries
 * @param[out] res array of size \p num_countries that will hold the result
 */
void distribute_population_uniform(unsigned long population,
                                   unsigned int num_countries,
                                   unsigned long res[]) {
  unsigned long frac = population / num_countries;
  for (int c = 0; c < num_countries; c++) {
    res[c] = frac;
  }
  res[0] += population % num_countries;
}

/**
 * @brief Decode a flag in which each bit represents the cardinal points.
 *
 * Assuming only one of NSWE are set and that bits are indexed according to
 * \c enum \c cardinal_point , returns the composite cardinal point (ex. if both
 * NORTH and WEST are set, it returns NORTH_WEST)
 *
 * @param[in] flag
 * @return the corresponding cardinal point, -1 if error
 */
int decode_cardinal_point_flag(unsigned char flag) {
  switch (flag) {
    case (1 << NORTH): {
      return NORTH;
    }
    case (1 << NORTH | 1 << EAST): {
      return NORTH_EAST;
    }
    case (1 << EAST): {
      return EAST;
    }
    case (1 << SOUTH | 1 << EAST): {
      return SOUTH_EAST;
    }
    case (1 << SOUTH): {
      return SOUTH;
    }
    case (1 << SOUTH | 1 << WEST): {
      return SOUTH_WEST;
    }
    case (1 << WEST): {
      return WEST;
    }
    case (1 << NORTH | 1 << WEST): {
      return NORTH_WEST;
    }
    default:
      return -1;
  }
}