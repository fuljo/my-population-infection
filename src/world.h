#pragma once

/**
 * @brief Cardinal points
 *
 */
typedef enum cardinal_point {
  NORTH,
  NORTH_EAST,
  EAST,
  SOUTH_EAST,
  SOUTH,
  SOUTH_WEST,
  WEST,
  NORTH_WEST,
} cardinal_point_t;

/**
 * @brief Limits of a country
 *
 * The lower bound is considered inclusive, the upper exclusive
 *
 */
typedef struct limits {
  unsigned long xmin, xmax;
  unsigned long ymin, ymax;
} limits_t;

#define NEIGHBOR_COUNT 8
#define NEIGHBOR_OPPOSITE_DISTANCE 4

void distribute_population_uniform(unsigned long population,
                                   unsigned int num_countries,
                                   unsigned long res[]);

int decode_cardinal_point_flag(unsigned char flag);
