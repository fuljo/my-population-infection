#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "utils.h"

/**
 * @brief Status of an individual
 *
 */
typedef enum individual_status {
  NOT_EXPOSED,
  EXPOSED,
  INFECTED,
  IMMUNE
} individual_status_t;

/**
 * @brief Returns a string representation of the given status code
 *
 * @param status
 * @return char*
 */
char *individual_status_string(int status);

typedef SLIST_HEAD(individual_list, individual) individual_list_t;

/**
 * @brief Represents an individual
 *
 */
typedef struct individual {
  unsigned long id; /**< Unique id for this individual in the world */
  double pos[2];    /**< (x,y) position */
  double displ[2]; /**< (dx, dy) displacement vector (velocity * t_step) applied
                      at each step */
  individual_status_t status; /**< current status of the individual */
  unsigned long t_status;     /**< Time passed since the individual entered the
                                current status */
  SLIST_ENTRY(individual) individuals;
} individual_t;

/**
 * @brief Print a human-readable representation of an individual
 *
 * @param ind individual
 */
void print_individual(individual_t *ind);

/**
 * @brief Creates an individual
 *
 * The individual will have the specified id, <tt>pos, displ = {0,0}</tt>
 * <tt>status = NOT_EXPOSED</tt>, <tt>t_status = 0</tt>.
 * The list pointer is uninitialized.
 *
 * @param id
 * @return individual_t
 */
individual_t *create_individual(unsigned long id);

/**
 * @brief Create an empty list of individuals
 *
 * @return individual_list_t
 */
individual_list_t create_individual_list();

#define insert_individual(head, ind) SLIST_INSERT_HEAD(head, ind, individuals)