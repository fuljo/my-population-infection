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

void print_individual(individual_t *ind);

individual_t *create_individual(unsigned long id);

individual_list_t create_individual_list();

#define insert_individual(head, ind) SLIST_INSERT_HEAD(head, ind, individuals)