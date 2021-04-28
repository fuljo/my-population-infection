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

/**
 * @brief Represents the summary of the number of individuals for each status
 *
 */
typedef struct summary {
  unsigned long subsceptible;
  unsigned long infected;
  unsigned long immune;
} summary_t;

void print_individual(individual_t *ind);

individual_t *create_individual(unsigned long id);

individual_list_t create_individual_list();

#define INDIVIDUAL_DISTANCE(ind1, ind2)      \
  sqrt(pow(ind1->pos[0] - ind2->pos[0], 2) + \
       pow(ind1->pos[1] - ind2->pos[1], 2))

#define INDIVIDUAL_INSERT(head, ind) SLIST_INSERT_HEAD(head, ind, individuals)

#define INDIVIDUAL_EMPTY(head) SLIST_EMPTY(head)

#define INDIVIDUAL_FIRST(head) SLIST_FIRST(head)

#define INDIVIDUAL_NEXT(ind) SLIST_NEXT(ind, individuals)

#define INDIVIDUAL_FOREACH(var, head) SLIST_FOREACH(var, head, individuals)

#define INDIVIDUAL_REMOVE_HEAD(head) SLIST_REMOVE_HEAD(head, individuals)

#define INDIVIDUAL_REMOVE_AFTER(ind) SLIST_REMOVE_AFTER(ind, individuals)

#define INDIVIDUAL_COUNT(head, count) \
  SLIST_COUNT(head, individuals, individual, count)

#define INDIVIDUAL_REMOVE(head, ind) \
  SLIST_REMOVE(head, ind, individual, individuals)
