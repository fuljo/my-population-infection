#pragma once

#include <argp.h>
#include <argz.h>
#include <stdio.h>
#include <stdlib.h>

#include "utils.h"

#define T_INFECTION_DEFAULT 10 * MINUTE
#define T_RECOVERY_DEFAULT 10 * DAY
#define T_IMMUNITY_DEFAULT 3 * 30 * DAY

/* Configuration parameters */
struct global_config {
  unsigned long num_individuals, inf_individuals;
  unsigned long world_w, world_l, country_w, country_l;
  double velocity;
  double spreading_distance;
  unsigned long t_infection, t_recovery, t_immunity;
  unsigned long t_step;      /* Simulation step in seconds */
  unsigned long t_target; /* Stop simulation after this timestamp */
};

/* Argument parser structures */
struct arguments {
  struct global_config config;
  char *argz;
  size_t argz_len;
};

/**
 * @brief Handler for argp options and arguments.
 *
 * It fills in the global configuration in state->input->config with the
 * provided options and validates them at the end.
 * The arguments are instead encoded in the argz vector.
 *
 * @param key argp key of the current option / argument
 * @param arg string value of the option / argument
 * @param state argp state
 * @return int status
 */
int parse_opt(int key, char *arg, struct argp_state *state);

/**
 * @brief Prints the configuration.
 * 
 * @param cfg configuration
 */
void print_config(struct global_config *cfg);

/**
 * @brief Validates the arguments and configuration in state->input.
 * 
 * @param state argp state at the end of parsing (ARGP_KEY_END)
 * @return int status
 */
int validate_config(struct global_config *cfg, int world_size);