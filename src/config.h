#pragma once

#include <argp.h>
#include <argz.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "utils.h"

#define LOG_DEFAULT LOG_INFO

#define T_INFECTION_DEFAULT 10 * MINUTE
#define T_RECOVERY_DEFAULT 10 * DAY
#define T_IMMUNITY_DEFAULT 3 * 30 * DAY

/* Configuration parameters */
typedef struct {
  unsigned long num_individuals, inf_individuals;
  unsigned long world_w, world_l, country_w, country_l;
  double velocity;
  double spreading_distance;
  unsigned long t_infection, t_recovery, t_immunity;
  unsigned long t_step;   /**< Simulation step in seconds */
  unsigned long t_target; /**< Stop simulation after this timestamp */
  unsigned int rand_seed;
  int log_level;
  bool write_trace; /**< Write a file with details of each ind. at each step */
} global_config_t;

/* Argument parser structures */
struct arguments {
  global_config_t config;
  char *argz;
  size_t argz_len;
};

int parse_opt(int key, char *arg, struct argp_state *state);

void init_config_default(global_config_t *cfg);

void print_config(global_config_t *cfg);

int validate_config(global_config_t *cfg, int world_size);