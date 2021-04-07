#include "config.h"

/* Argument parser variables */
const char *argp_program_bug_address = "alessandro.fulgini@mail.polimi.it";
const char *argp_program_version = "0.1-alpha";

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
int parse_opt(int key, char *arg, struct argp_state *state) {
  struct arguments *a = state->input;
  struct global_config *cfg = &a->config;
  /* Remove leading equal sign from the string */
  if (arg && arg[0] == '=') {
    arg++;
  }
  switch (key) {
    case 'N': {
      cfg->num_individuals = strtoul(arg, NULL, 10);
      break;
    }
    case 'I': {
      cfg->inf_individuals = strtoul(arg, NULL, 10);
      break;
    }
    case 'W': {
      cfg->world_w = strtoul(arg, NULL, 10);
      break;
    }
    case 'L': {
      cfg->world_l = strtoul(arg, NULL, 10);
      break;
    }
    case 'w': {
      cfg->country_w = strtoul(arg, NULL, 10);
      break;
    }
    case 'l': {
      cfg->country_l = strtoul(arg, NULL, 10);
      break;
    }
    case 'v': {
      cfg->velocity = strtod(arg, NULL);
      break;
    }
    case 'd': {
      cfg->spreading_distance = strtod(arg, NULL);
      break;
    }
    case 333: {
      cfg->t_infection = arg ? strtoul(arg, NULL, 10) : T_INFECTION_DEFAULT;
      break;
    }
    case 444: {
      cfg->t_recovery = arg ? strtoul(arg, NULL, 10) : T_RECOVERY_DEFAULT;
      break;
    }
    case 555: {
      cfg->t_immunity = arg ? strtoul(arg, NULL, 10) : T_IMMUNITY_DEFAULT;
      break;
    }
    case 666: {
      cfg->t_step = strtoul(arg, NULL, 10);
      break;
    }
    case 777: {
      cfg->t_target = strtoul(arg, NULL, 10) * DAY;
      break;
    }
    case 888: {
      cfg->rand_seed = arg ? strtoul(arg, NULL, 10) : time(NULL);
      break;
    }
    case ARGP_KEY_INIT: {
      a->argz = 0;
      a->argz_len = 0;
      break;
    }
    case ARGP_KEY_ARG: {
      argz_add(&a->argz, &a->argz_len, arg);
      break;
    }
    case ARGP_KEY_END: {
      /* This happens after all options and arguments */
      break;
    }
  }
  return 0;
}

/**
 * @brief Validates configuration and prints error to stderr.
 *
 * @param cfg global configuration
 * @param world_size number of MPI processes in world
 * @return int status (0: ok, 1: error)
 */
int validate_config(struct global_config *cfg, int world_size) {
  /* Infected */
  if (cfg->inf_individuals > cfg->num_individuals) {
    fprintf(stderr,
            "[ERR ] Infected individuals exceed population size (%lu > %lu)\n",
            cfg->inf_individuals, cfg->num_individuals);
    return 1;
  }
  /* World dimensions */
  if (cfg->world_w < cfg->country_w || cfg->world_l < cfg->country_l) {
    fprintf(stderr,
            "[ERR ] Country dimensions cannot exceed world dimensions\n");
    return 1;
  }
  if (cfg->world_w % cfg->country_w != 0) {
    fprintf(stderr, "[ERR ] Country width must divide world width\n");
    return 1;
  }
  if (cfg->world_l % cfg->country_l != 0) {
    fprintf(stderr, "[ERR ] Country length must divide world length\n");
    return 1;
  }
  int num_countries =
      (cfg->world_w / cfg->country_w) * (cfg->world_l / cfg->country_l);
  if (world_size != num_countries) {
    fprintf(stderr,
            "[ERR ] Number of processes does not match number of countries. "
            "Expected %d, got %d\n",
            num_countries, world_size);
    return 1;
  }
  /* Velocity */
  if (cfg->velocity <= 0.) {
    fprintf(stderr, "[ERR ] Velocity must be non-negative\n");
    return 1;
  }
  /* Spreading distance */
  if (cfg->spreading_distance <= 0.) {
    fprintf(stderr, "[ERR ] Spreading distance must be non-negative\n");
    return 1;
  }
  /* Simulation step */
  if (cfg->t_step > DAY) {
    fprintf(stderr, "[ERR ] Simulation step cannot be longer than one day\n");
    return 1;
  }
  /* Relation between movement and time step */
  if (cfg->t_step * cfg->velocity > MIN(cfg->country_w, cfg->country_l)) {
    fprintf(stderr,
            "[ERR ] The movement at each step is larger than a country: t * v "
            "= %f > %lu\n",
            cfg->t_step * cfg->velocity, MIN(cfg->country_w, cfg->country_l));
    return 1;
  }

  /* If we got here the configuration is valid */
  return 0;
}

/**
 * @brief Prints the configuration.
 *
 * @param cfg configuration
 */
void print_config(struct global_config *cfg) {
  printf(
      "--------------------\nGlobal configuration\n--------------------\n "
      "num_individuals %lu\n "
      "inf_individuals %lu\n world_w %lu\n world_l %lu\n country_w %lu\n "
      "country_l %lu\n velocity %f\n spreading_distance %f\n t_infection "
      "%lu\n t_recovery %lu\n t_immunity %lu\n t_step %lu\n t_target "
      "%lu\n--------------------\n",
      cfg->num_individuals, cfg->inf_individuals, cfg->world_w, cfg->world_l,
      cfg->country_w, cfg->country_l, cfg->velocity, cfg->spreading_distance,
      cfg->t_infection, cfg->t_recovery, cfg->t_immunity, cfg->t_step,
      cfg->t_target);
}