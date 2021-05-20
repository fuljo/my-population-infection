#include "config.h"

/* Argument parser variables */
const char *argp_program_bug_address = "alessandro.fulgini@mail.polimi.it";
const char *argp_program_version = "0.1-alpha";

/**
 * @brief Decodes log level code from string
 *
 * @param[in] arg log level string, case-insensitive, not null
 * @return int log level code, -1 if unknown
 */
int decode_log_level(char *arg) {
  if (strncasecmp(arg, "TRACE", 5) == 0) {
    return LOG_TRACE;
  }
  if (strncasecmp(arg, "DEBUG", 5) == 0) {
    return LOG_DEBUG;
  }
  if (strncasecmp(arg, "INFO", 5) == 0) {
    return LOG_INFO;
  }
  if (strncasecmp(arg, "WARN", 5) == 0) {
    return LOG_WARN;
  }
  if (strncasecmp(arg, "ERROR", 5) == 0) {
    return LOG_ERROR;
  }
  if (strncasecmp(arg, "FATAL", 5) == 0) {
    return LOG_FATAL;
  }
  return -1;
}

/**
 * @brief Handler for argp options and arguments.
 *
 * It fills in the global configuration in state->input->config with the
 * provided options and validates them at the end.
 * The arguments are instead encoded in the argz vector.
 *
 * @param[in] key argp key of the current option / argument
 * @param[in] arg string value of the option / argument
 * @param[in,out] state argp state
 * @return int status
 */
int parse_opt(int key, char *arg, struct argp_state *state) {
  struct arguments *a = state->input;
  global_config_t *cfg = &a->config;
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
    case 999: {
      if (arg == NULL || strlen(arg) == 0) {
        cfg->log_level = LOG_DEFAULT;
      } else {
        cfg->log_level = decode_log_level(arg);
      }
      break;
    }
    case 101010: {
      cfg->write_trace = true;
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
 * @brief Initialize configuration with default values
 *
 * @param[out] cfg configuration to be initialized
 */
void init_config_default(global_config_t *cfg) {
  cfg->t_infection = T_INFECTION_DEFAULT;
  cfg->t_recovery = T_RECOVERY_DEFAULT;
  cfg->t_immunity = T_IMMUNITY_DEFAULT;
  cfg->rand_seed = time(NULL);
  cfg->log_level = LOG_DEFAULT;
  cfg->write_trace = false;
}

/**
 * @brief Validates configuration and logs any errors.
 *
 * @param[in] cfg global configuration
 * @param[in] world_size number of MPI processes in world
 * @return int status (0: ok, 1: error)
 */
int validate_config(global_config_t *cfg, int world_size) {
  /* Infected */
  if (cfg->inf_individuals > cfg->num_individuals) {
    log_error("Infected individuals exceed population size (%lu > %lu)",
              cfg->inf_individuals, cfg->num_individuals);
    return 1;
  }
  /* World dimensions */
  if (cfg->world_w < cfg->country_w || cfg->world_l < cfg->country_l) {
    log_error("Country dimensions cannot exceed world dimensions");
    return 1;
  }
  if (cfg->world_w % cfg->country_w != 0) {
    log_error("Country width must divide world width");
    return 1;
  }
  if (cfg->world_l % cfg->country_l != 0) {
    log_error("Country length must divide world length");
    return 1;
  }
  int num_countries =
      (cfg->world_w / cfg->country_w) * (cfg->world_l / cfg->country_l);
  if (world_size != num_countries) {
    log_error(
        "Number of processes does not match number of countries. "
        "Expected %d, got %d",
        num_countries, world_size);
    return 1;
  }
  /* Velocity */
  if (cfg->velocity <= 0.) {
    log_error("Velocity must be non-negative");
    return 1;
  }
  /* Spreading distance */
  if (cfg->spreading_distance <= 0.) {
    log_error("Spreading distance must be non-negative");
    return 1;
  }
  /* Simulation step */
  if (cfg->t_step > DAY) {
    log_error("Simulation step cannot be longer than one day");
    return 1;
  }
  /* Relation between movement and time step */
  if (cfg->t_step * cfg->velocity > MIN(cfg->country_w, cfg->country_l)) {
    log_error(
        "The movement at each step is larger than a country: t * v "
        "= %f > %lu",
        cfg->t_step * cfg->velocity, MIN(cfg->country_w, cfg->country_l));
    return 1;
  }

  /* If we got here the configuration is valid */
  return 0;
}

/**
 * @brief Logs the configuration with level INFO.
 *
 * @param[in] cfg configuration
 */
void log_config(global_config_t *cfg) {
  log_info(
      "\n--------------------\nGlobal configuration\n--------------------\n "
      "num_individuals %lu\n "
      "inf_individuals %lu\n world_w %lu\n world_l %lu\n country_w %lu\n "
      "country_l %lu\n velocity %f\n spreading_distance %f\n t_infection "
      "%lu\n t_recovery %lu\n t_immunity %lu\n t_step %lu\n t_target "
      "%lu\n rand_seed %u\n log_level %s\n write_trace "
      "%d\n--------------------\n",
      cfg->num_individuals, cfg->inf_individuals, cfg->world_w, cfg->world_l,
      cfg->country_w, cfg->country_l, cfg->velocity, cfg->spreading_distance,
      cfg->t_infection, cfg->t_recovery, cfg->t_immunity, cfg->t_step,
      cfg->t_target, cfg->rand_seed, log_level_string(cfg->log_level),
      cfg->write_trace);
}