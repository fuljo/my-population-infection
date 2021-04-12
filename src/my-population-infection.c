#include <argp.h>
#include <argz.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "config.h"
#include "csv.h"
#include "mpi-datatypes.h"
#include "utils.h"
#include "world.h"

/* Function prototypes */
void initialize_individuals(global_config_t *cfg, int num_countries,
                            individual_list_t *subsceptible_individuals,
                            individual_list_t *infected_individuals);

int main(int argc, char **argv) {
  /* -------------------------------------------------------------------------*/
  /* Initialization                                                           */
  /* -------------------------------------------------------------------------*/

  /* Initialize MPI */
  MPI_Init(&argc, &argv);

  /* Get information about MPI environment */
  int world_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  log_info("Hello from node %d out of %d", rank, world_size);

  /* Create custom MPI datatypes */
  MPI_Datatype mpi_global_config = create_type_mpi_global_config();
  MPI_Datatype mpi_individual = create_type_mpi_individual();

  /* Read and parse command-line configuration */
  global_config_t cfg;
  if (rank == ROOT_RANK) {
    /* Define command-line options */
    struct argp_option options[] = {
        {0, 0, 0, 0, "Population options", 1},
        {"num-individuals", 'N', "INT", 0, "Number of individuals"},
        {"inf-individuals", 'I', "INT", 0,
         "Number of initially infected individuals"},
        {0, 0, 0, 0, "World options (lengths in meters)", 2},
        {"world-width", 'W', "INT", 0, "Width of the world rectangle"},
        {"world-length", 'L', "INT", 0, "Length of the world rectangle"},
        {"country-width", 'w', "INT", 0,
         "Width of a single country, must divide W"},
        {"country-length", 'l', "INT", 0,
         "Length of a single country, must divide L"},
        {0, 0, 0, 0, "Individual options (times in seconds)", 3},
        {"velocity", 'v', "FLOAT", 0, "Moving speed for and individual in m/s"},
        {"spreading-distance", 'd', "FLOAT", 0,
         "Maximum spreading distance in meters"},
        {"t-infection", 333, "INT", OPTION_ARG_OPTIONAL,
         "Minimum continuous exposure time for getting infected"},
        {"t-recovery", 444, "INT", OPTION_ARG_OPTIONAL,
         "Time needed for recovering"},
        {"t-immunity", 555, "INT", OPTION_ARG_OPTIONAL,
         "Duration of immunity period after recovering"},
        {0, 0, 0, 0, "Simulation options", 4},
        {"sim-step", 666, "INT", 0, "Simulation step in seconds"},
        {"sim-length", 777, "INT", 0, "Length of the simulation in days"},
        {"rand-seed", 888, "INT", OPTION_ARG_OPTIONAL,
         "Seed for PRNG. By default initialized to time(0)"},
        {0, 0, 0, 0, "Logging options", 5},
        {"log-level", 999, "[TRACE|DEBUG|INFO|WARN|ERROR|FATAL]", 0,
         "Logging level"},
        {0},
    };
    /* Define program description */
    struct argp argp = {options, parse_opt, NULL,
                        "A simple model for virus spreading."};

    /* Read command-line options and arguments */
    struct arguments arguments;
    init_config_default(&arguments.config);
    if (argp_parse(&argp, argc, argv, 0, 0, &arguments) == 0) {
      const char *prev = NULL;
      char *word;
      while ((word = argz_next(arguments.argz, arguments.argz_len, prev))) {
        printf("%s", word);
        prev = word;
      }
      printf("\n");
      free(arguments.argz);
    } else {
      MPI_Finalize();
      log_error("Error while parsing CLI arguments\n");
      return 1;
    }

    cfg = arguments.config;
    print_config(&cfg);

    /* Validate configuration */
    if (validate_config(&cfg, world_size) != 0) {
      MPI_Finalize();
      return 1; /* Errors are printed inside the function */
    }
  }

  /* Broadcast the configuration to all processes */
  MPI_Bcast(&cfg, 1, mpi_global_config, 0, MPI_COMM_WORLD);
  /* Initialize random number generator */
  /* NOTE: It is important to give variability between countries */
  srand(cfg.rand_seed + rank);

  /* Calculate the total number of countries */
  unsigned int num_countries =
      (cfg.world_w / cfg.country_w) * (cfg.world_l / cfg.country_l);

  /* Create empty lists of individuals */
  individual_list_t subsceptible_individuals = create_individual_list();
  individual_list_t infected_individuals = create_individual_list();
  individual_list_t immune_individuals = create_individual_list();
  individual_list_t gc_individuals = create_individual_list(); /* recycle bin */

  /* Distribute individuals between countries and initialize them */
  initialize_individuals(&cfg, num_countries, &subsceptible_individuals,
                         &infected_individuals);

  /* -------------------------------------------------------------------------*/
  /* Cleanup                                                                  */
  /* -------------------------------------------------------------------------*/
  MPI_Type_free(&mpi_global_config);
  MPI_Finalize();
  return 0;
}

/**
 * @brief Distributes individuals between countries and initializes them
 *
 * The population (individuals and infected) defined in the configuration is
 * uniformly distributed among countries. Then each country generates its
 * individuals and assign to each of them:
 *  - a random position
 *  - a displacement vector with random direction
 *  - status \c NOT_EXPOSED or \c INFECTED according to the distribution
 *  - <tt>t_status = 0</tt>
 *
 * They are inserted in the correct list according to their status.
 *
 * @param cfg[in] global configuration
 * @param num_countries[in] number of countries
 * @param subsceptible_individuals[out] head of the list where \c NOT_EXPOSED
 * individuals will be inserted
 * @param infected_individuals[out] head of the list where \c INFECTED
 * individuals will be inserted
 */
void initialize_individuals(global_config_t *cfg, int num_countries,
                            individual_list_t *subsceptible_individuals,
                            individual_list_t *infected_individuals) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int cols = (cfg->world_w / cfg->country_w);
  int col = rank % cols;
  int row = rank / cols;

  /* Distribute individuals and infected between countries */
  unsigned long num_individuals, num_infected;
  if (rank == ROOT_RANK) {
    unsigned long num_individuals_by_country[num_countries];
    unsigned long num_infected_by_country[num_countries];
    distribute_population_uniform(cfg->num_individuals, num_countries,
                                  num_individuals_by_country);
    distribute_population_uniform(cfg->inf_individuals, num_countries,
                                  num_infected_by_country);
    /* Scatter them among all countries */
    MPI_Scatter(num_individuals_by_country, 1, MPI_UNSIGNED_LONG,
                &num_individuals, 1, MPI_UNSIGNED_LONG, ROOT_RANK,
                MPI_COMM_WORLD);
    MPI_Scatter(num_infected_by_country, 1, MPI_UNSIGNED_LONG, &num_infected, 1,
                MPI_UNSIGNED_LONG, ROOT_RANK, MPI_COMM_WORLD);
  } else {
    /* Receive the scattered values */
    MPI_Scatter(NULL, 0, NULL, &num_individuals, 1, MPI_UNSIGNED_LONG,
                ROOT_RANK, MPI_COMM_WORLD);
    MPI_Scatter(NULL, 0, NULL, &num_infected, 1, MPI_UNSIGNED_LONG, ROOT_RANK,
                MPI_COMM_WORLD);
  }

  /* Determine the id of the first individual for each country */
  /* The ids are assigned incrementally, so we can perform an exclusive scan */
  unsigned long initial_id = 0; /* otherwise rank 0 receives undefined */
  MPI_Exscan(&num_individuals, &initial_id, 1, MPI_UNSIGNED_LONG, MPI_SUM,
             MPI_COMM_WORLD);

  log_debug("Rank %d -- individuals=%lu, infected=%lu, initial id=%lu", rank,
            num_individuals, num_infected, initial_id);

  if (num_individuals == 0) {
    return;
  }

  /* Initialize each individual */
  unsigned long x_min = col * cfg->country_w;
  unsigned long y_min = row * cfg->country_l;
  individual_t *ind;
  unsigned long id = initial_id;
  double theta;
  /* The first individual must not be made infected if (N - 1) % I = 0 to have
   * the right number of infected people */
  if (num_infected > 0 && num_individuals % num_infected == 1) {
    ind = malloc(sizeof(individual_t));
    ind->id = id;
    ind->t_status = 0;

    /* Random position inside the country */
    ind->pos[0] = RAND_DOUBLE(x_min, cfg->country_w);
    ind->pos[1] = RAND_DOUBLE(y_min, cfg->country_l);

    /* Random direction to compute displacement */
    theta = RAND_DOUBLE(0, 2 * M_PI);
    ind->displ[0] = cfg->t_step * cfg->velocity * cos(theta);
    ind->displ[1] = cfg->t_step * cfg->velocity * sin(theta);

    /* Not infected */
    ind->status = NOT_EXPOSED;
    insert_individual(subsceptible_individuals, ind);

    id++;
  }
  while (id < initial_id + num_individuals) {
    individual_t *ind = malloc(sizeof(individual_t));
    ind->id = id;
    ind->t_status = 0;
    ind->pos[0] = RAND_DOUBLE(x_min, cfg->country_w);
    ind->pos[1] = RAND_DOUBLE(y_min, cfg->country_l);
    theta = RAND_DOUBLE(0, 2 * M_PI);
    ind->displ[0] = cfg->t_step * cfg->velocity * cos(theta);
    ind->displ[1] = cfg->t_step * cfg->velocity * sin(theta);

    /* Set status and insert into correct list */
    if (num_infected > 0 && id % num_infected == 0) {
      ind->status = INFECTED;
      insert_individual(infected_individuals, ind);
    } else {
      ind->status = NOT_EXPOSED;
      insert_individual(subsceptible_individuals, ind);
    }

    id++;
  }
}
