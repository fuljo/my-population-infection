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

void update_exposure(double spreading_distance,
                     individual_list_t *subsceptible_individuals,
                     individual_list_t *infected_individuals);

void update_status(global_config_t *cfg,
                   individual_list_t *subsceptible_individuals,
                   individual_list_t *infected_individuals,
                   individual_list_t *immune_individuals);

void update_position(global_config_t *cfg, individual_list_t *individuals,
                     individual_list_t *gc_individuals,
                     individual_t *migrated_out[], size_t migrated_out_len[],
                     size_t migrated_out_capacity[], limits_t *limits,
                     int neighbors[]);

limits_t calculate_country_limits(global_config_t *cfg, int num_countries,
                                  int rank);

void calculate_neighbors(int neighbors[], global_config_t *cfg,
                         int num_countries, int rank);

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
  /* Set log level */
  log_set_level(cfg.log_level);
  /* Initialize random number generator */
  /* NOTE: It is important to give variability between countries */
  srand(cfg.rand_seed + rank);

  /* Calculate the total number of countries */
  const unsigned int num_countries =
      (cfg.world_w / cfg.country_w) * (cfg.world_l / cfg.country_l);

  /* Calculate country limits */
  limits_t limits = calculate_country_limits(&cfg, num_countries, rank);
  /* Calculate indices of neighbors (-1 if none) */
  int neighbors[NEIGHBOR_COUNT];
  calculate_neighbors(neighbors, &cfg, num_countries, rank);

  /* Create empty lists of individuals */
  individual_list_t subsceptible_individuals = create_individual_list();
  individual_list_t infected_individuals = create_individual_list();
  individual_list_t immune_individuals = create_individual_list();
  individual_list_t gc_individuals = create_individual_list(); /* recycle bin */

  /* Create buffers to move individuals from/to neighbor countries */
  individual_t *migrated_out[NEIGHBOR_COUNT] = {NULL};
  individual_t *migrated_in[NEIGHBOR_COUNT] = {NULL};
  size_t migrated_out_len[NEIGHBOR_COUNT] = {0};
  size_t migrated_out_capacity[NEIGHBOR_COUNT] = {0};
  size_t migrated_in_len[NEIGHBOR_COUNT] = {0};
  size_t migrated_in_capacity[NEIGHBOR_COUNT] = {0};

  /* Distribute individuals between countries and initialize them */
  initialize_individuals(&cfg, num_countries, &subsceptible_individuals,
                         &infected_individuals);

  /* -------------------------------------------------------------------------*/
  /* Main loop                                                                */
  /* -------------------------------------------------------------------------*/
  FILE *detail_csv = create_detail_csv("./results", rank);
  for (unsigned long t = 0; t < cfg.t_target; t += cfg.t_step) {
    log_debug("Rank %d -- t = %lu", rank, t);
    /* Update exposure of subsceptible individuals */
    update_exposure(cfg.spreading_distance, &subsceptible_individuals,
                    &infected_individuals);

    /* Update the status of all individuals based on t_status and move them
       into the correct list */
    update_status(&cfg, &subsceptible_individuals, &infected_individuals,
                  &immune_individuals);

    /* Move the individuals according to the displacement, perform bouncing and
     * populate the migrated_out buffers */
    update_position(&cfg, &subsceptible_individuals, &gc_individuals,
                    migrated_out, migrated_out_len, migrated_out_capacity,
                    &limits, neighbors);
    update_position(&cfg, &infected_individuals, &gc_individuals, migrated_out,
                    migrated_out_len, migrated_out_capacity, &limits,
                    neighbors);
    update_position(&cfg, &immune_individuals, &gc_individuals, migrated_out,
                    migrated_out_len, migrated_out_capacity, &limits,
                    neighbors);
  }
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
 * @param[in] cfg global configuration
 * @param[in] num_countries number of countries
 * @param[out] subsceptible_individuals  head of the list where \c NOT_EXPOSED
 * individuals will be inserted
 * @param[out] infected_individuals  head of the list where \c INFECTED
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
  const unsigned long x_min = col * cfg->country_w;
  const unsigned long y_min = row * cfg->country_l;
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
    INDIVIDUAL_INSERT(subsceptible_individuals, ind);

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
      INDIVIDUAL_INSERT(infected_individuals, ind);
    } else {
      ind->status = NOT_EXPOSED;
      INDIVIDUAL_INSERT(subsceptible_individuals, ind);
    }

    id++;
  }
}

/**
 * @brief Compute the exposure status of subsceptible individuals
 *
 * @pre All subsceptible individuals have <tt>status = NOT_EXPOSED<\tt>
 * @post Each subsceptible individual is flagged as \c EXPOSED if there is at
 * least one \c INFECTED individual in a \c spreading_distance radius from him;
 * otherwise it remains \c NOT_EXPOSED . No items are inserted or removed from
 * the lists.
 *
 * @param[in] spreading_distance inclusive distance to be considered exposed
 * @param[in,out] subsceptible_individuals list of all subsceptible individuals
 * @param[in] infected_individuals list of all infected individuals
 */
void update_exposure(double spreading_distance,
                     individual_list_t *subsceptible_individuals,
                     individual_list_t *infected_individuals) {
  individual_t *i, *j;
  /* We check each subsceptible individual against infected individual */
  INDIVIDUAL_FOREACH(i, subsceptible_individuals) {
    INDIVIDUAL_FOREACH(j, infected_individuals) {
      if (INDIVIDUAL_DISTANCE(i, j) <= spreading_distance) {
        /* As soon as one match is found, we can go on to the next i */
        i->status = EXPOSED;
        break;
      }
    }
  }
}

/**
 * @brief Updates the status of each individual in the given lists and moves it
 * to the correct list
 *
 * @pre All indivuduals are in the correct list according to their status.
 * All \c subsceptible_individuals that are actually exposed have
 * <tt>status == EXPOSED</tt>
 *
 * @post All individuals are in the correct list according to their status.
 * All \c subsceptible_individuals have \c status reset to \c EXPOSED.
 * If the individual was \c NOT_EXPOSED , \c t_status is reset to zero.
 * In all other cases \c t_status is incremented by \c t_step , except if there
 * is a status change, where it is reset to zero.
 *
 * @param[in] cfg global configuration
 * @param[in,out] subsceptible_individuals list of all subsceptible individuals
 * @param[in,out] infected_individuals list of all infected individuals
 * @param[in,out] immune_individuals list of all immune individuals
 */
void update_status(global_config_t *cfg,
                   individual_list_t *subsceptible_individuals,
                   individual_list_t *infected_individuals,
                   individual_list_t *immune_individuals) {
  /* Save the first element of each list, so we don't re-process elements that
   * are inserted in the meantime */
  individual_t *sub = INDIVIDUAL_FIRST(subsceptible_individuals);
  individual_t *inf = INDIVIDUAL_FIRST(infected_individuals);
  individual_t *imm = INDIVIDUAL_FIRST(immune_individuals);
  individual_t *prev, *next;

  /* Subsceptible -> Infected */
  prev = NULL;
  while (sub) {
    next = INDIVIDUAL_NEXT(sub);
    if (sub->status == EXPOSED) { /* EXPOSED */
      sub->t_status += cfg->t_step;
      if (sub->t_status >= cfg->t_infection) {
        /* The individual becomes infected */
        sub->status = INFECTED;
        sub->t_status = 0;
        /* Remove it from the current list */
        if (prev) {
          INDIVIDUAL_REMOVE_AFTER(prev);
        } else {
          INDIVIDUAL_REMOVE_HEAD(subsceptible_individuals);
        }
        /* Put it in the other list */
        INDIVIDUAL_INSERT(infected_individuals, sub);
      } else {
        /* If we didn't move the current element, we can advance prev */
        prev = sub;
      }
    } else { /* NOT_EXPOSED */
      sub->t_status = 0;
      prev = sub;
    }
    sub = next;
  }

  /* Infected -> Immune */
  prev = NULL;
  while (inf) {
    next = INDIVIDUAL_NEXT(inf);
    inf->t_status += cfg->t_step;
    if (inf->t_status >= cfg->t_infection) {
      /* The individual becomes immune */
      inf->status = IMMUNE;
      inf->t_status = 0;
      /* Remove it from the current list */
      if (prev) {
        INDIVIDUAL_REMOVE_AFTER(prev);
      } else {
        INDIVIDUAL_REMOVE(infected_individuals, inf);
      }
      /* Put it in the other list */
      INDIVIDUAL_INSERT(immune_individuals, inf);
    } else {
      /* If we didn't move the current element, we can advance prev */
      prev = inf;
    }
    inf = next;
  }

  /* Immune -> Subsceptible */
  prev = NULL;
  while (imm) {
    next = INDIVIDUAL_NEXT(imm);
    imm->t_status += cfg->t_step;
    if (imm->t_status >= cfg->t_immunity) {
      /* The individual becomes immune */
      imm->status = NOT_EXPOSED;
      imm->t_status = 0;
      /* Remove it from the current list */
      if (prev) {
        INDIVIDUAL_REMOVE_AFTER(prev);
      } else {
        INDIVIDUAL_REMOVE(immune_individuals, imm);
      }
      /* Put it in the other list */
      INDIVIDUAL_INSERT(subsceptible_individuals, imm);
    } else {
      /* If we didn't move the current element, we can advance prev */
      prev = imm;
    }
    imm = next;
  }
}

/**
 * @brief Updates the position of each individual in a list and moves
 * individuals that exited the country to an outbound buffer
 *
 * @param[in] cfg global configuration
 * @param[in,out] individuals list of individuals to be processed
 * @param[in,out] gc_individuals list of garbage-collected individuals
 * @param[in,out] migrated_out buffer indexed by cardinal point where to put
 * outbound individuals, each position must be dynamically allocated
 * @param[in,out] migrated_out_len currently used size of the buffer (in number
 * of individuals) for each position
 * @param[in,out] migrated_out_capacity current capacity of each position of the
 * buffer
 * @param[in] limits limits of the country
 * @param[in] neighbors array of indinces of neighbor countries, one for each of
 * the eight cardinal directions
 */
void update_position(global_config_t *cfg, individual_list_t *individuals,
                     individual_list_t *gc_individuals,
                     individual_t *migrated_out[], size_t migrated_out_len[],
                     size_t migrated_out_capacity[], limits_t *limits,
                     int neighbors[]) {
  individual_t *ind, *prev, *next;

  /* Out of bound flag: the bits are indexed according to cardinal_point_t */
  unsigned char out_flag;

  /* Cardinal point of the destination country */
  int dest;

  /* Iterate over the list */
  ind = INDIVIDUAL_FIRST(individuals);
  prev = NULL;
  while (ind) {
    out_flag = 0;
    next = INDIVIDUAL_NEXT(ind);

    /* Move of the given displacement */
    ind->pos[0] += ind->displ[0];
    ind->pos[1] += ind->displ[1];

    /* Calculate residuals w.r.t the boundaries */
    double res_xmin = ind->pos[0] - limits->xmin;
    double res_xmax = ind->pos[0] - limits->xmax;
    double res_ymin = ind->pos[1] - limits->ymin;
    double res_ymax = ind->pos[1] - limits->ymax;

    /* Check if out-of-bound horizontally */
    if (res_xmin < 0) { /* West */
      if (neighbors[WEST] < 0) {
        /* Out of world => bounce */
        ind->pos[0] += -2 * res_xmin;
        ind->displ[0] = -ind->displ[0];
      } else {
        out_flag += 1 << WEST;
      }
    } else if (res_xmax >= 0) { /* East */
      if (neighbors[EAST] < 0) {
        /* Out of world => bounce */
        ind->pos[0] += -2 * res_xmax;
        ind->displ[0] = -ind->displ[0];
      } else {
        out_flag += 1 << EAST;
      }
    }

    /* Check if out-of-bound vertically */
    if (res_ymin < 0) { /* South */
      if (neighbors[SOUTH] < 0) {
        /* Out of world => bounce */
        ind->pos[1] += -2 * res_ymin;
        ind->displ[1] = -ind->displ[1];
      } else {
        out_flag += 1 << SOUTH;
      }
    } else if (res_ymax >= 0) { /* North */
      if (neighbors[NORTH] < 0) {
        /* Out of world => bounce */
        ind->pos[1] += -2 * res_ymax;
        ind->displ[1] = -ind->displ[1];
      } else {
        out_flag += 1 << NORTH;
      }
    }

    /* Send out-of-bound individuals to another country */
    if (out_flag) {
      /* Determine destionation */
      dest = decode_cardinal_point_flag(out_flag);
      /* Remove from local list */
      if (prev) {
        INDIVIDUAL_REMOVE_AFTER(prev);
      } else {
        INDIVIDUAL_REMOVE_HEAD(individuals);
      }
      /* Copy to migration buffer */
      DYN_ARRAY_APPEND(*ind, migrated_out[dest], migrated_out_len[dest],
                       migrated_out_capacity[dest], individual_t);
      /* Insert old list element into garbage collection list to be reused */
      INDIVIDUAL_INSERT(gc_individuals, ind);
    } else {
      /* If we didn't remove the current element, we can advance prev */
      prev = ind;
    }

    /* Prepare for next iteration */
    ind = next;
  }

/**
 * @brief Calculate the min and max x and y values for a country.
 *
 * @param[in] cfg global configuration
 * @param[in] num_countries total number of countries
 * @param[in] rank rank of this country
 * @return limits_t a struct with the calculated limits
 */
limits_t calculate_country_limits(global_config_t *cfg, int num_countries,
                                  int rank) {
  /* Determine row and column of this country */
  int cols = (cfg->world_w / cfg->country_w);
  int col = rank % cols;
  int row = rank / cols;
  /* Build the struct with the values */
  limits_t limits = {
      col * cfg->country_w,       /* xmin */
      (col + 1) * cfg->country_w, /* xmax */
      row * cfg->country_l,       /* ymin */
      (row + 1) * cfg->country_l, /* ymax */
  };
  return limits;
}

/**
 * @brief Calculates the indices of the neighbors of this country
 *
 * @param[out] neighbors array where the results will be stored
 * @param[in] cfg global configuration
 * @param[in] num_countries total number of countries
 * @param[in] rank rank of this country
 */
void calculate_neighbors(int neighbors[], global_config_t *cfg,
                         int num_countries, int rank) {
  /* Determine row and column of this country */
  int cols = (cfg->world_w / cfg->country_w);
  int rows = (cfg->world_l / cfg->country_l);
  int col = rank % cols;
  int row = rank / cols;

  /* Set indices as if this were an internal node */
  neighbors[SOUTH_EAST] = rank - cols - 1;
  neighbors[SOUTH] = rank - cols;
  neighbors[SOUTH_WEST] = rank - cols + 1;
  neighbors[NORTH_EAST] = rank - cols - 1;
  neighbors[NORTH] = rank - cols;
  neighbors[NORTH_WEST] = rank - cols + 1;
  neighbors[WEST] = rank - 1;
  neighbors[EAST] = rank + 1;

  /* Correct the indices by considering the world boundaries */
  if (row == 0) { /* bottom row */
    neighbors[SOUTH_EAST] = -1;
    neighbors[SOUTH] = -1;
    neighbors[SOUTH_WEST] = -1;
  }
  if (row == rows - 1) { /* top row */
    neighbors[NORTH_EAST] = -1;
    neighbors[NORTH] = -1;
    neighbors[NORTH_WEST] = -1;
  }

  if (col == 0) { /* first column */
    neighbors[SOUTH_WEST] = -1;
    neighbors[WEST] = -1;
    neighbors[NORTH_WEST] = -1;
  }
  if (col == cols - 1) { /* last column */
    neighbors[NORTH_EAST] = -1;
    neighbors[EAST] = -1;
    neighbors[NORTH_EAST] = -1;
  }
}