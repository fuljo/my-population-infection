#include <argp.h>
#include <argz.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "utils.h"

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
  printf("[INFO] Hello from node %d out of %d\n", rank, world_size);

  /* Read and parse command-line configuration */
  if (rank == HEAD_RANK) {
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
        {0},
    };
    /* Define program description */
    struct argp argp = {options, parse_opt, NULL,
                        "A simple model for virus spreading."};

    /* Read command-line options and arguments */
    struct arguments arguments;
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
      fprintf(stderr, "[ERR ] Error while parsing CLI arguments\n");
      return 1;
    }

    print_config(&arguments.config);

    /* Validate configuration */
    if (validate_config(&arguments.config, world_size) != 0) {
      MPI_Finalize();
      return 1; /* Errors are printed inside the function */
    }
  }

  MPI_Finalize();
  return 0;
}
