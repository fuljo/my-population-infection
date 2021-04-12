#include "mpi-datatypes.h"

/**
 * @brief Create the MPI version of the global_config_t datatype.
 *
 * The type is both created and committed, but needs to be freed after use.
 *
 * @return MPI_Datatype
 */
MPI_Datatype create_type_mpi_global_config() {
  MPI_Datatype mpi_global_config;
  global_config_t cfg;
  /**
   * We use four blocks:
   * - MPI_UNSIGNED_LONG (6 elements)
   * - MPI_DOUBLE (2 elements)
   * - MPI_UNSIGNED_LONG (9 elements)
   * - MPI_INT (1 element)
   */
  int num_blocks = 4;
  const int block_lengths[] = {6, 2, 9, 1};
  const MPI_Aint displacements[] = {
      0,
      (size_t) & (cfg.num_individuals) - (size_t) & (cfg),
      (size_t) & (cfg.velocity) - (size_t) & (cfg),
      (size_t) & (cfg.t_infection) - (size_t) & (cfg),
      (size_t) & (cfg.log_level) - (size_t) & (cfg),
  };
  MPI_Datatype block_types[] = {
      MPI_UNSIGNED_LONG,
      MPI_DOUBLE,
      MPI_UNSIGNED_LONG,
      MPI_INT,
  };
  MPI_Type_create_struct(num_blocks, block_lengths, displacements, block_types,
                         &mpi_global_config);
  MPI_Type_commit(&mpi_global_config);
  return mpi_global_config;
}

/**
 * @brief Create the MPI version of the individual_t datatype.
 *
 * The type is both created and committed, but needs to be freed after use.
 *
 * @return MPI_Datatype
 */
MPI_Datatype create_type_mpi_individual() {
  MPI_Datatype mpi_individual;
  individual_t ind;
  /**
   * We use five blocks:
   * - MPI_UNSIGNED_LONG (1 element)
   * - MPI_DOUBLE (4 elements)
   * - MPI_INT (1 element)
   * - MPI_UNSIGNED_LONG (1 element)
   * - MPI_AINT (1 element)
   */
  int num_blocks = 5;
  const int block_lengths[] = {1, 4, 1, 1, 1};
  const MPI_Aint displacements[] = {
      0,
      (size_t) & (ind.pos) - (size_t) & (ind),
      (size_t) & (ind.status) - (size_t) & (ind),
      (size_t) & (ind.t_status) - (size_t) & (ind),
      (size_t) & (ind.individuals) - (size_t) & (ind),
  };
  MPI_Datatype block_types[] = {
      MPI_UNSIGNED_LONG,
      MPI_DOUBLE,
      MPI_INT,
      MPI_UNSIGNED_LONG,
      MPI_AINT,
  };
  MPI_Type_create_struct(num_blocks, block_lengths, displacements, block_types,
                         &mpi_individual);
  MPI_Type_commit(&mpi_individual);

  return mpi_individual;
}