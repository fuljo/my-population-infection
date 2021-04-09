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