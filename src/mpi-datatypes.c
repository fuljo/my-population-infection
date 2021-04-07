#include "mpi-datatypes.h"

/**
 * @brief Create the MPI version of the struct global config datatype.
 *
 * The type is both created and committed, but needs to be freed after use.
 *
 * @return MPI_Datatype
 */
MPI_Datatype create_type_mpi_global_config() {
  MPI_Datatype mpi_global_config;
  struct global_config cfg;
  /**
   * We use three blocks:
   * - MPI_UNSIGNED_LONG (6 elements)
   * - MPI_DOUBLE (2 elements)
   * - MPI_UNSIGNED_LONG (9 elements)
   */
  int num_blocks = 3;
  const int block_lengths[] = {6, 2, 9};
  const MPI_Aint displacements[] = {
      0,
      (size_t) & (cfg.num_individuals) - (size_t) & (cfg),
      (size_t) & (cfg.velocity) - (size_t) & (cfg),
      (size_t) & (cfg.t_infection) - (size_t) & (cfg),
  };
  MPI_Datatype block_types[] = {
      MPI_UNSIGNED_LONG,
      MPI_DOUBLE,
      MPI_UNSIGNED_LONG,
  };
  MPI_Type_create_struct(num_blocks, block_lengths, displacements, block_types,
                         &mpi_global_config);
  MPI_Type_commit(&mpi_global_config);
  return mpi_global_config;
}