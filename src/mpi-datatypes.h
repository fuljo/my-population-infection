#include <mpi.h>

#include "config.h"

/**
 * @brief Create the MPI version of the struct global config datatype.
 *
 * The type is both created and committed, but needs to be freed after use.
 *
 * @return MPI_Datatype
 */
MPI_Datatype create_type_mpi_global_config();