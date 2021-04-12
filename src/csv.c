#include "csv.h"

/**
 * @brief Create a csv file for individual's details and write header
 *
 * @param[in] directory path of the directory where to store the file, not NULL
 * @param[in] country country of the calling process
 * @return FILE* file pointer with write access, NULL if error
 */
FILE *create_detail_csv(char *directory, int country) {
  char *path = malloc(PATH_MAX * sizeof(char));
  /* Determine the filename and open the file */
  sprintf(path, "%s/detail_%d.csv", directory, country);
  FILE *csv = fopen(path, "w");
  if (csv) {
    /* Write the header */
    fprintf(csv, "country,t,id,pos_x,pos_y,displ_x,displ_y,status,t_status\n");
  } else {
    log_error("Cannot open file \"%s\" for writing", path);
  }
  free(path);
  return csv;
}

/**
 * @brief Write an individual's detail in the given csv file
 *
 * @param[in] csv csv file pointer, not NULL
 * @param[in] ind individual to be printed
 * @param[in] country country of the calling process
 * @param[in] t current time
 */
void detail_csv_write_line(FILE *csv, individual_t *ind, int country,
                           unsigned long t) {
  fprintf(csv, "%d,%lu,%lu,%.3f,%.3f,%.3f,%.3f,%s,%lu\n", country, t, ind->id,
          ind->pos[0], ind->pos[1], ind->displ[0], ind->displ[1],
          individual_status_string(ind->status), ind->t_status);
}