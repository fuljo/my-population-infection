#include "csv.h"

/**
 * @brief Create a csv file for individual's details and write header
 *
 * @param[in] directory path of the directory where to store the file, not NULL
 * @param[in] country country of the calling process
 * @return FILE* file pointer with write access, NULL if error
 */
FILE *create_trace_csv(char *directory, int country) {
  char *path = malloc(PATH_MAX * sizeof(char));
  /* Determine the filename and open the file */
  sprintf(path, "%s/trace_%d.csv", directory, country);
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
 * @brief Write in the given csv file the details of a list of individuals
 *
 * @param[in] csv csv file pointer, not NULL
 * @param[in] individuals list of individuals to be printed
 * @param[in] country country of the calling process
 * @param[in] t current time
 */
void trace_csv_write_step(FILE *csv, individual_list_t *individuals,
                          int country, unsigned long t) {
  individual_t *ind;
  INDIVIDUAL_FOREACH(ind, individuals) {
    fprintf(csv, "%d,%lu,%lu,%.3f,%.3f,%.3f,%.3f,%s,%lu\n", country, t, ind->id,
            ind->pos[0], ind->pos[1], ind->displ[0], ind->displ[1],
            individual_status_string(ind->status), ind->t_status);
  }
}

/**
 * @brief Create a csv file for summary and write header
 *
 * @param[in] directory path of the directory where to store the file, not NULL
 * @return FILE* file pointer with write access, NULL if error
 */
FILE *create_summary_csv(char *directory) {
  char *path = malloc(PATH_MAX * sizeof(char));
  /* Determine the filename and open the file */
  sprintf(path, "%s/summary.csv", directory);
  FILE *csv = fopen(path, "w");
  if (csv) {
    /* Write the header */
    fprintf(csv, "day,country,subsceptible,infected,immune\n");
  } else {
    log_error("Cannot open file \"%s\" for writing", path);
  }
  free(path);
  return csv;
}

/**
 * @brief Write the given collection of summaries in the summary file
 *
 * @param csv csv file pointer, not null
 * @param summaries array of summaries to write
 * @param len length of the array
 * @param day number of the day that has just ended (zero-based by convention)
 */
void summary_csv_write_day(FILE *csv, summary_t summaries[], size_t len,
                           unsigned long day) {
  summary_t *s;
  for (size_t i = 0; i < len; i++) {
    s = &summaries[i];
    fprintf(csv, "%lu,%lu,%lu,%lu,%lu\n", day, i, s->subsceptible, s->infected,
            s->immune);
  }
}