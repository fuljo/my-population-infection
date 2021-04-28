#pragma once

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "individual.h"
#include "utils.h"

FILE *create_detail_csv(char *directory, int country);

void detail_csv_write_line(FILE *csv, individual_t *ind, int country,
                           unsigned long t);

FILE *create_summary_csv(char *directory);

void summary_csv_write_day(FILE *csv, summary_t summaries[], size_t len,
                           unsigned long day);