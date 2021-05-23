#pragma once

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "individual.h"
#include "utils.h"

FILE *create_trace_csv(const char *directory, int country);

void trace_csv_write_step(FILE *csv, individual_list_t *individuals,
                          int country, unsigned long t);

FILE *create_summary_csv(const char *directory);

void summary_csv_write_day(FILE *csv, summary_t summaries[], size_t len,
                           unsigned long day);