#pragma once

#include <math.h>
#include <stdlib.h>
#include <sys/queue.h>

#include "log.h"

#define SECOND 1
#define MINUTE 60
#define DAY 60 * 60 * 24

#define ROOT_RANK 0

#define RAND_DOUBLE(offset, range) \
  (offset) + ((double)rand() * (range) / RAND_MAX)

#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)
