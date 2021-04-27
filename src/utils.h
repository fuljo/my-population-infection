#pragma once

#include <math.h>
#include <stdlib.h>
#include <sys/queue.h>

#include "log.h"

#define SECOND 1
#define MINUTE 60
#define DAY 60 * 60 * 24

/* How many elements will be allocated when a dynamic array reaches its capacity
 */
#define DYN_ARRAY_CHUNK 64

#define ROOT_RANK 0

#define RAND_DOUBLE(offset, range) \
  (offset) + ((double)rand() * (range) / RAND_MAX)

#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)

#define SLIST_REMOVE_AFTER(elm, field)                             \
  do {                                                             \
    (elm)->field.sle_next = (elm)->field.sle_next->field.sle_next; \
  } while (0)

#define DYN_ARRAY_APPEND(elm, arr, len, capacity, type) \
  do {                                                  \
    if (len >= capacity) {                              \
      capacity += DYN_ARRAY_CHUNK;                      \
      arr = realloc(arr, sizeof(type) * capacity);    \
    }                                                   \
    arr[len] = elm;                                     \
    len++;                                              \
  } while (0)