#include "individual.h"

/**
 * @brief Returns a string representation of the given status code
 *
 * @param status
 * @return char*
 */
char *individual_status_string(int status) {
  switch (status) {
    case NOT_EXPOSED:
      return "NOT_EXPOSED";
    case EXPOSED:
      return "EXPOSED";
    case INFECTED:
      return "INFECTED";
    case IMMUNE:
      return "IMMUNE";
    default:
      return "UNKNOWN";
  }
}

/**
 * @brief Print a human-readable representation of an individual
 *
 * @param ind individual
 */
void print_individual(individual_t *ind) {
  printf(
      "--------------------\n     Individual\n--------------------\n id "
      "%lu\n pos   %.3f, %.3f\n displ %.3f, %.3f\n status %s\n t_status %lu\n",
      ind->id, ind->pos[0], ind->pos[1], ind->displ[0], ind->displ[1],
      individual_status_string(ind->status), ind->t_status);
}

/**
 * @brief Creates an individual
 *
 * The individual will have the specified id, <tt>pos, displ = {0,0}</tt>
 * <tt>status = NOT_EXPOSED</tt>, <tt>t_status = 0</tt>.
 * The list pointer is uninitialized.
 *
 * @param id
 * @return individual_t
 */
individual_t *create_individual(unsigned long id) {
  individual_t *ind = malloc(sizeof(individual_t));
  ind->id = id;
  ind->pos[0] = ind->pos[1] = ind->displ[0] = ind->displ[1] = 0.;
  ind->status = NOT_EXPOSED;
  ind->t_status = 0;
  return ind;
}

/**
 * @brief Create an empty list of individuals
 *
 * @return individual_list_t
 */
individual_list_t create_individual_list() {
  individual_list_t head = SLIST_HEAD_INITIALIZER(head);
  SLIST_INIT(&head);
  return head;
}