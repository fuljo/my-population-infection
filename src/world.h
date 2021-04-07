
/**
 * @brief Uniformly distributes a population between countries.
 *
 * All countries will get <tt>floor(population / num_countries)</tt>
 * individuals, and the last country will also get the remaining <tt>population
 * % num_countries</tt> individuals
 *
 * @param population number of individuals to be distributed
 * @param num_countries number of countries
 * @param res array of size \p num_countries that will hold the result
 */
void distribute_population_uniform(unsigned long population,
                                   unsigned int num_countries,
                                   unsigned long res[]);