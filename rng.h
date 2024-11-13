#ifndef RNG_H
#define RNG_H

#include <random>

std::mt19937_64 rng;
std::uniform_real_distribution<double> unif;

/**
 * @brief ran - 0.0 to 1.0
 * @return
 */
double ran() {
    return unif(rng);
}

#endif // RNG_H
