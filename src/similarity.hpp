#ifndef SIMILARITY_HPP__
#define SIMILARITY_HPP__

#include "spectrum.hpp"

#include <algorithm>
#include <cmath>
#include <memory>


enum class Similarities {
    WASSERSTEIN,
    EUCLIDEAN,
    JACCARD_SCORE
};


double wassersteinDistance(Spectrum &S1, Spectrum &S2);


double euclideanDistance(Spectrum &S1, Spectrum &S2, double resolution);


double jaccardScore(Spectrum &S1, Spectrum &S2, double resolution);


#endif // SIMILARITY_HPP__