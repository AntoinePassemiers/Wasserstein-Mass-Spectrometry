/**
    similarity.hpp
    (Dis)similarity measures for mass spectra
    
    @author Antoine Passemiers
    @version 1.0 30/03/2019
*/

#ifndef SIMILARITY_HPP__
#define SIMILARITY_HPP__

#include "spectrum.hpp"

#include <algorithm>
#include <cmath>
#include <memory>

namespace wassersteinms {

enum class Similarities {
    WASSERSTEIN,
    EUCLIDEAN,
    JACCARD_SCORE
};


double wassersteinDistance(Spectrum &S1, Spectrum &S2);


double euclideanDistance(Spectrum &S1, Spectrum &S2, double resolution);


double jaccardScore(Spectrum &S1, Spectrum &S2, double resolution, double eps);

} // namespace wassersteinms

#endif // SIMILARITY_HPP__