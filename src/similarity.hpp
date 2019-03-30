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

    // First Wasserstein distance
    WASSERSTEIN,

    // Euclidean distance
    EUCLIDEAN,

    // Jaccard score
    JACCARD_SCORE
};

/**
    Computes Wasserstein distance between two mass spectra.

    @param S1 First spectrum
    @param S2 Second spectrum
    @return First Wasserstein distance
*/
double wassersteinDistance(Spectrum &S1, Spectrum &S2);

/**
    Computes Euclidean distance between two mass spectra.

    @param S1 First spectrum
    @param S2 Second spectrum
    @param resolution Resolution to use for identifying matching peaks
    @return Euclidean
*/
double euclideanDistance(Spectrum &S1, Spectrum &S2, double resolution);

/**
    Computes Jaccard score between two mass spectra.

    @param S1 First spectrum
    @param S2 Second spectrum
    @param resolution Resolution to use for identifying matching peaks
    @param eps Relative intensity threshold below which peaks are matched
    @return First Wasserstein distance
*/
double jaccardScore(Spectrum &S1, Spectrum &S2, double resolution, double eps);

} // namespace wassersteinms

#endif // SIMILARITY_HPP__