/**
    similarity.cpp
    (Dis)similarity measures for mass spectra
    
    @author Antoine Passemiers
    @version 1.0 30/03/2019
*/

#include "similarity.hpp"

namespace wassersteinms {

double wassersteinDistance(Spectrum &S1, Spectrum &S2) {
    Spectrum N(S1), M(S2);
    auto i = N.begin(), j = M.begin();
    double distance = 0.0;
    while ((i != N.end()) && (j != M.end())) {
        intensity_t d = std::min(i->second, j->second);
        distance = distance + d * std::abs(i->first - j->first);
        i->second = i->second - d;
        j->second = j->second - d;
        if (i->second == 0.0) {
            ++i;
        } else {
            ++j;
        }
    }
    return distance;
}


double euclideanDistance(Spectrum &S1, Spectrum &S2, double resolution) {
    // TODO: normalize with (x - x_min) / (x_max - x_min)
    Spectrum N = S1.normalize().changeResolution(resolution);
    Spectrum M = S2.normalize().changeResolution(resolution);
    auto i = N.begin(), j = M.begin();
    double distance = 0.0;
    while ((i != N.end()) && (j != M.end())) {
        if (i->first == j->first) {
            distance += std::pow(i->second - j->second, 2.0);
            i++, j++;
        } else if (i->first < j->first) {
            i++;
        } else {
            j++;
        }
    }
    return std::sqrt(distance);
}


double jaccardScore(Spectrum &S1, Spectrum &S2, double resolution, double eps) {
    // TODO: normalize with (x - x_min) / (x_max - x_min)
    Spectrum N = S1.changeResolution(resolution);
    Spectrum M = S2.changeResolution(resolution);
    auto i = N.begin(), j = M.begin();
    unsigned int nPeaks = 0, nMatchingPeaks = 0;
    while ((i != N.end()) && (j != M.end())) {
        if (i->first == j->first) {
            double min = std::min(i->second, j->second);
            double max = std::max(i->second, j->second);
            double delta = (max - min) / max;
            if (delta < eps) nMatchingPeaks += 2;
            nPeaks += 2;
            i++, j++;
        } else if (i->first < j->first) {
            nPeaks++;
            i++;
        } else {
            nPeaks++;
            j++;
        }
    }
    return static_cast<double>(nMatchingPeaks) / nPeaks;
}

} // namespace wassersteinms