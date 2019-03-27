#include "similarity.hpp"


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