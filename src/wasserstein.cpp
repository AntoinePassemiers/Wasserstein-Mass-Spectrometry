#include "wasserstein.hpp"


double wassersteinDistance(std::unique_ptr<Spectrum> &S1, std::unique_ptr<Spectrum> &S2) {
    Spectrum N(*S1), M(*S2);
    N.sort();
    M.sort();
    size_t i = 0, j = 0;
    double distance = 0.0;
    size_t n = N.length(), m = M.length();
    while ((i < n) && (j < m)) {
        intensity_t d = std::min(N.getIntensity(i), M.getIntensity(j));
        distance = distance + d * std::abs(N.getRatio(i) - M.getRatio(j));
        N.setIntensity(i, N.getIntensity(i) - d);
        M.setIntensity(j, M.getIntensity(j) - d);
        if (N.getIntensity(i) == 0.0) {
            i++;
        } else {
            j++;
        }
    }
    return distance;
}