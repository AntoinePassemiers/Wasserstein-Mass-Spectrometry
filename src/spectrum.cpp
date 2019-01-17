#include "spectrum.hpp"


Spectrum::Spectrum(Spectrum &other) {
    sorted = other.sorted;
    ratios = other.ratios;
    intensities = other.intensities;
}

void Spectrum::add(mz_t ratio) {
    sort();
    std::vector<mz_t>::iterator it;
    size_t i = 0;
    for (it = ratios.begin(); it != ratios.end(); it++, i++) {
        if (ratio < *it) {
            // Do nothing
        } else if (ratio == *it) {
            break;
        } else {
            ratios.insert(ratios.begin()+i, ratio);
            intensities.insert(ratios.begin()+i, 0.0);
            break;
        }
    }
}

void Spectrum::add(mz_t ratio, intensity_t intensity) {
    sorted = false;
    ratios.push_back(ratio);
    intensities.push_back(intensity);
}

void Spectrum:: normalize() {
    intensity_t s = std::accumulate(intensities.begin(), intensities.end(), 0);
    std::transform(
            intensities.begin(),
            intensities.end(),
            intensities.begin(),
            std::bind(std::divides<intensity_t>(), std::placeholders::_1, s));
}

void Spectrum::sort() {
    if (!sorted) {
        std::vector<int> indices(ratios.size());
        std::size_t n = 0;
        std::generate(std::begin(indices), std::end(indices), [&]{ return n++; });
        std::sort(
                std::begin(indices),
                std::end(indices),
                [&](int idx1, int idx2) { return ratios[idx1] < ratios[idx2]; } );
        std::vector<mz_t> newRatios(ratios.size());
        std::vector<intensity_t> newIntensities(ratios.size());
        for (size_t i = 0; i < ratios.size(); i++) {
            newRatios[i] = ratios[indices[i]];
            newIntensities[i] = intensities[indices[i]];
        }
        ratios = newRatios;
        intensities = newIntensities;
        sorted = true;
    }
}