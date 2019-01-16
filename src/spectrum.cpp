#include "spectrum.hpp"


Spectrum::Spectrum(Spectrum &other) {
    sorted = other.sorted;
    ratios = other.ratios;
    intensities = other.intensities;
}


void Spectrum::add(mz_t ratio, intensity_t intensity) {
    sorted = false;
    ratios.push_back(ratio);
    intensities.push_back(intensity);
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
        for (int i = 0; i < ratios.size(); i++) {
            newRatios[i] = ratios[indices[i]];
            newIntensities[i] = intensities[indices[i]];
        }
        ratios = newRatios;
        intensities = newIntensities;
        sorted = true;
    }
}