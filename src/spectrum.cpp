#include "spectrum.hpp"


void Spectrum::add(mz_t ratio, intensity_t intensity, intensity_t relIntensity) {
    sorted = false;
    ratios.push_back(ratio);
    intensities.push_back(intensity);
    relIntensities.push_back(relIntensity);
}