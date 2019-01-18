#include "spectrum.hpp"


Spectrum::Spectrum(Spectrum &other) {
    this->sorted = other.sorted;
    this->ratios = other.ratios;
    this->intensities = other.intensities;
}

void Spectrum::addRatio(mz_t ratio) {
    this->sort();
    if (this->ratios.size() == 0) {
        this->ratios.push_back(ratio);
        this->intensities.push_back(0.0);
    } else {
        std::vector<mz_t>::iterator it;
        size_t i = 0;
        for (it = this->ratios.begin(); it != this->ratios.end(); it++, i++) {
            if (ratio == *it) {
                break;
            } else if (ratio > *it) {
                // Do nothing
            } else {
                this->ratios.insert(this->ratios.begin()+i, ratio);
                this->intensities.insert(this->intensities.begin()+i, 0.0);
                break;
            }
        }
        if (i == this->ratios.size()) {
            this->ratios.push_back(ratio);
            this->intensities.push_back(0.0);
        }
    }
}

void Spectrum::add(mz_t ratio, intensity_t intensity) {
    this->sorted = false;
    this->ratios.push_back(ratio);
    this->intensities.push_back(intensity);
}

void Spectrum:: normalize() {
    intensity_t s = std::accumulate(this->intensities.begin(), this->intensities.end(), 0);
    std::transform(
            this->intensities.begin(),
            this->intensities.end(),
            this->intensities.begin(),
            std::bind(std::divides<intensity_t>(), std::placeholders::_1, s));
}

void Spectrum::sort() {
    if (!sorted) {
        std::vector<int> indices(this->ratios.size());
        std::size_t n = 0;
        std::generate(std::begin(indices), std::end(indices), [&]{ return n++; });
        std::sort(
                std::begin(indices),
                std::end(indices),
                [&](int idx1, int idx2) { return this->ratios[idx1] < this->ratios[idx2]; } );
        std::vector<mz_t> newRatios(this->ratios.size());
        std::vector<intensity_t> newIntensities(this->ratios.size());
        for (size_t i = 0; i < this->ratios.size(); i++) {
            newRatios[i] = this->ratios[indices[i]];
            newIntensities[i] = this->intensities[indices[i]];
        }
        this->ratios = newRatios;
        this->intensities = newIntensities;
        this->sorted = true;
    }
}