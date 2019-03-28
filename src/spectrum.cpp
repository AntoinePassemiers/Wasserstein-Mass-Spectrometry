#include "spectrum.hpp"


intensity_t& Spectrum::operator[](const mz_t key) {
    return this->peaks[key];
}


Spectrum Spectrum::changeResolution(const double resolution) {
    Spectrum newSpectrum(*this);
    for (auto it = newSpectrum.begin(); it != newSpectrum.end(); it++) {
        mz_t key = it->first;
        intensity_t value = it->second;
        newSpectrum.remove(key);
        mz_t roundedKey = resolution * std::round(key / resolution);
        if (this->peaks.find(roundedKey) != this->peaks.end()) {
            newSpectrum[roundedKey] += value;
        } else {
            newSpectrum[roundedKey] = value;
        }
    }
    return newSpectrum;
}


Spectrum Spectrum::normalize() {
    Spectrum newSpectrum(*this);
    intensity_t s = 0.0;
    for (auto it = newSpectrum.begin(); it != newSpectrum.end(); it++) s += it->second;
    if (s != 0.0) {
        for (auto it = newSpectrum.begin(); it != newSpectrum.end(); it++) it->second /= s;
    }
    return newSpectrum;
}


void Spectrum::addKeys(Spectrum &other) {
    for (std::map<mz_t, intensity_t>::iterator it = other.begin(); it != other.end(); it++) {
        if (this->peaks.find(it->first) == this->peaks.end()) {
            this->peaks[it->first] = 0.0;
        }
    }
}