#include "spectrum.hpp"


intensity_t& Spectrum::operator[](const mz_t key) {
    return this->peaks[key];
}

void Spectrum::normalize() {
    intensity_t s = 0.0;
    for (auto it = this->begin(); it != this->end(); it++) s += it->second;
    if (s != 0.0) {
        for (auto it = this->begin(); it != this->end(); it++) it->second /= s;
    }
}

void Spectrum::addKeys(Spectrum &other) {
    for (std::map<mz_t, intensity_t>::iterator it = other.begin(); it != other.end(); it++) {
        if (this->peaks.find(it->first) == this->peaks.end()) {
            this->peaks[it->first] = 0.0;
        }
    }
}