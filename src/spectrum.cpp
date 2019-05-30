/**
    spectrum.cpp
    Mass spectra
    
    @author Antoine Passemiers
    @version 1.0 30/03/2019
*/

#include "spectrum.hpp"

namespace wassersteinms {

Spectrum Spectrum::changeResolution(const double resolution) {
    Spectrum newSpectrum;
    for (auto it = this->begin(); it != this->end(); it++) {
        mz_t key = it->first;
        intensity_t value = it->second;
        mz_t roundedKey = resolution * std::round(key / resolution);
        if (newSpectrum.peaks.find(roundedKey) != newSpectrum.end()) {
            it->second += value;
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


Spectrum Spectrum::scale() {
    Spectrum newSpectrum(*this);
    intensity_t max = -MAX_INTENSITY;
    intensity_t min = MAX_INTENSITY;
    for (auto it = newSpectrum.begin(); it != newSpectrum.end(); it++) {
        if (it->second < min) {
            min = it->second;
        } else if (it->second > max) {
            max = it->second;
        }
    }
    for (auto it = newSpectrum.begin(); it != newSpectrum.end(); it++) {
        it->second = (it->second - min) / (max - min);
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


std::ostream& operator<<(std::ostream& out, const Spectrum& arr) {
    out << "[ ";
    if (arr.size() < 1000) {
        for (auto it = arr.begin(); it != arr.end(); ++it) {
            out << ' ' << it->second << ' ';
        }
    } else {
        auto it = arr.begin();
        for (int i = 0; i < 15; ++i) {
            out << ' ' << it->second << ' ';
            it++;
        }
        out << "...";
        for (auto it = std::prev(arr.end(), 5); it != arr.end(); it++) {
            out << ' ' << it->second << ' ';
        }
    }
    out << " ] " << std::endl;
    return out;
}

} // namespace wassersteinms
