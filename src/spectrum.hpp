/**
    spectrum.hpp
    Mass spectra
    
    @author Antoine Passemiers
    @version 1.0 30/03/2019
*/

#ifndef SPECTRUM_HPP_
#define SPECTRUM_HPP_

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <vector>
#include <numeric>
#include <iostream>

namespace wassersteinms {

typedef double mz_t;
typedef double intensity_t;
const double MAX_INTENSITY = std::numeric_limits<double>::infinity();


class Spectrum {
    friend std::ostream& operator<<(std::ostream&, const Spectrum&);
private:
    std::map<mz_t, intensity_t> peaks;

public:
    explicit Spectrum(): peaks() {}
    explicit Spectrum(Spectrum &other): peaks(other.peaks) {}
    Spectrum(const Spectrum &other): peaks(other.peaks) {}
    Spectrum& operator=(const Spectrum& other) { this->peaks = other.peaks; return *this; }
    ~Spectrum() = default;

    intensity_t& operator[](const mz_t key) { return this->peaks[key]; }
    const intensity_t& operator[](const mz_t key) const { return this->peaks.at(key); }

    Spectrum normalize();
    Spectrum scale();
    Spectrum changeResolution(const double resolution);
    void remove(mz_t key) { this->peaks.erase(key); }
    void addKeys(Spectrum &other);
    std::map<mz_t, intensity_t>::iterator begin() { return this->peaks.begin(); }
    std::map<mz_t, intensity_t>::iterator end() { return this->peaks.end(); }
    std::map<mz_t, intensity_t>::const_iterator begin() const { return this->peaks.begin(); }
    std::map<mz_t, intensity_t>::const_iterator end() const { return this->peaks.end(); }
    const int size() const { return this->peaks.size(); }
};

} // namespace wassersteinms

#endif // SPECTRUM_HPP_