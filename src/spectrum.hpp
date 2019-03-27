#ifndef SPECTRUM_HPP_
#define SPECTRUM_HPP_

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <vector>
#include <numeric>
#include <iostream>


typedef double mz_t;
typedef double intensity_t;


class Spectrum {
private:
    std::map<mz_t, intensity_t> peaks;

public:
    explicit Spectrum(): peaks() {}
    explicit Spectrum(Spectrum &other): peaks(other.peaks) {}
    Spectrum(const Spectrum &other): peaks(other.peaks) {}
    Spectrum& operator=(const Spectrum& other) { this->peaks = other.peaks; return *this; }
    ~Spectrum() = default;

    intensity_t& operator[](const mz_t key);

    Spectrum normalize();
    Spectrum changeResolution(const double resolution);
    void remove(mz_t key) { this->peaks.erase(key); }
    void addKeys(Spectrum &other);
    std::map<mz_t, intensity_t>::iterator begin() { return this->peaks.begin(); }
    std::map<mz_t, intensity_t>::iterator end() { return this->peaks.end(); }
    int size() { return this->peaks.size(); }
};

#endif // SPECTRUM_HPP_