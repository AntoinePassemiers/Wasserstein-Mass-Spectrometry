#ifndef SPECTRUM_HPP_
#define SPECTRUM_HPP_

#include <algorithm>
#include <functional>
#include <vector>
#include <numeric>
#include <iostream>


typedef double mz_t;
typedef double intensity_t;


class Spectrum {
private:
    bool sorted;
    std::vector<mz_t> ratios;
    std::vector<intensity_t> intensities;
public:
    explicit Spectrum() : sorted(true) {}
    explicit Spectrum(Spectrum &other);
    Spectrum(const Spectrum &other) = delete;
    Spectrum& operator=(const Spectrum&) = delete;
    ~Spectrum() = default;
    void addRatio(mz_t ratio);
    void add(mz_t ratio, intensity_t intensity);
    void sort();
    void normalize();
    size_t length() { return this->ratios.size(); }
    mz_t getRatio(size_t i) { return this->ratios[i]; }
    intensity_t getIntensity(size_t i) { return this->intensities[i]; }
    void setRatio(size_t i, mz_t value) { this->ratios[i] = value; }
    void setIntensity(size_t i, intensity_t value) { this->intensities[i] = value; }
};

#endif // SPECTRUM_HPP_